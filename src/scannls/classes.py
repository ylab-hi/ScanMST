#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
import argparse
import os
import random
import re
import subprocess
import time
from multiprocessing import Process
from typing import Any
from typing import List
from typing import Tuple
from typing import Union

import psutil
from align import aligner
from Bio import SearchIO
from Bio.Seq import Seq
from loguru import logger

from .draft.helper import cigar_validity
from .draft.nls_inference import infer_nls_from_connected_reads
from .exception import ReadNotFoundError
from .utils import reverse_complement


class Path(object):
    """store chimeric reads as nodes in a path (directed acyclic graph)
    :param nodes: a list of Read as nodes
    :type nodes: List[Read]
    :param sms: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length)
    :type sms: tuple or None
    :param sequence: reads sequence of last added read
    :type sequence: str or None
    :param nm: summation of number-of-mismatches of Reads in the path
    :type nm: int
    :param mode: a dictionary of Reads-pair to mode in the path
    :type mode: dict
    .. note::
        We have to pay attention on 'sms':
            * every path only keep one 'sms' value (per path instead of per node).
            * every path only keep one 'sequence' value (per path).
            * once new node added to the path, update the value of 'sms' using the summed 'sms' value from the function 'test_is_connected'
            * once new node added to the path, update the value of 'sequence' using the summed 'sequence' value from the function 'test_is_connected'
    """

    __slots__ = ("nodes", "sms", "sequence", "nm", "mode")

    def __init__(self) -> None:
        self.nodes = []
        self.sms = None
        self.sequence = None
        self.nm = 0
        self.mode = {}

    def add(self, read) -> None:
        """add one chimeric read to the path"""
        self.nodes.append(read)
        self.nm = self.nm + read.nm

    def add_mode(self, read_pair_dict) -> None:
        """add read-pair=> mode to the path"""
        self.mode.update(read_pair_dict)

    def __len__(self) -> int:
        return len(self.nodes)

    def __lt__(self, other) -> bool:
        return len(self.nodes) < len(other.nodes)

    def __repr__(self) -> str:
        return ";".join(map(str, self.nodes))

    def __hash__(self) -> int:
        return hash(";".join(map(str, self.nodes)))

    def __eq__(self, other) -> bool:
        return ";".join(map(str, self.nodes)) == ";".join(map(str, other.nodes))


class Read(object):
    """build a read class for storing information of every junction read
    :param chrom: chromosome of genome
    :type chrom: str
    :param ref_start: start position of chimeric read
    :type ref_start: int
    :param strand: direction of chimeric read (-|+)
    :type strand: str
    :param cigarstring: cigar string of chimeric read (-|+)
    :type cigarstring: str
    :param mapq: MAPQ of chimeric read
    :type mapq: int
    :param nm: number of mismatches of chimeric read
    :type nm: int
    :param query_sequence: read sequence in the BAM file
    :type query_sequence: str
    :param linked_paths: linked paths for the read
    :type linked_paths: list
    :param lt_soft_len: softclipped segment length on the left side
    :type lt_soft_len: int
    :param rt_soft_len: softclipped segment length on the right side
    :type rt_soft_len: int
    :param read_match_size: M+I
    :type read_match_size: int
    :param reference_match_size: M+D+N
    :type reference_match_size: int
    :param indel_size: D+N-I
    :type indel_size: int
    :param cigartuples: cigarstring tuple version: [ (operation code, length) ]; operation code: {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    :type cigartuples: list
    :param cigartuples_without_soft: cigarstring tuple verion [exclude softclipping]: [(operation code, length)]; operation code: {'M':0,'I':1,'D':2,'N':3}
    :type cigartuples_without_soft: list
    :param query_length: length of the chimeric read
    :type query_length: int
    """

    __slots__ = (
        "chrom",
        "ref_start",
        "strand",
        "cigarstring",
        "mapq",
        "nm",
        "query_sequence",
        "linked_paths",
        "lt_soft_len",
        "rt_soft_len",
        "read_match_size",
        "reference_match_size",
        "indel_size",
        "cigartuples_without_soft",
        "cigartuples",
        "query_length",
        "adhocsms",
        "adhocseq",
        "mode",
    )

    def __init__(
        self,
        chrom,
        ref_start,
        strand,
        cigarstring,
        mapq,
        nm,
        query_sequence,
        lt_soft_len,
        rt_soft_len,
        read_match_size,
        reference_match_size,
        indel_size,
        cigartuples_without_soft,
        query_length,
        cigartuples,
    ) -> None:
        self.chrom = chrom
        self.ref_start = ref_start
        self.strand = strand
        self.cigarstring = cigarstring
        self.mapq = mapq
        self.nm = nm
        self.query_sequence = query_sequence
        self.linked_paths = []
        self.lt_soft_len = lt_soft_len
        self.rt_soft_len = rt_soft_len
        self.read_match_size = read_match_size
        self.reference_match_size = reference_match_size
        self.indel_size = indel_size
        self.cigartuples_without_soft = cigartuples_without_soft
        self.query_length = query_length
        self.cigartuples = cigartuples

        self.adhocsms = None
        self.adhocseq = None
        self.mode = None

    def __eq__(self, other) -> bool:
        if isinstance(other, Read):
            if (
                self.chrom == other.chrom
                and self.ref_start == other.ref_start
                and self.ref_end == other.ref_end
                and self.strand == other.strand
                and self.mapq == other.mapq
                and self.nm == other.nm
            ):
                return True
        return False

    def __hash__(self) -> int:
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.strand)
            ^ hash(self.mapq)
            ^ hash(self.nm)
        )

    def __lt__(self, other) -> bool:
        if self.chrom == other.chorm:
            return self.ref_start < other.ref_start
        else:
            chrm_dict = {"chrM": 0, "MT": 0, "chrX": 23, "chrY": 24, "X": 23, "Y": 24}
            for i in range(1, 23):
                chrm_dict.update({f"chr{i}": i})
                chrm_dict.update({f"{i}": i})
            return chrm_dict[self.chrom] < chrm_dict[other.chrom]

    # for debug purpose
    def __repr__(self) -> str:
        return fr"Read({self.chrom}, {self.ref_start}, {self.ref_end}, {self.strand}, {self.mapq}, {self.nm})"

    def __str__(self) -> str:
        return ",".join(
            map(
                str,
                [
                    self.chrom,
                    self.ref_start,
                    self.ref_end,
                    self.strand,
                    self.mapq,
                    self.nm,
                    f"{self.lt_soft_len}:{self.read_match_size}:{self.rt_soft_len}",
                ],
            ),
        )

    @staticmethod
    def _calculate_features(cigar_str):
        cigar_char_dict = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5}
        # 'length', 'operation char'
        len_type_tuple = re.findall(r"(\d+)(\w)", cigar_str)
        # (operation code, length)
        cigartuples = [(cigar_char_dict[j], int(i)) for i, j in len_type_tuple]

        query_length = 0
        indel_size = 0
        reference_match_size = 0
        read_match_size = 0
        cigar_without_soft = []
        for op_code, _len_ in cigartuples:
            if op_code == 0:  # M
                reference_match_size += _len_
                read_match_size += _len_
                query_length += _len_
                cigar_without_soft.append([0, _len_])
            elif op_code == 1:  # I
                indel_size += -_len_
                read_match_size += _len_
                query_length += _len_
                cigar_without_soft.append([1, _len_])
            elif op_code == 2:  # D
                indel_size += _len_
                reference_match_size += _len_
                cigar_without_soft.append([2, _len_])
            elif op_code == 3:  # N
                indel_size += _len_
                reference_match_size += _len_
                cigar_without_soft.append([3, _len_])
            elif op_code == 4:  # S
                query_length += _len_

        lt_soft_len = 0
        rt_soft_len = 0
        lt_op, lt_len = cigartuples[0]
        rt_op, rt_len = cigartuples[-1]
        if lt_op == 4:
            lt_soft_len = lt_len
        if rt_op == 4:
            rt_soft_len = rt_len
        return (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigar_without_soft,
            query_length,
            cigartuples,
        )

    @classmethod
    def init(cls, chrom, position, strand, cigar_str, mapq, nm, query_seq):

        (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigar_without_soft,
            query_length,
            cigartuples,
        ) = Read._calculate_features(cigar_str)

        return cls(
            chrom,
            position,
            strand,
            cigar_str,
            mapq,
            nm,
            query_seq,
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigar_without_soft,
            query_length,
            cigartuples,
        )

    @property
    def ref_end(self) -> int:
        return self.ref_start + self.reference_match_size

    @property
    def reference_span(self) -> int:
        """M+N+D"""
        return self.reference_match_size

    @property
    def sms(self) -> tuple:
        return self.lt_soft_len, self.read_match_size, self.rt_soft_len

    def add_path(self, path) -> None:
        """path is an instance of Path class"""
        self.linked_paths.append(path)

    def get_exons_and_introns(self) -> tuple:
        """get the coordinates for reads matched part (without softclipping)
        :return: exons coordinates and introns coordinates
        :rtype: tuple
        """
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start
        for op_code, _len_ in self.cigartuples_without_soft:
            if op_code in {0, 2}:  # M, D
                current_pos = current_pos + _len_
            elif op_code == 3:  # N
                exons.append((start_pos, current_pos))
                current_pos = current_pos + _len_
                start_pos = current_pos
        exons.append((start_pos, current_pos))

        # No 'N' in the cigar
        if len(exons) == 1:
            introns = []
        elif len(exons) > 1:
            _positions = []
            for i, j in exons:
                _positions.extend([i, j])
            _positions.sort()
            _positions.pop(0)
            _positions.pop(-1)
            introns = list(zip(_positions[::2], _positions[1::2]))
        return exons, introns

    def splice_site_checker(self, genome_fasta, fraction_cutoff=0.6) -> bool:
        """check whether the fraction of canonical splice site usage in read reference matched part is bigger than 'fraction_cutoff' or not
        :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
        :param fraction_cutoff: fraction of canonical splice sites used in the putative introns inferred from the CIGAR
        :type genome_fasta: pyfaidx.Fasta
        :type fraction_cutoff: float
        :return: using canonical splice sites OR not
        :rtype: bool
        """
        exons, introns = self.get_exons_and_introns()

        if len(introns) == 0:
            return True

        intron_count = 0
        can_count = 0
        can_sites = {"GT-AG", "GC-AG", "AT-AC"}
        for start, end in introns:
            if end - start >= 10:
                intron_count += 1
                if self.strand == "-":
                    left_site = genome_fasta[self.chrom][
                        end - 2 : end
                    ].reverse.complement.seq
                    right_site = genome_fasta[self.chrom][
                        start : start + 2
                    ].reverse.complement.seq
                else:
                    left_site = genome_fasta[self.chrom][start : start + 2].seq
                    right_site = genome_fasta[self.chrom][end - 2 : end].seq
                if f"{left_site}-{right_site}" in can_sites:
                    can_count += 1
        if can_count / intron_count >= fraction_cutoff:
            return True
        else:
            return False


class NoneInsertion(Read):
    def __init__(self, hit_num: int, query_sequence: str):
        self.query_sequence = query_sequence
        self.hit = hit_num

    # self.chrom = chrom
    # self.prev_breakpoint = prev_bp
    # self.next_breakpoint = next_bp
    # self.strand = strand
    # self.ref_start = ref_start
    # self.ref_end = ref_end
    # self.exons = exons
    # self.sv_type = sv_type
    # self.modes = modes
    # self.genes = genes
    # self.annotation_code = annot
    # self.splicing_code = canonical
    # self.sr = sr
    # self.insertion_info = insertion_info

    # self.prev_breakpoint = prev_bp
    # self.next_breakpoint = next_bp
    # self.exons = exons
    # self.sv_type = sv_type
    # self.modes = modes
    # self.genes = genes
    # self.annotation_code = annot
    # self.splicing_code = canonical
    # self.sr = sr
    # self.insertion_info = insertion_info
    #
    #     self.chrom = chrom
    #     self.ref_start = position
    #     self.strand = strand
    #     self.cigarstring = cigar_str
    #     self.mapq = mapq
    #     self.nm = nm
    #     self.query_sequence = query_seq
    #     self.linked_paths = []
    #     self.lt_soft_len = lt_soft_len
    #     self.rt_soft_len = rt_soft_len
    #     self.read_match_size = read_match_size
    #     self.reference_match_size = reference_match_size
    #     self.indel_size = indel_size
    #     self.cigartuples_without_soft = cigar_without_soft
    #     self.query_length = query_length
    #     self.cigartuples = cigartuples
    #
    #     self.adhocsms = None
    #     self.adhocseq = None
    #     self.mode = None


class Insertion(Read):
    def __init__(
        self,
        hit_num: int,
        chrom: str,
        ref_start: int,
        strand: str,
        cigarstring: str,
        mapq: int,
        nm: int,
        query_sequence: str,
    ):
        (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigar_without_soft,
            query_length,
            cigartuples,
        ) = Read._calculate_features(cigarstring)
        super().__init__(
            chrom,
            ref_start,
            strand,
            cigarstring,
            mapq,
            nm,
            query_sequence,
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigar_without_soft,
            query_length,
            cigartuples,
        )
        self.hit_num = hit_num
        self.sv_type = None

        self.prev_breakpoint = None
        self.next_breakpoint = None
        self.exons = None
        self.modes = None
        self.genes = None
        self.annotation_code = None
        self.splicing_code = None
        self.sr = None

    def update_cigarstring(self, sms, source_s):
        _ls, _m, _rs = sms
        if source_s == "left":
            ls = _ls - self.query_length
            rs = _rs + _m
        else:
            ls = _ls + _m
            rs = _rs - self.query_length

        self.cigarstring = cigar_validity(f"{ls}S{self.cigarstring}{rs}S")

    def reverse_completement_query(self):
        self.query_sequence = reverse_complement(self.query_sequence)


class LengthAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values <= 0:
            parser.error("Minimum length for {0} is 1".format(option_string))
        setattr(namespace, self.dest, values)


class Node(object):
    """build a breakpoint node class for storing information of every breakpoint
    :param prev_bp: breakpoint for the previous breakpoints connections
    :type prev_bp: str
    :param next_bp: breakpoint for the next breakpoints connections
    :type next_bp: str
    :param strand: direction of chimeric read (-|+)
    :type strand: str
    :param chrom: chromosome
    :type chrom: str
    :param ref_start: reference start position
    :type ref_start: int
    :param ref_end: reference end position
    :type ref_end: int
    :param exons: CIGAR inferred exons in the read. e.g., [(100, 200), (300, 500)]
    :type exons: list of tuple
    :param sv_type: one of the SV types (TDUP/INV/TRA)
    :type sv_type: str
    :param annot: gene annotation code
    :type annot: int
    :param canonical: canonical splice site code {1: canonical, 0: noncanonical}
    :type canonical: int
    :param modes: read modes of connected breakpoints
    :type modes: tuple
    :param genes: overlapped genes of connected breakpoints
    :type genes: tuple

    .. note::
        connection-level fields:
        * sv_type
        * modes
        * genes
        * annotation_code
        * splicing_code
    """

    __slots__ = (
        "next_breakpoint",
        "prev_breakpoint",
        "strand",
        "chrom",
        "ref_start",
        "ref_end",
        "exons",
        "sv_type",
        "modes",
        "genes",
        "annotation_code",
        "splicing_code",
        "sr",
        "insertion_info",
    )

    def __init__(
        self,
        prev_bp=None,
        next_bp=None,
        strand=None,
        chrom=None,
        ref_start=None,
        ref_end=None,
        exons=None,
        sv_type=None,
        annot=None,
        canonical=None,
        modes=None,
        genes=None,
        sr=None,
        insertion_info=None,
    ) -> None:
        self.chrom = chrom
        self.prev_breakpoint = prev_bp
        self.next_breakpoint = next_bp
        self.strand = strand
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.exons = exons
        self.sv_type = sv_type
        self.modes = modes
        self.genes = genes
        self.annotation_code = annot
        self.splicing_code = canonical
        self.sr = sr
        self.insertion_info = insertion_info

    def __eq__(self, other) -> bool:
        if isinstance(other, Node):
            if (
                self.chrom == other.chrom
                and self.ref_start == other.ref_start
                and self.ref_end == other.ref_end
                and self.sv_type == other.sv_type
                and self.prev_breakpoint == other.prev_breakpoint
                and self.next_breakpoint == other.next_breakpoint
                and self.strand == other.strand
            ):
                return True
        return False

    def __hash__(self) -> int:
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.sv_type)
            ^ hash(self.prev_breakpoint)
            ^ hash(self.next_breakpoint)
            ^ hash(self.strand)
        )

    # for debug purpose
    def __repr__(self) -> str:
        exons_repr = "|".join([f"{i}-{j}" for i, j in self.exons])
        return fr"Node({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}, {exons_repr}, {self.sv_type}, {self.prev_breakpoint}, {self.next_breakpoint})"

    def __str__(self) -> str:
        exons_repr = "|".join([f"{i}-{j}" for i, j in self.exons])
        return fr"Node({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}, {exons_repr}, {self.sv_type}, {self.prev_breakpoint}, {self.next_breakpoint})"

    def is_next_node(self, other) -> bool:
        if self.next_breakpoint == other.prev_breakpoint:
            return True
        else:
            return False

    def is_previous_node(self, other) -> bool:
        if self.prev_breakpoint == other.next_breakpoint:
            return True
        else:
            return False

    @classmethod
    def create_nodes(cls, number):
        return [cls() for _ in range(number)]

    @property
    def introns(self):
        if len(self.exons) <= 1:
            return []
        else:
            _positions = []
            for i, j in self.exons:
                _positions.extend([i, j])
            _positions.pop(0)
            _positions.pop(-1)
            _introns = list(zip(_positions[::2], _positions[1::2]))
            return _introns


class Event:
    def __init__(self, event):
        (
            sv_type,
            annot,
            canonical,
            _positions,
            read1_info,
            read2_info,
            insertion_info,
            strands,
            genes,
        ) = event

        self.sv_type = sv_type
        if self.sv_type != "NA":
            self.annotation_code = annot
            self.splicing_code = canonical
            self.genes = genes
            self.insertion_info = insertion_info
            self.bp1, self.bp2 = _positions[:2]
            self.mode1, self.mode2 = _positions[2:]
            self.strand1, self.strand2 = strands
            self.read1_ref_start, self.read1_ref_end, self.read1_exons = read1_info
            self.read2_ref_start, self.read2_ref_end, self.read2_exons = read2_info

    @property
    def modes(self):
        return [self.mode1, self.mode2]

    @property
    def chrom1(self):
        return self.bp1.split(":")[0]

    @property
    def chrom2(self):
        return self.bp2.split(":")[0]

    @property
    def insertion_seq1(self):
        return self.insertion_info[0][1:]

    @property
    def insertion_seq2(self):
        return self.insertion_info[1][1:]

    def is_NA(self):
        return True if self.sv_type == "NA" else False

    def has_insertion(self):
        return True if self.insertion_info[0].startswith("+") else False

    def is_same_strand(self):
        return self.strand1 == self.strand2

    def read1(self, read_chains):
        for read in read_chains:
            if read.ref_start == self.read1_ref_start:
                return read
        else:
            raise ReadNotFoundError

    def read2(self, read_chains):
        for read in read_chains:
            if read.ref_start == self.read2_ref_start:
                return read
        else:
            raise ReadNotFoundError

    def update_specific_info_within_event(self, node, info_key_list):
        """
        update node info from the event by the info_key_list

        :param node:  Node
        :param info_key_list: [key1, key2, ...]
        :return: Node with updated info
        """

        for key in info_key_list:
            setattr(node, key, getattr(self, key))

        return node

    def update_node_info(
        self, flag, new_node, insertion, is_update_insertion_info=True
    ):
        """
        update the common info the node in the front, and the common info includes

        sv_type, annot, canonical, genes, insertion_info, and the breakpoints, mode

        :param flag:
        :param new_node:
        :param insertion:
        :param is_update_insertion_info:
        :return:
        """
        new_node = self.update_specific_info_within_event(
            new_node, ["sv_type", "annotation", "splicing_code", "modes", "genes"]
        )
        if is_update_insertion_info:
            new_node.insertion_info = (flag, insertion)
        return new_node

    def update_insertion_info(self, insertion):
        insertion.prev_breakpoint = insertion.ref_start
        insertion.next_breakpoint = insertion.ref_end

        insertion.exons, _ = insertion.get_exons_and_introns()
        # TODO: the event between insertion and read has these attributes?
        return self.update_specific_info_within_event(
            insertion, ["sv_type", "annotation", "splicing_code", "modes", "genes"]
        )


class Series(object):
    """construct a sequence of Nodes for storing information of connected breakpoints
    :param nodes: sequence of Nodes
    :type nodes: list
    :param assemblied: The series is from assembly of reads (True) or a single read (False)
    :type assemblied: bool or None

    .. note::
        [('TDUP', 0, 1, ('chr17:7708250', 'chr17:7701656', 1, 2), ('+', '+'), ['INTERGENIC', 'INTERGENIC']),
        ('TRA', 0, 1, ('chr17:7702552', 'chr1:15872815', 1, 2), ('+', '+'), ['INTERGENIC', 'INTERGENIC']),
        ('TDUP', 0, 1, ('chr1:15876678', 'chr1:15777169', 1, 2), ('+', '+'), ['INTERGENIC', 'INTERGENIC'])]

        Node(TDUP, None, chr17:7708250, +);Node(TRA, chr17:7701656, chr17:7702552, +);Node(TDUP, chr1:15872815, chr1:15876678, +);Node(None, chr1:15777169, None, +)


                         bp1               bp2  bp3                bp4
                ---------|------    -------|----|------    --------|---------
                       Node1                 Node2                Node3
    prev_breakpoint:   None                 bp2                   bp4
    next_breakpoint:    bp1                 bp3                   None
    sv_type:        TDUP/INV/TRA        TDUP/INV/TRA              None
    """

    def __init__(self, blat) -> None:
        self.nodes = []
        self.assemblied = None
        self.blat = blat

    def add_node(self, node: Node) -> None:
        self.nodes.append(node)

    def init(
        self,
        event_list,
        read_chains,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        update_bps=False,
    ) -> None:
        """add event list as Node to self.nodes"""
        event_list = [
            Event(event)
            for event in self.order_events_by_trancription_direction(event_list)
            if event[0] != "NA"
        ]
        event_list_len = len(event_list)
        previous_breakpoint = None
        for index, event in enumerate(event_list):

            read1_node = Node(
                prev_bp=previous_breakpoint,
                next_bp=event.bp1,
                strand=event.strand1,
                chrom=event.chrom1,
                ref_start=event.read1_ref_start,
                exons=event.read1_exons,
            )
            previous_breakpoint = event.bp2

            # is insertions
            if event.has_insertion():

                insertion_seq = event.insertion_seq1  # pick from the first read
                flag, insertion = self.blat.query_insertion(insertion_seq)
                if flag:  # only one hit
                    # add first node and insertion node

                    # get type of insertion between first node and insertion node
                    read1 = event.read1(read_chains)
                    insertion.update_cigarstring(read1.cigarstring, source_s="right")

                    insertion_mode = 2 if event.mode1 == 1 else 1

                    read1_insertion_event = Event(
                        infer_nls_from_connected_reads(
                            read_lt=read1,
                            read_rt=insertion,
                            lt_mode=event.mode1,
                            rt_mode=insertion_mode,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            update_bps=update_bps,
                        )
                    )

                    # get type of insertion between insertion node and second node
                    read2 = event.read2(read_chains)
                    insertion_mode = 2 if event.mode2 == 1 else 1

                    if event.strand1 != event.strand2:
                        insertion.reverse_completement_query()

                    insertion_read2_event = Event(
                        infer_nls_from_connected_reads(
                            read_lt=insertion,
                            read_rt=read2,
                            lt_mode=insertion_mode,
                            rt_mode=read2,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            update_bps=update_bps,
                        )
                    )
                    if read1_insertion_event.is_NA() or insertion_read2_event.is_NA():
                        # only add read1
                        read1_node = event.update_node_info(flag, read1_node, insertion)
                        self.add_node(read1_node)
                    else:
                        # add read1 and insertion
                        read1_node = read1_insertion_event.update_node_info(
                            flag, read1_node, insertion
                        )
                        self.add_node(read1_node)
                        insertion = insertion_read2_event.update_insertion_info(
                            insertion
                        )
                        self.add_node(insertion)

                else:  # no hits or multiple hits
                    read1_node = event.update_node_info(flag, read1_node, insertion)
                    self.add_node(read1_node)
            # no insertion
            else:
                read1_node = event.update_node_info(False, read1_node, None, False)
                self.add_node(read1_node)

            # add final node
            if index == event_list_len - 1:
                # TODO check breakpoint
                final_node = Node(
                    prev_bp=previous_breakpoint,
                    strand=event.strand2,
                    chrom=event.chrom2,
                    ref_start=event.read2_ref_start,
                    exons=event.read2_exons,
                )

                self.add_node(final_node)

    @staticmethod
    def reorder_event(event):
        """
        order breakpoints pairs following the transcription direction using information of reads 'mode' and 'strand'
        +1;-1 => up;down
        +2;-2 => down;up
        """
        (
            sv_type,
            annot,
            canonical,
            _positions,
            read1_info,
            read2_info,
            strands,
            genes,
        ) = event
        bp1 = _positions[0]
        bp2 = _positions[1]
        mode1 = _positions[2]
        mode2 = _positions[3]
        strand1 = strands[0]
        strand2 = strands[1]
        if strand1 == "+" and strand2 == "-":
            if mode1 == 1 and mode2 == 1:
                is_bp1_upstream = True
            elif mode1 == 2 and mode2 == 2:
                is_bp1_upstream = False
        elif strand1 == "-" and strand2 == "+":
            if mode1 == 1 and mode2 == 1:
                is_bp1_upstream = False
            elif mode1 == 2 and mode2 == 2:
                is_bp1_upstream = True
        elif strand1 == "+" and strand2 == "+":
            if mode1 == 1 and mode2 == 2:
                is_bp1_upstream = True
            elif mode1 == 2 and mode2 == 1:
                is_bp1_upstream = False
        elif strand1 == "-" and strand2 == "-":
            if mode1 == 1 and mode2 == 2:
                is_bp1_upstream = False
            elif mode1 == 2 and mode2 == 1:
                is_bp1_upstream = True

        if not is_bp1_upstream:
            if annot == 1:
                annot = 2
            elif annot == 2:
                annot = 1
            _positions = (bp2, bp1, mode2, mode1)
            strands = (strand2, strand1)
            genes = list(reversed(genes))
            return (
                sv_type,
                annot,
                canonical,
                _positions,
                read2_info,
                read1_info,
                strands,
                genes,
            )
        else:
            return event

    @staticmethod
    def order_events_by_trancription_direction(event_list):
        """
        construct breakpoints order following transcription direction for multiple-hop events or one-hop events
                bp1                bp2   bp3               bp4
        ---------|------    -------|----|------    --------|---------
              Node1                 Node2                Node3
        ..note ::
               requirements
               * bp2 and bp3 at the same chromosome
               * if strand(+): bp3 > bp2
                 if strand(-): bp3 < bp2
        """
        kept_right_pos = None
        kept_right_chrm = None
        kept_right_strand = None

        keep_event_list_order = True
        output_event_list = []
        for evt in event_list:
            ordered_evt = Series.reorder_event(evt)
            output_event_list.append(ordered_evt)
            (
                sv_type,
                annot,
                canonical,
                _positions,
                read1_info,
                read2_info,
                strands,
                genes,
            ) = ordered_evt
            chrm1, _pos1 = _positions[0].split(":")
            chrm2, _pos2 = _positions[1].split(":")
            pos1 = int(_pos1)
            pos2 = int(_pos2)
            strand1 = strands[0]
            strand2 = strands[1]
            if not kept_right_pos:
                kept_right_pos = pos2
                kept_right_chrm = chrm2
                kept_right_strand = strand2
            else:
                if chrm1 == kept_right_chrm and strand1 == kept_right_strand:
                    if strand1 == "+" and pos1 > kept_right_pos:
                        kept_right_pos = pos2
                        kept_right_chrm = chrm2
                        kept_right_strand = strand2
                    elif strand1 == "-" and pos1 < kept_right_pos:
                        kept_right_pos = pos2
                        kept_right_chrm = chrm2
                        kept_right_strand = strand2
                    else:
                        keep_event_list_order = False
                else:
                    keep_event_list_order = False

        if not keep_event_list_order:
            output_event_list = list(reversed(output_event_list))

        return output_event_list

    def __getitem__(self, index):
        return self.nodes[index]

    def __hash__(self) -> int:
        return hash(";".join(map(str, self.nodes)))

    def __eq__(self, other) -> bool:
        return ";".join(map(str, self.nodes)) == ";".join(map(str, other.nodes))

    def __len__(self) -> int:
        return len(self.nodes)

    def __lt__(self, other) -> bool:
        return len(self.nodes) < len(other.nodes)

    def __repr__(self) -> str:
        return ";".join(map(str, self.nodes))

    def decompose(self) -> list:
        """Decompose the sequence of Nodes into Nodes pair"""
        paired_breakpoints = []
        for i, j in zip(self.nodes[::1], self.nodes[1::1]):
            paired_breakpoints.append(
                f"{i.sv_type}-{i.next_breakpoint}-{j.prev_breakpoint}-{i.strand}-{j.strand}"
            )
        return paired_breakpoints


class Blat(object):
    """
    the Blat class is used to integrate the blat service (gfServer and
    gfClient) so that we can query certain sequences from the genome shamelessly
    """

    def __init__(
        self, ref_2bit: str, logger: logger, port: int, output_dir: str
    ) -> None:
        """
        :param ref_2bit: the path of reference for blat alignment
        :param logger: the logger for logging
        :param port: the port of server service for blat alignment
        :param output_dir: the path for storing alignment result
        """
        self.port, self.ref_2bit = port, ref_2bit
        self.output_dir = output_dir
        self.ran_id = random.getrandbits(30)
        self.is_start_server = True
        self.logger = logger

    @property
    def ref_dir(self) -> str:
        """
        the property for ref_dir, which is the path of reference for blat

        :return: the absolute path of reference dir
        """
        if self.ref_2bit.startswith("~"):
            abs_2bit = os.path.join(
                os.path.expanduser("~"), self.ref_2bit.replace("~/", "")
            )
            ref_dir = os.path.dirname(abs_2bit)
        else:
            abs_2bit = os.path.abspath(self.ref_2bit)
            ref_dir = os.path.dirname(abs_2bit)
        return ref_dir

    @property
    def log_file(self) -> str:
        """
        the property for log_file, which is the path of log file for blat
        """
        return f"{self.ref_dir}/gfserver.temp.{self.ran_id}.log"

    def is_ready(self) -> bool:
        """
        the function for checking whether the blat server is ready or not
        after starting the server service
        :return: the boolean value of whether the server is ready or not
        """
        flag = False
        self.logger.debug("check if the server starts")
        if os.path.exists(self.log_file):
            with open(self.log_file) as f:
                for line in f:
                    if "Server ready" in line:
                        flag = True
        return flag

    def is_running(self) -> bool:
        """
        the function for checking whether the blat server is running or not

        :return: the boolean value of whether the server is running or not
        """
        return True if self._search_processing() else False

    def _search_processing(self) -> List:
        """
        the function for searching the process of blat server
        in current system
        :return: the list of process of blat server
        """
        result = []
        self.logger.debug("searching server service")
        for proc in psutil.process_iter(["pid", "name"]):
            if "gfServer".lower() == proc.name().lower():
                if proc.cmdline():
                    result.append(proc)
        return result

    def _run_cmd(self, cmd: str) -> None:
        """
        the function is used to run the command in the system
        :param cmd: the command to be run
        """
        subprocess.run(cmd.split(), check=True)

    def _start_server(self) -> Process:
        """gfServer should run at the directory where gfServer,
        gfClient and hg38.2bit located"""

        cwd = os.path.abspath(os.getcwd())

        # change to use_blat directory
        os.chdir(self.ref_dir)

        if os.path.exists(self.log_file):
            os.remove(self.log_file)

        cmd = f"gfServer -canStop -log={self.log_file} -stepSize=5 start localhost {self.port} {self.ref_2bit}"
        process = Process(target=self._run_cmd, args=[cmd])
        process.start()
        self.logger.debug("starting server service")
        os.chdir(cwd)
        return process

    def start_server(self) -> None:
        """
        the function for starting the server service, if the server is not running,
        we will start the server service
        """
        running_flag = self.is_running()
        if not running_flag:
            self._start_server()
        else:
            self.is_start_server = False

    def stop_server(self) -> None:
        """
        the function for stopping the server service, if the server is running,
        """
        procs = self._search_processing()
        self.logger.debug("stopping server service")
        for proc in procs:
            proc.kill()

    def _query(self, in_seq: str, miniIdentity: int = 90) -> str:
        """
        the function is help function in order to using gfClient
        to query 'in_seq' to generate alignment file (in PSL format).

        :param miniIdentity: the threshold of the identity for aligning
        :param in_seq: sequence of softclipped segment
        :return: the path for PSL file
        """
        self.logger.debug("querying the sequence")
        ran_id = random.getrandbits(30)
        in_fasta = os.path.join(self.output_dir, "{}.fasta".format(ran_id))
        with open(in_fasta, "w", buffering=1) as fasta_file:
            fasta_file.write(">{}\n".format(ran_id))
            fasta_file.write("{}\n".format(in_seq))

        out_psl = os.path.join(self.output_dir, "{}.psl".format(ran_id))

        cwd = os.path.abspath(os.getcwd())

        os.chdir(self.ref_dir)
        cmd = "gfClient -minScore=20 -minIdentity={} localhost {} {} {} {} > /dev/null".format(
            miniIdentity, self.port, self.ref_dir, in_fasta, out_psl
        )
        try:
            ret = subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as err:
            raise SystemExit(f"{err} {err.output}")

        os.chdir(cwd)
        self._remove(in_fasta)

        return out_psl

    def _wait_ready(self, interval: int = 30) -> None:
        """
        the function for waiting the server service to be ready,

        :param interval: the interval time for checking the server service
        """
        while not self.is_ready():
            time.sleep(interval)

    def query(self, in_seq: str, miniIdentity: int = 90) -> str:
        """
        the function for querying the sequence to the server service

        :param in_seq: the sequence of input sequence
        :param miniIdentity: the threshold of the identity for aligning
        :return: the path for PSL file
        """

        if self.is_start_server:
            if self.is_ready():
                out_psl = self._query(in_seq, miniIdentity)
            else:
                self._wait_ready()
                out_psl = self._query(in_seq, miniIdentity)
        else:
            out_psl = self._query(in_seq, miniIdentity)

        return out_psl

    def query_insertion(
        self,
        insert_seq: str,
        threshold_identity: float = 0.99,
        top: int = 3,
        align_len_threshold: int = 20,
    ) -> Any:
        """
        the function for querying the insertion sequence to the server service, and
        the function is a specific version of the function 'query'.

        :param insert_seq: insertion sequence
        :param threshold_identity: the threshold of the identity for aligning
        :param top: the top number of the alignments
        :param align_len_threshold: the threshold of the insertion sequence length
        :return: insertion sequence alignment in NamedTuple format
        """

        flag = False  # flag for checking the insertion  if its hit is only one

        if len(insert_seq) < align_len_threshold:
            return flag, NoneInsertion(hit_num=0, query_sequence=insert_seq)

        out_blat = self.query(in_seq=insert_seq)
        try:
            blat = SearchIO.read(out_blat, "blat-psl")
        except ValueError:
            return flag, NoneInsertion(hit_num=0, query_sequence=insert_seq)

        hsps = blat.hsps
        hsps.sort(key=lambda x: x.score, reverse=True)
        hsps = hsps[:top]
        keep_hsp = []
        for hsp in hsps:
            if sum(hsp.hit_span_all) / len(insert_seq) > threshold_identity:
                keep_hsp.append(hsp)
        hit = len(keep_hsp)

        if hit == 1:
            top_hsp = keep_hsp[0]
            flag = True
            # start_end = top_hsp.hit_range_all
            # strand = "+" if top_hsp.hit_strand_all[0] == 1 else "-"
            # chrom = top_hsp.hit_id

            ref_chrom, position, strand, cigar, num_of_mismatch = self.psl2sam(
                top_hsp, in_seq_len=len(insert_seq)
            )
            return flag, Insertion(
                hit_num=1,
                chrom=ref_chrom,
                ref_start=position,
                strand=strand,
                cigarstring=cigar,
                mapq=60,
                nm=num_of_mismatch,
                query_sequence=insert_seq,
            )
        else:
            return flag, NoneInsertion(hit_num=hit, query_sequence=insert_seq)

    @staticmethod
    def _remove(file):
        """
        the function for removing the file

        :param file: the path of the file
        """
        if os.path.exists(file):
            os.remove(file)

    @staticmethod
    def _calculate_mapq(hsps: Any, in_seq_len: int, threshold_identity: float) -> int:
        """
        the function is used to calculate map quality of the insertion.

        :param hsps: the list of hsp after aligning the insertion sequence
        :param in_seq_len: the length of the input sequence
        :param threshold_identity: the threshold of the identity for aligning
        :return: the map quality of the insertion
        """
        num_of_locations = 0

        for hsp in hsps:
            if (
                hsp.ident_pct / 100 >= threshold_identity
                and hsp.query_span / in_seq_len >= threshold_identity
            ):
                num_of_locations += 1
        if num_of_locations == 1:
            mapq = 60
        elif num_of_locations == 2:
            mapq = 3
        elif num_of_locations == 3:
            mapq = 2
        elif 4 <= num_of_locations <= 9:
            mapq = 1
        else:
            mapq = 0
        return mapq

    def fetch_mapq(self, in_seq: str, threshold_identity: float) -> Any:
        """
        the function is used to fetch the map quality of the insertion.

        :param in_seq: the input sequence
        :param threshold_identity: the threshold of the identity for aligning
        :return: the top hit of the insertion sequence, and the map quality of the insertion
        """
        psl_file = self.query(in_seq=in_seq)

        try:
            blat = SearchIO.read(psl_file, "blat-psl")
        except ValueError:
            self.logger.error(f"No Blat hit found {in_seq}")
            raise SystemExit
        else:
            hsps = blat.hsps
            hsps.sort(key=lambda k: k.score, reverse=True)
            top_hsp = hsps[0]
            Blat._remove(psl_file)
            mapq = Blat._calculate_mapq(hsps, len(in_seq), threshold_identity)
        return top_hsp, mapq

    def psl2sam(self, hsp: Any, in_seq_len: int) -> Tuple[str, int, str, str, int]:
        """
        Convert the top HSP in PSL file to SAM fields chrom, reference_start,
        strand, cigarstring, num_of_mismatch. The function try to implement
        the psl2sam.pl script and return the cigar and mapping position
        estimated from psl file

        :param hsp: the selected HSP form BLAT
        :param in_seq_len: the length of the input sequence
        :return: chrom, reference_start, strand, cigarstring, num_of_mismatch
        """

        cigar = ""
        query_start = hsp.query_start
        query_end = hsp.query_end

        _strand = hsp.query_strand_all[0]  # may need replace by query_strand
        ref_start, ref_end = hsp.hit_range
        ref_chrom = hsp.hit_id
        num_of_mismatch = hsp.mismatch_num

        soft_len = 0
        if _strand == -1:
            query_start = in_seq_len - hsp.query_end
            query_end = in_seq_len - hsp.query_start
        if query_start:
            # 5'-end clipping
            soft_len = query_start
            cigar += str(query_start) + "S"
        x = hsp.query_span_all
        if _strand == -1:
            y = [
                in_seq_len - item[1] for item in hsp.query_range_all
            ]  # may need replace by query_start_all when the bug is fixed in Biopython
        else:
            y = [
                item[0] for item in hsp.query_range_all
            ]  # may need replace by query_start_all when the bug is fixed in Biopython
        z = hsp.hit_start_all
        y0, z0 = y[0], z[0]
        for i in range(1, len(hsp)):
            ly = y[i] - y[i - 1] - x[i - 1]
            lz = z[i] - z[i - 1] - x[i - 1]
            if ly < lz:
                # del: the reference gap is longer
                cigar += str(y[i] - y0) + "M"
                if lz - ly >= 10:
                    cigar += str(lz - ly) + "N"
                else:
                    cigar += str(lz - ly) + "D"
                y0, z0 = y[i], z[i]
            elif lz < ly:
                # ins: the query gap is longer
                cigar += str(z[i] - z0) + "M"
                cigar += str(ly - lz) + "I"
                y0, z0 = y[i], z[i]

        cigar += str(query_end - y0) + "M"
        # print(cigar)
        # return cigar, soft_len
        if in_seq_len != query_end:
            # 3'-end clipping
            end3 = in_seq_len - query_end
            if end3 > soft_len:
                soft_len = end3
            cigar += str(end3) + "S"
        # return cigar, soft_len
        strand = "+" if _strand == 1 else "-"

        return ref_chrom, ref_start + 1, strand, cigar, num_of_mismatch


class ReadsConnecter(object):
    """
    the ReadsConnecter class is used to connect the reads and identify the mode of the reads
    """

    def __init__(
        self,
        aln_list: List[Read],
        blat: Blat,
        logger: logger,
        soft_len_cutoff: int = 30,
    ) -> None:
        self.reads_chain, self.candidate_nodes = [], []
        self.read_pair_mode_dict, self.insertion_dict = {}, {}
        self.aln_list = aln_list
        self.soft_len_cutoff = soft_len_cutoff
        self.logger = logger
        self.blat = blat

    @staticmethod
    def init_mode_judge(sms: Tuple[int, int, int]) -> int:
        _lt, _, _rt = sms
        # SM
        if _lt > _rt:
            return 2
        # MS
        else:
            return 1

    @staticmethod
    def conduct_glocal_alignment_forMS(
        query_seq: str,
        target_seq: str,
        same_strand: bool,
        is_align: bool,
        s_position: str,
        threshold: float = 0.7,
    ) -> Tuple[bool, Union[None, str]]:
        """query_seq: M  target_seq: S"""
        # do not conduct alignment

        insert_seq = None  # None means M is not consist with S
        match_flag = False
        if not is_align:
            #
            # insert_len = len(target_seq) - len(query_seq)
            #
            # insert_seq = (
            #     target_seq[-insert_len:]
            #     if s_position == "left"
            #     else target_seq[:insert_len]
            # )
            #
            # local_alignment_result = aligner(insert_seq, query_seq, method="local")[0]
            #
            # if (
            #         local_alignment_result.start2 == 0
            #         or local_alignment_result.end2 == len(query_seq)
            # ):
            #
            #     return True, None
            #
            # else:
            #     return True, insert_seq
            return True, insert_seq

        if not same_strand:
            target_seq = str(Seq(target_seq).reverse_complement())

        alignment_result = aligner(query_seq, target_seq, method="semi-global")[0]
        _query_seq = alignment_result.seq1.decode("utf-8")
        _target_seq = alignment_result.seq2.decode("utf-8")

        _query_seq_len, _target_seq_len = len(_query_seq), len(_target_seq)

        query_identity = 1 - (
            len(query_seq)
            - _query_seq_len
            + alignment_result.n_gaps1
            + alignment_result.n_gaps2
            + alignment_result.n_mismatches
        ) / len(query_seq)

        if query_identity > threshold:
            match_flag = True
            #
            # if len(query_seq) >= len(target_seq):
            #     return match_flag, insert_seq
            #
            # if s_position == "left":
            #     insert_len = len(_query_seq) - len(_query_seq.rstrip("-"))
            # else:
            #     insert_len = alignment_result.start2
            #
            # if insert_len > 0:
            #     local_len = int(0.25 * len(_query_seq.rstrip("-"))) + insert_len
            #     local_query_seq = _query_seq[-insert_len - local_len : -insert_len]
            #     local_target_seq = _target_seq[-local_len:]
            #     local_alignment_result = aligner(
            #         local_query_seq, local_target_seq, method="local"
            #     )
            #     if local_alignment_result[0].end2 < local_len:
            #         insert_seq = (
            #             target_seq[:insert_len]
            #             if s_position == "left"
            #             else target_seq[-insert_len:]
            #         )

        return match_flag, insert_seq

    def test_4case(
        self, start_read: Read, read: Read, is_align_for_ms: bool
    ) -> Tuple[bool, Read]:

        _lt_len_r1, _read_match_r1, _rt_len_r1 = start_read.adhocsms
        _lt_len_r2, _read_match_r2, _rt_len_r2 = read.sms

        self.logger.debug(f"{start_read.mode}, {read.mode}")

        self.logger.debug(f"{start_read.adhocsms}, {read.sms}")

        same_strand = True if start_read.adhocseq == read.query_sequence else False

        # first case
        self.logger.debug("testing first case M vs LS")
        match_flag, insertion_1_seq = ReadsConnecter.conduct_glocal_alignment_forMS(
            start_read.adhocseq[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
            read.query_sequence[:_lt_len_r2],
            same_strand,
            is_align_for_ms,
            "left",
        )

        if match_flag:  # may same
            # if insertion_1_seq is not None:  # insertion exist
            #     insertion = self.map_with_blat_for_genome(insertion_1_seq)
            #     self.insertion_dict[(start_read, read)] = insertion

            read.mode = 2
            self.logger.debug(f"{start_read.mode}, {read.mode}")
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)

            if read in self.candidate_nodes:
                self.candidate_nodes.remove(read)

            start_read = read
            start_read.adhocsms = 0, _lt_len_r2 + _read_match_r2, _rt_len_r2
            start_read.adhocseq = read.query_sequence

            return True, start_read

        self.logger.debug("testing second case M vs RS")
        # second case
        match_flag, insertion_2_seq = ReadsConnecter.conduct_glocal_alignment_forMS(
            start_read.adhocseq[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
            read.query_sequence[-_rt_len_r2:],
            same_strand,
            is_align_for_ms,
            "right",
        )

        if match_flag:

            read.mode = 1

            self.logger.debug(f"{start_read.mode}, {read.mode}")

            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)

            self.reads_chain.append(read)

            if read in self.candidate_nodes:
                self.candidate_nodes.remove(read)

            start_read = read
            start_read.adhocsms = (
                _lt_len_r2,
                _read_match_r2 + _rt_len_r2,
                0,
            )
            start_read.adhocseq = read.query_sequence

            return True, start_read

        self.logger.debug("testing third case LS vs M")
        # third case
        match_flag, insertion_3_seq = ReadsConnecter.conduct_glocal_alignment_forMS(
            read.query_sequence[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
            start_read.adhocseq[:_lt_len_r1],
            same_strand,
            is_align_for_ms,
            "left",
        )

        if match_flag:
            read, start_read = start_read, read

            read.mode = 1
            self.logger.debug(f"{read.mode}, {start_read.mode}")
            self.read_pair_mode_dict[(read, start_read)] = (read.mode, start_read.mode)

            self.reads_chain.append(start_read)

            if start_read in self.candidate_nodes:
                self.candidate_nodes.remove(start_read)

            read.adhocsms = 0, _lt_len_r1 + _read_match_r1, _rt_len_r1

            return True, read

        self.logger.debug("testing fourth case RS vs M")
        # fourth case
        match_flag, insertion_4_seq = ReadsConnecter.conduct_glocal_alignment_forMS(
            read.query_sequence[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
            start_read.adhocseq[-_rt_len_r1:],
            same_strand,
            is_align_for_ms,
            "right",
        )

        if match_flag:
            read, start_read = start_read, read

            start_read.mode = 2

            self.logger.debug(f"{read.mode}, {start_read.mode}")
            self.read_pair_mode_dict[(read, start_read)] = (read.mode, start_read.mode)

            self.reads_chain.append(start_read)

            if start_read in self.candidate_nodes:
                self.candidate_nodes.remove(start_read)

            read.adhocsms = (
                _lt_len_r1,
                _rt_len_r1 + _read_match_r1,
                0,
            )

            return True, read

    def run(self) -> None:
        """Find the best connected paths for a list of chimeric alignments
        .. note::
            Read-to-Read chain scenarios
            * [[Read1, Read2, Read3]]
            * [[Read1, Read2, Read3],[Read4,Read5]]

            Dictionary of Read-pair scenarios
            * (Read1, Read2) => mode-of-Read1, mode-of-Read2
            * (Read2, Read1) => mode-of-Read2, mode-of-Read1
        """
        start_nodes = []

        # find start node and end node
        for read in self.aln_list:
            if (
                read.lt_soft_len < self.soft_len_cutoff
                or read.rt_soft_len < self.soft_len_cutoff
            ):
                start_nodes.append(read)
            else:
                self.candidate_nodes.append(read)

        start_read = start_nodes[0]
        end_read = start_nodes[1]

        start_read.adhocsms = start_read.sms
        start_read.adhocseq = start_read.query_sequence

        self.reads_chain.append(start_read)

        if not self.candidate_nodes:  # []

            self.logger.debug("ReadConnecter: candidate_nodes is []")
            start_read.mode, end_read.mode = (
                ReadsConnecter.init_mode_judge(start_read.adhocsms),
                ReadsConnecter.init_mode_judge(end_read.sms),
            )

            _, start_read = self.test_4case(start_read, end_read, is_align_for_ms=False)

        else:
            for read in self.candidate_nodes:
                start_read.mode, read.mode = (
                    ReadsConnecter.init_mode_judge(start_read.adhocsms),
                    ReadsConnecter.init_mode_judge(read.sms),
                )
                _, start_read = self.test_4case(start_read, read, is_align_for_ms=True)

            start_read.mode, end_read.mode = (
                ReadsConnecter.init_mode_judge(start_read.adhocsms),
                ReadsConnecter.init_mode_judge(end_read.sms),
            )
            _, start_read = self.test_4case(start_read, end_read, is_align_for_ms=True)
