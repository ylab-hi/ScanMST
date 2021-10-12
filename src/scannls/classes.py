#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
from typing import Iterable
from pyfaidx import Fasta
import argparse
import re

class Path(object):
    '''store chimeirc reads as nodes in a path (directed acyclic graph)
    :param nodes: a list of Read as nodes
    :type nodes: Read
    :param sms: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length)
    :type sms: tuple
    :param nm: summation of number-of-mismatches of Reads in the path
    :type nm: int
    :param mode: a dictionary of Reads-pair to mode in the path
    :type mode: dict
    .. note::
        We have to pay attention on 'sms':
            * every path only have keep one 'sms' value (per path instead of per node).
            * once new node added to the path, update the value of 'sms' using the summed 'sms' value from the function 'test_is_connected'
    '''
    __slots__ = (
        "nodes",
        "sms",
        "nm",
        "mode"
    )

    def __init__(self) -> None:
        self.nodes = []
        self.sms = None
        self.nm = 0
        self.mode = {}

    def add(self, read) -> None:
        ''' add one chimeric read to the path
        '''
        self.nodes.append(read)
        self.nm = self.nm + read.nm

    def add_mode(self, read_pair_dict) -> None:
        ''' add read-pair=> mode to the path
        '''
        self.mode.update(read_pair_dict)

    def __len__(self) -> int:
        return len(self.nodes)

    def __lt__(self, other) -> bool:
        return len(self.nodes) < len(other.nodes)

    def __repr__(self) -> str:
        return ";".join( map(str, self.nodes) )

    def __hash__(self) -> int:
        return hash( ";".join( map(str, self.nodes) ) )

    def __eq__(self, other) -> bool:
        return ";".join( map(str, self.nodes) ) == ";".join( map(str, other.nodes) )


class Read(object):
    """build a read class for storing information of every junction read
    :param chrom: chromosome of genome
    :type chrom: str
    :param position: start position of chimeric read
    :type start: int
    :param strand: direction of chimeric read (-|+)
    :type strand: str
    :param cigar_str: cigar string of chimeric read (-|+)
    :type cigar_str: str
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
        "query_length"
    )

    def __init__(self, chrom, position, strand, cigar_str, mapq, nm, query_seq, lt_soft_len, rt_soft_len, read_match_size, reference_match_size, indel_size, cigar_without_soft, query_length, cigartuples) -> None:
        self.chrom = chrom
        self.ref_start = position
        self.strand = strand
        self.cigarstring = cigar_str
        self.mapq = mapq
        self.nm = nm
        self.query_sequence = query_seq
        self.linked_paths = []
        self.lt_soft_len = lt_soft_len
        self.rt_soft_len = rt_soft_len
        self.read_match_size = read_match_size
        self.reference_match_size = reference_match_size
        self.indel_size = indel_size
        self.cigartuples_without_soft = cigar_without_soft
        self.query_length = query_length
        self.cigartuples = cigartuples

    def __eq__(self, other) -> bool:
        if isinstance(other, Read):
            if self.chrom == other.chrom and self.ref_start == other.ref_start \
                    and self.ref_end == other.ref_end and self.strand == other.strand \
                    and self.mapq == other.mapq and self.nm == other.nm:
                return True
        return False

    def __hash__(self) -> int:
        return hash(self.chrom) ^ hash(self.ref_start) ^ hash(self.ref_end) ^ hash(self.strand) ^ hash(self.mapq) ^ hash(self.nm)

    def __lt__(self, other) -> bool:
        if self.chrom == other.chorm:
            return self.ref_start < other.ref_start
        else:
            chrm_dict = {'chrM':0, 'MT':0, 'chrX':23, 'chrY':24, 'X':23, 'Y':24}
            for i in range(1,23):
                chrm_dict.update({f'chr{i}':i})
                chrm_dict.update({f'{i}':i})
            return chrm_dict[self.chrom] < chrm_dict[other.chrom]

    # for debug purpose
    def __repr__(self) -> str:
        return (
            fr"Read({self.chrom}, {self.ref_start}, {self.ref_end}, {self.strand}, {self.mapq}, {self.nm})"
        )

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
                    f'{self.lt_soft_len}:{self.read_match_size}:{self.rt_soft_len}'
                ],
            ),
        )


    @classmethod
    def init(cls, chrom, position, strand, cigar_str, mapq, nm, query_seq):
        cigar_char_dict = {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
        # 'length', 'operation char'
        len_type_tuple = re.findall(r'(\d+)(\w)', cigar_str)
        # (operation code, length)
        cigartuples = [(cigar_char_dict[j], int(i)) for i,j in len_type_tuple]

        query_length = 0
        indel_size = 0
        reference_match_size = 0
        read_match_size = 0
        cigar_without_soft = []
        for op_code, _len_ in cigartuples:
            if op_code == 0:# M
                reference_match_size += _len_
                read_match_size += _len_
                query_length += _len_
                cigar_without_soft.append([0, _len_])
            elif op_code == 1:# I
                indel_size += -_len_
                read_match_size += _len_
                query_length += _len_
                cigar_without_soft.append([1, _len_])
            elif op_code == 2:# D
                indel_size += _len_
                reference_match_size += _len_
                cigar_without_soft.append([2, _len_])
            elif op_code == 3:# N
                indel_size += _len_
                reference_match_size += _len_
                cigar_without_soft.append([3, _len_])
            elif op_code == 4:# S
                query_length += _len_

        lt_soft_len = 0
        rt_soft_len = 0
        lt_op, lt_len = cigartuples[0]
        rt_op, rt_len = cigartuples[-1]
        if lt_op == 4:
            lt_soft_len = lt_len
        if rt_op == 4:
            rt_soft_len = rt_len
        return cls(chrom, position, strand, cigar_str, mapq, nm, query_seq, lt_soft_len, rt_soft_len, read_match_size, reference_match_size, indel_size, cigar_without_soft, query_length, cigartuples)


    @property
    def ref_end(self) -> int:
        return self.ref_start + self.reference_match_size

    @property
    def reference_span(self) -> int:
        ''' M+N+D
        '''
        return self.reference_match_size

    @property
    def sms(self) -> tuple:
        #return f'{self.lt_soft_len}\t{self.read_match_size}\t{self.rt_soft_len}'
        return (self.lt_soft_len, self.read_match_size, self.rt_soft_len)

    def add_path(self, path) -> None:
        '''path is an instance of Path class
        '''
        self.linked_paths.append(path)

    def splice_site_checker(self, genome_fasta, fraction_cutoff=0.6) -> bool:
        ''' check whether the fraction of canonical splice site usage in read reference matched part is bigger than 'fraction_cutoff' or not
        :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
        :param fraction_cutoff: fraction of canonical splice sites used in the putative introns inferred from the CIGAR
        :type genome_fasta: pyfaidx.Fasta
        :type fraction_cutoff: float
        :return: using canonical splice sites OR not
        :rtype: bool
        '''
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start
        for op_code, _len_ in self.cigartuples_without_soft:
            if op_code in {0, 2}: # M, D
                current_pos = current_pos + _len_
            elif op_code == 3:# N
                exons.append((start_pos,  current_pos))
                current_pos = current_pos + _len_
                start_pos = current_pos
        exons.append((start_pos,  current_pos))

        # No 'N' in the cigar
        if len(exons) == 1:
            return True
        elif len(exons) > 1:
            _positions = []
            for i,j in exons:
                _positions.extend([i, j])
            _positions.sort()
            _positions.pop(0)
            _positions.pop(-1)
            intron_positions = zip(_positions[::2], _positions[1::2])

            intron_count = 0
            can_count = 0
            can_sites = {'GT-AG','GC-AG','AT-AC'}
            for start,end in intron_positions:
                if end - start >= 10:
                    intron_count += 1
                    if self.strand == '-':
                        left_site = genome_fasta[self.chrom][end-2:end].reverse.complement.seq
                        right_site = genome_fasta[self.chrom][start:start+2].reverse.complement.seq
                    else:
                        left_site = genome_fasta[self.chrom][start:start+2].seq
                        right_site = genome_fasta[self.chrom][end-2:end].seq
                    if f'{left_site}-{right_site}' in can_sites:
                        can_count += 1
            if can_count/intron_count >= fraction_cutoff:
                return True
            else:
                return False

class LengthAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values <= 0:
            parser.error("Minimum length for {0} is 1".format(option_string))
        setattr(namespace, self.dest, values)

