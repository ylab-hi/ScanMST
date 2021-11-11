#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
import argparse
import re


class Path(object):
    """store chimeric reads as nodes in a path (directed acyclic graph)
    :param nodes: a list of Read as nodes
    :type nodes: Read
    :param sms: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length)
    :type sms: tuple
    :param sequence: reads sequence of last added read
    :type sequence: str
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
        "query_length",
    )

    def __init__(
        self,
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
    ) -> None:
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

    @classmethod
    def init(cls, chrom, position, strand, cigar_str, mapq, nm, query_seq):
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
        # return f'{self.lt_soft_len}\t{self.read_match_size}\t{self.rt_soft_len}'
        return (self.lt_soft_len, self.read_match_size, self.rt_soft_len)

    def add_path(self, path) -> None:
        """path is an instance of Path class"""
        self.linked_paths.append(path)

    def get_exons_and_introns(self) -> tuple:
        """get the coordiantes for reads matched part (without softclipping)
        :return: exons coordiantes and introns coordiantes
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
    ) -> None:
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
        return [cls() for i in range(number)]


class Series(object):
    """construct a sequence of Nodes for storing information of connected breakpoints
    :param nodes: sequence of Nodes
    :type nodes: list
    :param assemblied: The series is from assembly of reads (True) or a single read (False)
    :type assemblied: bool

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

    __slots__ = ("nodes", "assemblied")

    def __init__(self) -> None:
        self.nodes = []
        self.assemblied = None

    def init(self, event_list) -> None:
        """add event list as Node to self.nodes"""
        event_list = self.order_events_by_trancription_direction(event_list)
        hop_number = len(event_list)
        # print('hop_number:', hop_number)
        self.nodes = Node.create_nodes(hop_number + 1)

        for i in range(hop_number):
            (
                sv_type,
                annot,
                canonical,
                _positions,
                read1_info,
                read2_info,
                strands,
                genes,
            ) = event_list[i]
            _bp1 = _positions[0]
            _bp2 = _positions[1]
            _chrm1 = _bp1.split(":")[0]
            _chrm2 = _bp2.split(":")[0]
            _mode1 = _positions[2]
            _mode2 = _positions[3]
            _strand1 = strands[0]
            _strand2 = strands[1]
            read1_ref_start, read1_ref_end, read1_exons = read1_info
            read2_ref_start, read2_ref_end, read2_exons = read2_info

            self.nodes[i].next_breakpoint = _bp1
            self.nodes[i + 1].prev_breakpoint = _bp2

            self.nodes[i].strand = _strand1
            self.nodes[i + 1].strand = _strand2

            self.nodes[i].chrom = _chrm1
            self.nodes[i + 1].chrom = _chrm2

            self.nodes[i].ref_start = read1_ref_start
            self.nodes[i + 1].ref_start = read2_ref_start

            self.nodes[i].ref_end = read1_ref_end
            self.nodes[i + 1].ref_end = read2_ref_end

            self.nodes[i].exons = read1_exons
            self.nodes[i + 1].exons = read2_exons

            self.nodes[i].sv_type = sv_type
            self.nodes[i].annotation_code = annot
            self.nodes[i].splicing_code = canonical
            self.nodes[i].modes = (_mode1, _mode2)
            self.nodes[i].genes = genes

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

    def order_events_by_trancription_direction(self, event_list):
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

    # def __hash__(self) -> int:
    #    hash_val = 0
    #    for i in self.nodes:
    #        hash_val ^= hash(i)
    #    return hash_val

    def __len__(self) -> int:
        return len(self.nodes)

    def __lt__(self, other) -> bool:
        return len(self.nodes) < len(other.nodes)

    def __repr__(self) -> str:
        return ";".join(map(str, self.nodes))

    # def reversed(self) -> None:
    #    """reverse the sequence of Nodes"""
    #    self.nodes = self.nodes[::-1]

    def decompose(self) -> list:
        """Decompose the sequence of Nodes into Nodes pair"""
        paired_breakpoints = []
        for i, j in zip(self.nodes[::1], self.nodes[1::1]):
            paired_breakpoints.append(
                f"{i.sv_type}-{i.next_breakpoint}-{j.prev_breakpoint}-{i.strand}-{j.strand}"
            )
        return paired_breakpoints

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
