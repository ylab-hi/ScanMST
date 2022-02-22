"""Cpp extension for BAM file parser"""
from __future__ import annotations
import typing


class Aligner:
    def SetReferenceSequence(self, arg0: str, arg1: int) -> int: ...

    @typing.overload
    def __init__(self) -> None: ...

    @typing.overload
    def __init__(self, arg0: int, arg1: int, arg2: int, arg3: int) -> None: ...

    pass


class Alignment:
    def Clear(self) -> None: ...

    def __init__(self) -> None: ...

    @property
    def best_score(self) -> int:
        """
        :type: int
        """

    @best_score.setter
    def best_score(self, arg0: int) -> None:
        pass

    @property
    def best_score2(self) -> int:
        """
        :type: int
        """

    @best_score2.setter
    def best_score2(self, arg0: int) -> None:
        pass

    @property
    def cigar(self) -> typing.List[int]:
        """
        :type: typing.List[int]
        """

    @cigar.setter
    def cigar(self, arg0: typing.List[int]) -> None:
        pass

    @property
    def cigar_string(self) -> str:
        """
        :type: str
        """

    @cigar_string.setter
    def cigar_string(self, arg0: str) -> None:
        pass

    @property
    def mismatches(self) -> int:
        """
        :type: int
        """

    @mismatches.setter
    def mismatches(self, arg0: int) -> None:
        pass

    @property
    def query_begin(self) -> int:
        """
        :type: int
        """

    @query_begin.setter
    def query_begin(self, arg0: int) -> None:
        pass

    @property
    def query_end(self) -> int:
        """
        :type: int
        """

    @query_end.setter
    def query_end(self, arg0: int) -> None:
        pass

    @property
    def ref_end_next_best(self) -> int:
        """
        :type: int
        """

    @ref_end_next_best.setter
    def ref_end_next_best(self, arg0: int) -> None:
        pass

    @property
    def reference_begin(self) -> int:
        """
        :type: int
        """

    @reference_begin.setter
    def reference_begin(self, arg0: int) -> None:
        pass

    @property
    def reference_end(self) -> int:
        """
        :type: int
        """

    @reference_end.setter
    def reference_end(self, arg0: int) -> None:
        pass

    pass


class Filter:
    @typing.overload
    def __init__(self) -> None: ...

    @typing.overload
    def __init__(self, arg0: bool, arg1: bool, arg2: int, arg3: int) -> None: ...

    @property
    def distance_filter(self) -> int:
        """
        :type: int
        """

    @distance_filter.setter
    def distance_filter(self, arg0: int) -> None:
        pass

    @property
    def report_begin_position(self) -> bool:
        """
        :type: bool
        """

    @report_begin_position.setter
    def report_begin_position(self, arg0: bool) -> None:
        pass

    @property
    def report_cigar(self) -> bool:
        """
        :type: bool
        """

    @report_cigar.setter
    def report_cigar(self, arg0: bool) -> None:
        pass

    @property
    def score_filter(self) -> int:
        """
        :type: int
        """

    @score_filter.setter
    def score_filter(self, arg0: int) -> None:
        pass

    pass


class Rescuer:
    @typing.overload
    def __init__(self) -> None: ...

    @typing.overload
    def __init__(self, arg0: str, arg1: int, arg2: int, arg3: int, arg4: float) -> None: ...

    def __repr__(self) -> str: ...

    @staticmethod
    @typing.overload
    def calculate_sr(chrom, start, end, mode, current_names, names_in_graph) -> int: ...


class parseCigarResult:
    def __repr__(self) -> str: ...

    @property
    def cigartuples(self) -> typing.List[int]:
        """
        :type: typing.List[int]
        """

    @property
    def cigartuples_without_soft(self) -> typing.List[int]:
        """
        :type: typing.List[int]
        """

    @property
    def indel_len(self) -> int:
        """
        :type: int
        """

    @property
    def lt_soft_len(self) -> int:
        """
        :type: int
        """

    @property
    def query_len(self) -> int:
        """
        :type: int
        """

    @property
    def read_match(self) -> int:
        """
        :type: int
        """

    @property
    def ref_match(self) -> int:
        """
        :type: int
        """

    @property
    def rt_soft_len(self) -> int:
        """
        :type: int
        """

    pass


def parseCigar(*args, **kwargs) -> parseCigarResult:
    """
    parse cigar string
    """
