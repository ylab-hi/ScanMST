"""BAM file parser"""
import typing

__all__ = [
    "parseCigar",
    "parseCigarResult"
]


class parseCigarResult():
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


def parseCigar(arg0: str) -> parseCigarResult:
    """
    parse cigar string
    """
