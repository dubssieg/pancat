"Tools to represent GFA format"
from enum import Enum
from re import sub


class Orientation(Enum):
    "Describes the way a node is read"
    FORWARD = '+'
    REVERSE = '-'


class GfaStyle(Enum):
    "Describes the different possible formats"
    RGFA = 'rGFA'
    GFA1 = 'GFA1'
    GFA1_1 = 'GFA1.1'
    GFA1_2 = 'GFA1.2'
    GFA2 = 'GFA2'


class LineType(Enum):
    "Modelizes the line type in GFA format by the meaning of first char of sequence"
    HEADER = 'H'
    SEGMENT = 'S'
    LINE = 'L'
    CONTAINMENT = 'C'
    PATH = 'P'
    WALK = 'W'
    JUMP = 'J'

    def __repr__(self) -> str:
        return self.name

    def __str__(self) -> str:
        return self.value


class Record():
    """Modelizes the data within a GFA-like line.

    Raises:
        ValueError: if line is not intended for targeted format
    """

    class Segment():
        """
        Modelizes a small sequence.
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            self.name = sub('\D', '', datas[1])
            # self.seq = datas[2]
            self.length = len(datas[2])
            if gfa_style == GfaStyle.RGFA:
                self.origin = int(datas[6][5:])

    class Line():
        """
        Modelizes a link between two segments
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style == GfaStyle.RGFA:
                self.origin = int(datas[6][5:])
            self.start = sub('\D', '', datas[1])
            self.end = sub('\D', '', datas[3])
            self.orientation = f"{datas[2]}/{datas[4]}"

    class Containment():
        """
        Modelizes a line overlap
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style == GfaStyle.RGFA:
                raise ValueError(
                    f"Incompatible version format, C-lines vere added in GFA1 and were absent from {gfa_style}.")

    class Path():
        """
        Modelizes a path within the graph
        Path is of shape [ ( (a,+) , (b,+) ) ... ]
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style == GfaStyle.RGFA:
                raise ValueError(
                    f"Incompatible version format, P-lines vere added in GFA1 and were absent from {gfa_style}.")
            self.name = datas[1]
            self.path = [
                (
                    node[:-1],
                    Orientation(node[-1])
                )
                for node in datas[2].split(',')
            ]

    class Header():
        """
        Head of file to put version in it
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style == GfaStyle.RGFA:
                raise ValueError(
                    f"Incompatible version format, H-lines vere added in GFA1 and were absent from {gfa_style}.")
            self.version = datas[1][5:]

    class Jump():
        """
        J-lines shows jumps within paths/walks
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style in [GfaStyle.RGFA, GfaStyle.GFA1, GfaStyle.GFA1_1]:
                raise ValueError(
                    f"Incompatible version format, J-lines vere added in GFA1.2 and were absent from {gfa_style}.")
            raise NotImplementedError

    class Walk():
        """
        W-lines describes paths within the graph.
        """

        def __init__(self, datas: list, gfa_style: GfaStyle) -> None:
            if gfa_style in [GfaStyle.RGFA, GfaStyle.GFA1]:
                raise ValueError(
                    f"Incompatible version format, W-lines vere added in GFA1.1 and were absent from {gfa_style}.")
            self.idf = datas[1]
            self.origin = int(datas[2])
            self.name = datas[3]
            self.target = datas[4]
            self.length = datas[5]
            self.path = [
                (
                    node[1:],
                    Orientation(node[0])
                )
                for node in datas[6].replace('>', ',+').replace('<', ',-')[1:].split(',')
            ]

    def __init__(self, line: str, gfa_type: str) -> None:
        datas: list = line.split()
        self.gfastyle: GfaStyle = GfaStyle(gfa_type)
        self.linetype: LineType = LineType(line[0])
        if self.linetype == LineType.SEGMENT:
            self.line = self.Segment(datas, self.gfastyle)
        elif self.linetype == LineType.LINE:
            self.line = self.Line(datas, self.gfastyle)
        elif self.linetype == LineType.CONTAINMENT:
            self.line = self.Containment(datas, self.gfastyle)
        elif self.linetype == LineType.PATH:
            self.line = self.Path(datas, self.gfastyle)
        elif self.linetype == LineType.HEADER:
            self.line = self.Header(datas, self.gfastyle)
        elif self.linetype == LineType.WALK:
            self.line = self.Walk(datas, self.gfastyle)
        elif self.linetype == LineType.JUMP:
            self.line = self.Jump(datas, self.gfastyle)
        else:
            raise ValueError('Unknown line for GFA standard')

    def w_to_p_line(self) -> None:
        if not isinstance(self.line, self.Walk):
            raise ValueError(
                f'Bad call of method on {type(self.line)} input type')
        datas = ['P', self.line.name, ','.join(self.line.walk)]
        self.linetype: LineType = LineType(datas[0])
        self.line = self.Path(datas, self.gfastyle)
