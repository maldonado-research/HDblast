"""PROPOSED_UNVERIFIED_FOLLOWON: exact folded-cone geometry only.

No source, field, waveform, sign certificate, or physics oracle is evaluated.
Coordinates and physical areas are in (s,u), so no coordinate Jacobian is
required. The separate wave-kernel weights are supplied as exact constants.
"""

from dataclasses import dataclass
from fractions import Fraction


ZERO = Fraction(0)
HALF = Fraction(1, 2)
ONE = Fraction(1)
Point = tuple[Fraction, Fraction]
Polygon = tuple[Point, ...]


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("coordinates and mesh_step must be Fraction or int")
    return Fraction(value)


def _clean(vertices):
    result = []
    for vertex in vertices:
        if not result or vertex != result[-1]:
            result.append(vertex)
    if len(result) > 1 and result[0] == result[-1]:
        result.pop()
    return tuple(result)


def _clip(polygon, normal_s, normal_u, bound):
    """Intersect a convex polygon with normal_s*s+normal_u*u <= bound."""
    if not polygon:
        return ()
    out = []
    previous = polygon[-1]
    previous_gap = bound - normal_s * previous[0] - normal_u * previous[1]
    for current in polygon:
        current_gap = bound - normal_s * current[0] - normal_u * current[1]
        if (previous_gap >= 0) != (current_gap >= 0):
            t = previous_gap / (previous_gap - current_gap)
            out.append((previous[0] + t * (current[0] - previous[0]),
                        previous[1] + t * (current[1] - previous[1])))
        if current_gap >= 0:
            out.append(current)
        previous, previous_gap = current, current_gap
    return _clean(out)


def _measure(polygon, first_moments):
    """Shoelace area and optional integrals of s and u over the polygon."""
    twice_area = ZERO
    moment_s_times_six = ZERO
    moment_u_times_six = ZERO
    for p, q in zip(polygon, polygon[1:] + polygon[:1]):
        cross = p[0] * q[1] - q[0] * p[1]
        twice_area += cross
        if first_moments:
            moment_s_times_six += (p[0] + q[0]) * cross
            moment_u_times_six += (p[1] + q[1]) * cross
    if twice_area < 0:
        raise ValueError("internal polygon orientation is clockwise")
    moments = ((moment_s_times_six / 6, moment_u_times_six / 6)
               if first_moments else None)
    return twice_area / 2, moments


@dataclass(frozen=True)
class Cell:
    mesh_index: tuple[int, int]
    region: str
    polygon: Polygon
    area: Fraction
    bbox: tuple[Fraction, Fraction, Fraction, Fraction]
    first_moments: tuple[Fraction, Fraction] | None
    bbox_exterior_error_eligible: bool

    @property
    def lambda_wyy_weight(self):
        return ONE if self.region == "Dminus" else HALF

    @property
    def all_dirichlet_weight(self):
        return ZERO if self.region == "Dminus" else HALF


@dataclass(frozen=True)
class FoldedGeometry:
    a: Fraction
    b: Fraction
    mesh_step: Fraction
    dminus: tuple[Cell, ...]
    strip: tuple[Cell, ...]
    boundary_interval: tuple[Fraction, Fraction]
    boundary_cells: tuple[tuple[Fraction, Fraction], ...]


def _cell(polygon, mesh_index, region, first_moments):
    if len(polygon) < 3:
        return None
    area, moments = _measure(polygon, first_moments)
    if area == 0:
        return None
    ss, us = zip(*polygon)
    bbox = (min(ss), max(ss), min(us), max(us))
    return Cell(mesh_index, region, polygon, area, bbox, moments,
                bbox[2] >= bbox[1])


def _ceil_ratio(numerator, denominator):
    ratio = numerator / denominator
    return -(-ratio.numerator // ratio.denominator)


def build_geometry(a, b, mesh_step, *, first_moments=True):
    """Return positive-area cells in Dminus and Dplus minus Dminus separately.

    All inputs must be Fraction or int; a,b>=0 and mesh_step>0. The uniform
    square mesh starts at (0,0). Half-plane clipping uses exact Fractions.
    Shared edges have zero area; the two returned cell lists are disjoint up
    to such edges. A zero-time cone returns no area cells. Dminus is empty
    for a<=b. Boundary data use the separately returned interval [0,max(a-b,0)]
    and positive-length one-dimensional mesh intervals.

    Folded vector weights are I/2 on the strip and P_N on Dminus. For lambda
    W_yy (a Robin row), the scalar weights are 1/2 and 1. For an all-Dirichlet
    vector they are 1/2 and 0. Weights never modify the geometric area.

    bbox_exterior_error_eligible means exactly u_lo>=s_hi on the WHOLE bbox.
    Polygon membership alone is insufficient for that error criterion. This
    generic geometry routine makes no oracle-domain claim for arbitrary a,b.
    """
    a, b, mesh_step = map(_rational, (a, b, mesh_step))
    if a < 0 or b < 0:
        raise ValueError("a and b must be nonnegative")
    if mesh_step <= 0:
        raise ValueError("mesh_step must be positive")
    d = a - b
    minus = []
    strip = []
    for i in range(_ceil_ratio(a, mesh_step)):
        s0, s1 = i * mesh_step, min((i + 1) * mesh_step, a)
        for j in range(_ceil_ratio(a + b, mesh_step)):
            u0, u1 = j * mesh_step, min((j + 1) * mesh_step, a + b)
            polygon = ((s0, u0), (s1, u0), (s1, u1), (s0, u1))
            # Mesh rectangles already impose s>=0 and u>=0.
            polygon = _clip(polygon, ONE, ONE, a + b)
            polygon = _clip(polygon, ONE, -ONE, d)
            if d > 0:
                overlap = _cell(_clip(polygon, ONE, ONE, d), (i, j),
                                "Dminus", first_moments)
                remainder = _cell(_clip(polygon, -ONE, -ONE, -d), (i, j),
                                  "strip", first_moments)
                if overlap is not None:
                    minus.append(overlap)
            else:
                remainder = _cell(polygon, (i, j), "strip", first_moments)
            if remainder is not None:
                strip.append(remainder)
    boundary_end = max(d, ZERO)
    boundary_cells = tuple((i * mesh_step,
                            min((i + 1) * mesh_step, boundary_end))
                           for i in range(_ceil_ratio(boundary_end, mesh_step)))
    return FoldedGeometry(a, b, mesh_step, tuple(minus), tuple(strip),
                          (ZERO, boundary_end), boundary_cells)
