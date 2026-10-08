#!/usr/bin/env python3
"""Exact finite-positive-weight incoming-state transport; no file/data loader.

Only caller-supplied exact rational rows are accepted. This module has no source,
NPZ, network, callback, phase evaluator, quadrature regeneration or state projector.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from typing import Iterable, Mapping

PI = F(14488038916154245685, 4611686018427387904)
EPSILON = F(3777893186295716171, 37778931862957161709568)
LA, LB = F(2, 9), F(2, 7)
ZERO = F(0)
CAPSULES = tuple(f"{s}/{g}" for s in ("positive_B", "signed_uB") for g in ("coarse", "fine"))
PREFIXES = {"coarse": {64: 2048, 128: 4096, 256: 8192},
            "fine": {64: 4096, 128: 8192, 256: 16384}}


def rational(x: object, name: str = "value") -> F:
    if type(x) is not F:
        raise TypeError(f"{name} must be fractions.Fraction, never float")
    return x


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self) -> None:
        rational(self.lo, "interval lower")
        rational(self.hi, "interval upper")
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(x: F) -> "Interval":
        return Interval(x, x)

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def scale(self, x: F) -> "Interval":
        rational(x)
        if x >= 0:
            return Interval(x*self.lo, x*self.hi)
        return Interval(x*self.hi, x*self.lo)

    def shift(self, x: F) -> "Interval":
        return self + Interval.point(x)

    def intersect(self, other: "Interval") -> "Interval":
        return Interval(max(self.lo, other.lo), min(self.hi, other.hi))

    def absmax(self) -> F:
        return max(abs(self.lo), abs(self.hi))

    def contains(self, x: F) -> bool:
        return self.lo <= x <= self.hi


@dataclass(frozen=True)
class NodeError:
    k: F
    weight: F
    # d=2*k*c=(2*k*Re(u_1)-Im(w_1))/epsilon; target c=0 proved externally.
    d: F
    delta_u_re: Interval
    delta_u_im: Interval
    delta_w_re: Interval
    delta_w_im: Interval

    def __post_init__(self) -> None:
        for name in ("k", "weight", "d"):
            rational(getattr(self, name), name)
        if self.k <= 0 or self.weight <= 0:
            raise ValueError("k and inherited weight must be strictly positive")
        for name in ("delta_u_re", "delta_u_im", "delta_w_re", "delta_w_im"):
            if type(getattr(self, name)) is not Interval:
                raise TypeError(f"{name} must be Interval")
        # This is a consistency check, not a proof of the caller's enclosures.
        # It never alters/project/replaces the exact saved canonical residual.
        self.delta_w_im.shift(self.d).scale(F(1)/(2*self.k)).intersect(self.delta_u_re)

    @property
    def c(self) -> F:
        return self.d/(2*self.k)


def from_saved_target(k: F, weight: F, u_saved: tuple[F, F],
                      w_saved: tuple[F, F], target: Mapping,
                      exponent: int = -96) -> NodeError:
    """Bridge to entire_target_engine integer endpoint output.

    u_saved,w_saved are ALREADY normalized exact U=u_1/epsilon,W=w_1/epsilon.
    target is {U:{real:[int,int],imag:[int,int]}, W:{...}} on 2**exponent.
    No assumption about target center normalization replaces its analytic proof.
    """
    if type(exponent) is not int or exponent > 0 or exponent < -16384:
        raise ValueError("unsupported target endpoint exponent")
    scale = F(1, 1 << -exponent)
    def difference(saved: F, part: str, axis: str) -> Interval:
        rational(saved)
        endpoints = target[part][axis]
        if len(endpoints) != 2 or any(type(x) is not int for x in endpoints):
            raise TypeError("target endpoints must be two exact integers")
        lo, hi = endpoints
        if lo > hi:
            raise ValueError("reversed target endpoints")
        return Interval(saved-hi*scale, saved-lo*scale)
    return NodeError(k, weight, 2*k*u_saved[0]-w_saved[1],
                     difference(u_saved[0], "U", "real"),
                     difference(u_saved[1], "U", "imag"),
                     difference(w_saved[0], "W", "real"),
                     difference(w_saved[1], "W", "imag"))


def node_terms(n: NodeError) -> dict:
    """Fused exact terms: no reciprocal-k denominator enters a weighted sum."""
    k, k2, d = n.k, n.k*n.k, n.d
    q = n.weight/(4*PI*PI)
    mu = 2*q*k2
    dw = n.delta_w_re.absmax()+n.delta_w_im.absmax()
    du = n.delta_u_re.absmax()+n.delta_u_im.absmax()
    ac = abs(d)
    rca = q*(k2+3*LA*LA/2)*d
    rcb = q*(k2+3*LB*LB/2)*d
    pca = q*(k2/3-LA*LA/2)*d
    pcb = q*(k2/3-LB*LB/2)*d
    rwa = n.delta_w_im.scale(q*3*LA*LA/2) + n.delta_w_re.scale(-q*LA*k)
    pwa = n.delta_w_im.scale(-q*(2*k2/3+LA*LA/2)) + n.delta_w_re.scale(-q*LA*k)
    ra = q*(LA*k+3*LA*LA/2)*dw
    rb = q*(LB*k+3*LB*LB/2)*dw
    workc = q*3*(LB*LB-LA*LA)/2*d
    timec = q*(k2/3-(LB-LA)/2)*d
    return {
        "weight_sum": n.weight,
        "mu_sum": mu,
        "weight_k1_sum": n.weight*k,
        "weight_k2_sum": n.weight*k2,
        "weight_k3_sum": n.weight*k2*k,
        "weight_k4_sum": n.weight*k2*k2,
        "anchor_R_c": rca, "anchor_P_c": pca,
        "endpoint_R_c": rcb, "endpoint_P_c": pcb,
        "anchor_R_W": rwa, "anchor_P_W": pwa,
        "R_c_abs_sum": q*(k2+3*LB*LB/2)*ac,
        "R_A_abs_sum": rb,
        "P_c_abs_sum": q*(k2/3+LB*LB/2)*ac,
        "P_A_abs_sum": q*(2*k2/3+5*LB*LB/4)*dw,
        "work_c_signed": workc,
        "work_c_abs_sum": q*3*(LB*LB-LA*LA)/2*ac,
        "work_A_abs_sum": ra+rb,
        "int_pressure_c_signed": timec,
        "int_pressure_c_abs_sum": q*(k2/3+(LB-LA)/2)*ac,
        "int_pressure_A_abs_sum": q*(2*k/3+(LA+LB)/2)*dw,
        "delta_U_weighted_L1_upper": mu*du,
        "delta_W_weighted_L1_upper": mu*dw,
        "kA_weighted_L1_upper": mu*dw/2,
        "canonical_d_weighted_q_L1": q*ac,
        "canonical_d_weighted_qk2_L1": q*k2*ac,
    }


class PrefixAccumulator:
    """One monotonically ordered capsule, consumed once, with copied snapshots."""
    def __init__(self) -> None:
        self.count = 0
        self.previous_k = ZERO
        self.sums: dict = {}
        self.maxima = {key: ZERO for key in (
            "delta_U_L1_upper", "delta_W_L1_upper", "c_abs", "kA_L1_upper",
            "first_order_Wronskian_defect_abs")}

    def add(self, node: NodeError) -> None:
        if type(node) is not NodeError:
            raise TypeError("expected NodeError")
        if node.k <= self.previous_k:
            raise ValueError("k nodes must be strictly increasing, with no duplicates")
        terms = node_terms(node)
        for key, val in terms.items():
            if key not in self.sums:
                self.sums[key] = val
            else:
                self.sums[key] = self.sums[key]+val
        cabs = abs(node.c)
        ul1 = node.delta_u_re.absmax()+node.delta_u_im.absmax()
        wl1 = node.delta_w_re.absmax()+node.delta_w_im.absmax()
        values = {"delta_U_L1_upper": ul1, "delta_W_L1_upper": wl1,
                  "c_abs": cabs, "kA_L1_upper": wl1/2,
                  "first_order_Wronskian_defect_abs": 2*EPSILON*cabs}
        for key, val in values.items():
            self.maxima[key] = max(self.maxima[key], val)
        self.count += 1
        self.previous_k = node.k

    def snapshot(self, capsule: str, cutoff: F) -> dict:
        if self.count == 0:
            raise ValueError("empty prefix")
        s = dict(self.sums)
        s["anchor_R_interval"] = s["anchor_R_W"].shift(s["anchor_R_c"])
        s["anchor_P_interval"] = s["anchor_P_W"].shift(s["anchor_P_c"])
        # Canonical aggregate is affine in L^2: its absolute maximum on the
        # full time interval occurs at one of the endpoints, exactly.
        s["R_uniform_abs_upper"] = max(abs(s["anchor_R_c"]), abs(s["endpoint_R_c"]))+s["R_A_abs_sum"]
        s["P_uniform_abs_upper"] = max(abs(s["anchor_P_c"]), abs(s["endpoint_P_c"]))+s["P_A_abs_sum"]
        s["work_abs_upper"] = abs(s["work_c_signed"])+s["work_A_abs_sum"]
        s["int_pressure_abs_upper"] = abs(s["int_pressure_c_signed"])+s["int_pressure_A_abs_sum"]
        s["integral_absolute_pressure_abs_upper"] = s["P_uniform_abs_upper"]  # duration=1
        return {"case": f"{capsule}/{cutoff}", "capsule": capsule,
                "K": cutoff, "node_count": self.count, "last_k": self.previous_k,
                "sums": s, "suprema": dict(self.maxima)}


def aggregate_capsule(capsule: str, nodes: Iterable[NodeError],
                      expected_count: int, prefixes: Mapping[int, int]) -> list[dict]:
    if not isinstance(capsule, str) or not capsule:
        raise ValueError("capsule identity required")
    if type(expected_count) is not int or expected_count <= 0:
        raise ValueError("positive exact expected node count required")
    plan = list(prefixes.items())
    if (not plan or any(type(k) is not int or type(n) is not int or k <= 0 or n <= 0 for k, n in plan)
            or plan != sorted(plan) or any(plan[i][1] >= plan[i+1][1] for i in range(len(plan)-1))
            or plan[-1][1] != expected_count):
        raise ValueError("ordered increasing cutoff/count prefixes ending at full count required")
    acc, results = PrefixAccumulator(), []
    bycount = {count: cutoff for cutoff, count in plan}
    for index, node in enumerate(nodes):
        if index >= expected_count:
            raise ValueError("too many capsule nodes")
        # Check the positional prefix definition against every exact k. No
        # value is silently discarded or replaced by a regenerated node.
        for cutoff, count in plan:
            if (index < count and node.k >= cutoff) or (index >= count and node.k <= cutoff):
                raise ValueError("positional prefix and exact k cutoff mask disagree")
        acc.add(node)
        if acc.count in bycount:
            results.append(acc.snapshot(capsule, F(bycount[acc.count])))
    if acc.count != expected_count:
        raise ValueError("incomplete capsule")
    return results


def aggregate_registered(capsules: Mapping[str, Iterable[NodeError]]) -> dict:
    if set(capsules) != set(CAPSULES):
        raise ValueError("all four registered capsules required separately")
    rows = []
    for capsule in CAPSULES:
        plan = PREFIXES[capsule.rsplit("/", 1)[1]]
        rows.extend(aggregate_capsule(capsule, capsules[capsule], plan[256], plan))
    return {"scope": "conditional finite weighted discrete incoming-state difference",
            "rows": rows, "row_count": len(rows), "distinct_capsule_nodes": 49152,
            "prefix_node_appearances": 86016,
            "hypotheses": ["input rectangles enclose saved-minus-BD target at incoming anchor",
                           "target canonical c=0 independently proved",
                           "identical subsequent real forcing and complete contacts",
                           "fixed direct linear operators; eta in [-9/2,-7/2]"],
            "full_continuous_integrals": "UNRESOLVED",
            "full_twelve_case_pressure_contact_gate_2e_minus_8": "UNRESOLVED",
            "contact_source_time_quadrature_UV": "UNRESOLVED",
            "historical_metric_calibration": "FAIL_UNCHANGED",
            "higher_dimensional_big_bang_cause": "NOT_ESTABLISHED",
            "external_mathematical_novelty": "NOT_ASSESSED"}


def encode_exact(value):
    """Exact JSON-friendly bounded summary; no per-node output is created."""
    if isinstance(value, F):
        return {"numerator": str(value.numerator), "denominator": str(value.denominator)}
    if isinstance(value, Interval):
        return {"lo": encode_exact(value.lo), "hi": encode_exact(value.hi)}
    if isinstance(value, dict):
        return {k: encode_exact(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode_exact(v) for v in value]
    return value
