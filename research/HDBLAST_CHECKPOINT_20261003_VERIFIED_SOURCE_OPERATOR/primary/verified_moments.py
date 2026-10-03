"""Candidate Arb moment machinery. Only fabricated inputs may run before freeze.

This module does not decode checkpoint arrays. The real-source Taylor callback
requires a caller-owned public-registration authorization callback and is never
called by the fabricated test or benchmark entry points.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
from typing import Callable, Iterable

from flint import acb, acb_poly, arb, arb_series, ctx, fmpq

SOURCE_DEGREE = 24
PHASE_DEGREE = 96
PANELS = 64
PRECISION_BITS = 256
HALF_WIDTH = Fraction(1, 128)
ANALYTIC_RADIUS = Fraction(1, 8)
SOURCE_MAJORANT = Fraction(64)
SOURCE_REMAINDER = Fraction(1, 15 * 2**90)
SOURCE_WORK_MAJORANT = Fraction(32)
SOURCE_WORK_REMAINDER = Fraction(1, 30 * 2**90)
MAX_LOCAL_PHASE = Fraction(4)


def rational_ball(x: Fraction | int) -> arb:
    x = Fraction(x)
    return arb(fmpq(x.numerator, x.denominator))


def inflate_real(x: arb, radius: Fraction) -> arb:
    if radius < 0:
        raise ValueError("A negative remainder is forbidden")
    return x + arb(0, rational_ball(radius).upper())


def inflate_complex(x: acb, radius: Fraction) -> acb:
    # A complex disk |error| <= R is safely enclosed by this square.
    return acb(inflate_real(x.real, radius), inflate_real(x.imag, radius))


def rational_endpoints(x: arb) -> dict[str, str]:
    """Exact rational endpoints; no midpoint decimal is presented as a bound."""
    if not x.is_finite():
        raise ValueError("Nonfinite enclosure")
    return {"lower": str(x.lower().fmpq()), "upper": str(x.upper().fmpq())}


def complex_endpoints(x: acb) -> dict[str, dict[str, str]]:
    return {"real": rational_endpoints(x.real), "imag": rational_endpoints(x.imag)}


def j_monomial(j: int) -> Fraction:
    return Fraction(2, j + 1) if j % 2 == 0 else Fraction(0)


def mixed_monomial(j: int, m: int) -> Fraction:
    """Exactly integral[-1,1] x^j (1-x)^m dx, using y=(1-x)/2."""
    return 2 ** (m + 1) * sum(
        (Fraction(comb(j, ell) * (-2) ** ell, m + ell + 1) for ell in range(j + 1)),
        Fraction(0),
    )


@dataclass(frozen=True)
class PolynomialPanel:
    center: Fraction
    half_width: Fraction
    coefficients: tuple[arb, ...]  # Normalized x=(eta-center)/half_width.
    source_remainder: Fraction

    def __post_init__(self) -> None:
        if self.half_width <= 0 or self.source_remainder < 0:
            raise ValueError("Invalid panel geometry/remainder")
        if not self.coefficients or any(not x.is_finite() for x in self.coefficients):
            raise ValueError("Invalid coefficient balls")


class EntireKernelFamily:
    """Entire exponential and integrated-exponential kernels, without 1/k.

    E_j(z) = integral[-1,1] x^j exp(z(1-x)) dx
    Q_j(z) = integral[-1,1] x^j (exp(z(1-x))-1)/z dx

    The second expression is an entire function, interpreted continuously at
    z=0. We evaluate its entire power series directly at every k, including 0.
    acb_poly evaluates fixed rational coefficient polynomials in C.
    """

    def __init__(self, degree: int = SOURCE_DEGREE, phase_degree: int = PHASE_DEGREE,
                 phase_cap: Fraction = MAX_LOCAL_PHASE):
        if degree < 0 or phase_degree < 0 or phase_cap < 0:
            raise ValueError("Invalid kernel parameters")
        if ctx.prec < PRECISION_BITS:
            raise ValueError("Candidate kernels require at least256-bit Arb precision")
        if phase_cap != MAX_LOCAL_PHASE:
            raise ValueError("Prototype proved exponential majorant only for cap4")
        self.degree = degree
        self.phase_degree = phase_degree
        self.phase_cap = phase_cap
        self.exp_polynomials = []
        self.drift_polynomials = []
        for j in range(degree + 1):
            e = [rational_ball(mixed_monomial(j, m) / factorial(m))
                 for m in range(phase_degree + 1)]
            q = [rational_ball(mixed_monomial(j, m + 1) / factorial(m + 1))
                 for m in range(phase_degree + 1)]
            self.exp_polynomials.append(acb_poly(e))
            self.drift_polynomials.append(acb_poly(q))
        # exp(1)<=49/18<11/4 by its positive series (METHOD.md), and
        # (11/4)^8=214358881/65536<4096.
        self.e_tail = Fraction(2 * 4096 * 8 ** (phase_degree + 1), factorial(phase_degree + 1))
        self.q_tail = Fraction(4 * 4096 * 8 ** (phase_degree + 1), factorial(phase_degree + 2))

    def at(self, k: Fraction | int, half_width: Fraction = HALF_WIDTH) -> tuple[tuple[acb, ...], tuple[acb, ...]]:
        k = Fraction(k)
        if k < 0 or half_width <= 0 or 2 * k * half_width > self.phase_cap:
            raise ValueError("Momentum outside proved local phase cap")
        z = acb(0, rational_ball(2 * k * half_width))
        e = tuple(inflate_complex(p(z), self.e_tail) for p in self.exp_polynomials)
        q = tuple(inflate_complex(p(z), self.q_tail) for p in self.drift_polynomials)
        return e, q


def panel_moments(panel: PolynomialPanel, kernels: tuple[tuple[acb, ...], tuple[acb, ...]]) -> dict[str, acb | arb | Fraction]:
    e, q = kernels
    if len(e) != len(panel.coefficients) or len(q) != len(panel.coefficients):
        raise ValueError("Source/kernel degree mismatch")
    h = rational_ball(panel.half_width)
    m0 = h * sum((a * rational_ball(j_monomial(j)) for j, a in enumerate(panel.coefficients)), arb(0))
    me = h * sum((a * value for a, value in zip(panel.coefficients, e)), acb(0))
    mu = h * h * sum((a * value for a, value in zip(panel.coefficients, q)), acb(0))
    r0 = 2 * panel.half_width * panel.source_remainder
    ru = 2 * panel.half_width ** 2 * panel.source_remainder
    return {
        "M0": inflate_real(m0, r0),
        "Mexp": inflate_complex(me, r0),
        "Mdrift": inflate_complex(mu, ru),
        "source_M0_remainder": r0,
        "source_Mexp_remainder": r0,
        "source_Mdrift_remainder": ru,
    }


def stable_drift(k: Fraction | int, distance: Fraction) -> acb:
    """(exp(2ikd)-1)/(2ik), entire treatment near k=0.

    |phase|<=1 uses a separate degree96 entire series. Else division is by an
    exact nonzero represented rational whose modulus is bounded away from 0
    relative to the requested distance. This is not a 1/k recurrence.
    """
    k = Fraction(k)
    if k < 0 or distance < 0:
        raise ValueError("Invalid drift arguments")
    z = acb(0, rational_ball(2 * k * distance))
    if 2 * k * distance <= 1:
        coefficients = [rational_ball(Fraction(1, factorial(m + 1)))
                        for m in range(PHASE_DEGREE + 1)]
        p = acb_poly(coefficients)
        # exp(1)<3 and tail(phi1)<=3/(N+2)! at |phase|<=1.
        tail = distance * Fraction(3, factorial(PHASE_DEGREE + 2))
        return inflate_complex(rational_ball(distance) * p(z), tail)
    omega = acb(0, rational_ball(2 * k))
    return (z.exp() - 1) / omega


def forced_response(panels: Iterable[PolynomialPanel], k: Fraction | int,
                    incoming_u: acb, incoming_w: acb, epsilon: Fraction,
                    kernel_family: EntireKernelFamily) -> dict[str, acb | arb | Fraction]:
    """Additive global Duhamel sum; no repeated uncertain-state wrapping.

    The caller must authenticate physical input bytes before invoking with real
    panels/modes. This candidate function sees only its explicit objects.
    """
    panels = tuple(panels)
    if not panels:
        raise ValueError("Empty panel collection")
    h = panels[0].half_width
    if any(p.half_width != h for p in panels):
        raise ValueError("This prototype requires a fixed uniform half-width")
    for left, right in zip(panels, panels[1:]):
        if left.center + left.half_width != right.center - right.half_width:
            raise ValueError("Panels do not tile the exact interval")
    a = panels[0].center - h
    b = panels[-1].center + h
    kernels = kernel_family.at(k, h)
    weighted_exp, weighted_drift = acb(0), acb(0)
    source_m0 = arb(0)
    for panel in panels:
        r = panel.center + h
        distance = b - r
        rotation = acb(0, rational_ball(2 * Fraction(k) * distance)).exp()
        moments = panel_moments(panel, kernels)
        weighted_exp += rotation * moments["Mexp"]
        weighted_drift += stable_drift(k, distance) * moments["M0"] + rotation * moments["Mdrift"]
        source_m0 += moments["M0"]
    phase = acb(0, rational_ball(2 * Fraction(k) * (b - a))).exp()
    eps = rational_ball(epsilon)
    return {
        "u": incoming_u + stable_drift(k, b - a) * incoming_w - eps * weighted_drift,
        "w": phase * incoming_w - eps * weighted_exp,
        "M0": source_m0,
        "Mexp": weighted_exp,
        "Mdrift": weighted_drift,
        # Positive source-only bounds, distinct from arithmetic coefficient balls.
        "source_w_remainder": abs(epsilon) * sum((2 * p.half_width * p.source_remainder for p in panels), Fraction(0)),
        "source_u_remainder": abs(epsilon) * sum(
            ((2 * p.half_width * (b - p.center - p.half_width) + 2 * p.half_width ** 2)
             * p.source_remainder for p in panels), Fraction(0)),
    }


def forcing_series_coefficients(eta: arb_series, source_series: arb_series,
                                half_width: Fraction = HALF_WIDTH,
                                degree: int = SOURCE_DEGREE) -> tuple[tuple[arb, ...], tuple[arb, ...]]:
    """Pure jet algebra, callable with fabricated source series before freeze.

    Caller must set ctx.cap>=degree+3 and provide that many source coefficients;
    this function does not construct/evaluate a registered physical source.
    """
    if ctx.cap < degree + 3 or source_series.prec < degree + 3 or eta.prec < degree + 3:
        raise ValueError("Insufficient coefficient/cap precision")
    l = -eta.inv()
    hp = source_series.derivative() / rational_ball(half_width)
    hpp = hp.derivative() / rational_ball(half_width)
    g = 4 * l * l * source_series - 2 * l * hp - hpp
    work = l * g
    if g.prec < degree + 1 or work.prec < degree + 1:
        raise ValueError("Insufficient differentiated series precision")
    return tuple(g[j] for j in range(degree + 1)), tuple(work[j] for j in range(degree + 1))


def future_registered_source_bundle(center: Fraction, source: str,
                                    authorize: Callable[[str, dict], None] | None = None) -> dict[str, PolynomialPanel]:
    """Candidate real-source callback, inaccessible without caller freeze guard.

    The guard is owned by the future complete driver and must check the public
    freeze, source hashes, opaque input pins and every registration condition.
    This prototype does not implement that guard or authorize a physical run.
    """
    if authorize is None:
        raise RuntimeError("No physical source evaluation before public freeze")
    if source not in ("positive_B", "signed_uB"):
        raise ValueError("Unexpected exact source")
    if ctx.prec < PRECISION_BITS:
        raise ValueError("Source coefficients require at least256-bit Arb precision")
    if not Fraction(-9, 2) <= center <= Fraction(-7, 2):
        raise ValueError("Taylor center outside proved domain")
    authorize("real_source_taylor", {"source": source, "center": str(center)})
    prec = SOURCE_DEGREE + 3
    # python-flint clamps series arithmetic by ctx.cap independently of the
    # constructor's prec argument. Set it explicitly and restore the caller's
    # cap; otherwise omitted high coefficients could be misclassified as zero.
    previous_cap = ctx.cap
    ctx.cap = prec
    try:
        eta = arb_series([rational_ball(center), rational_ball(HALF_WIDTH)], prec=prec)
        z = eta + 4
        bump = (1 - (1 - z * z).inv()).exp()
        source_series = bump if source == "positive_B" else z * bump
        coefficients, work_coefficients = forcing_series_coefficients(eta, source_series)
    finally:
        ctx.cap = previous_cap
    return {
        "g": PolynomialPanel(center, HALF_WIDTH, coefficients, SOURCE_REMAINDER),
        "Lg": PolynomialPanel(center, HALF_WIDTH, work_coefficients, SOURCE_WORK_REMAINDER),
    }


def future_registered_source_panel(center: Fraction, source: str,
                                   authorize: Callable[[str, dict], None] | None = None) -> PolynomialPanel:
    return future_registered_source_bundle(center, source, authorize)["g"]
