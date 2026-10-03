#!/usr/bin/env python3
"""Adversarial fabricated-only checks; never calls the physical source callback."""
from fractions import Fraction as F
from pathlib import Path
import sys

from flint import acb, arb, arb_series, ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verified_moments as vm


def expect_contains(label, actual, expected):
    if not actual.contains(expected):
        raise AssertionError(f"{label}: candidate ball does not contain reference ball")


def fabricated_jet_cap_check():
    """Use h(eta)=2+3eta+5eta^2 at c=2, outside the physical source domain."""
    ctx.prec = 256
    c, width = F(2), F(1, 128)
    ncoeff = 25

    def fake_g(cap):
        ctx.cap = cap
        eta = arb_series([vm.rational_ball(c), vm.rational_ball(width)], prec=ncoeff + 2)
        source = 2 + 3 * eta + 5 * eta * eta
        hp = source.derivative() / vm.rational_ball(width)
        hpp = hp.derivative() / vm.rational_ball(width)
        lapse = -eta.inv()
        return 4 * lapse * lapse * source - 2 * lapse * hp - hpp

    hs = (2 + 3 * c + 5 * c * c, width * (3 + 10 * c), 5 * width * width)
    hps = (3 + 10 * c, 10 * width)
    expected = []
    for n in range(ncoeff):
        l = lambda j: (-1) ** (j + 1) * width**j / c ** (j + 1)
        m = lambda j: (j + 1) * (-1) ** j * width**j / c ** (j + 2)
        value = (
            4 * sum((m(n - j) * hs[j] for j in range(min(n, 2) + 1)), F(0))
            - 2 * sum((l(n - j) * hps[j] for j in range(min(n, 1) + 1)), F(0))
            - (10 if n == 0 else 0)
        )
        expected.append(value)

    unsafe = fake_g(10)
    assert unsafe.prec < ncoeff
    assert expected[-1] != 0
    assert unsafe[ncoeff - 1] == 0
    assert not unsafe[ncoeff - 1].contains(vm.rational_ball(expected[-1]))
    safe = fake_g(ncoeff + 2)
    assert safe.prec >= ncoeff
    for n, value in enumerate(expected):
        expect_contains(f"normalized source derivative coefficient {n}", safe[n], vm.rational_ball(value))
    eta = arb_series([vm.rational_ball(c), vm.rational_ball(width)], prec=ncoeff + 2)
    source = 2 + 3 * eta + 5 * eta * eta
    computed_g, computed_work = vm.forcing_series_coefficients(eta, source, width, ncoeff - 1)
    for n, value in enumerate(expected):
        expect_contains(f"pure helper g coefficient {n}", computed_g[n], vm.rational_ball(value))
        expected_work = sum(
            ((-1) ** (j + 1) * width**j / c ** (j + 1) * expected[n - j] for j in range(n + 1)), F(0)
        )
        expect_contains(f"pure helper Lg coefficient {n}", computed_work[n], vm.rational_ball(expected_work))
    ctx.cap = 10
    try:
        vm.forcing_series_coefficients(eta, source, width, ncoeff - 1)
    except ValueError:
        pass
    else:
        raise AssertionError("Pure jet helper failed to reject insufficient global cap")
    ctx.cap = ncoeff + 2
    print(f"SERIES_CAP_HAZARD_REPRODUCED unsafe_prec={unsafe.prec}; normalized fake jet passes at cap={ncoeff + 2}")


def closed_polynomial_reference(k, epsilon, source_offset=F(0)):
    """Independent closed polynomial integrals, at much higher precision."""
    ctx.prec = 1024
    initial_u = acb(vm.rational_ball(F(1)), vm.rational_ball(F(1, 3)))
    initial_w = acb(vm.rational_ball(F(2)), vm.rational_ball(F(-1, 5)))
    coefficients = (F(2) + source_offset, F(3), F(5))
    source0 = sum((c / (j + 1) for j, c in enumerate(coefficients)), F(0))
    omega = acb(0, vm.rational_ball(2 * k))
    if k == 0:
        phase = acb(1)
        drift = acb(1)
        f = [acb(vm.rational_ball(F(1, j + 1))) for j in range(3)]
        d = [acb(vm.rational_ball(F(1, (j + 1) * (j + 2)))) for j in range(3)]
    else:
        phase = omega.exp()
        drift = (phase - 1) / omega
        f = [drift]
        for j in (1, 2):
            f.append((j * f[-1] - 1) / omega)
        d = [(value - vm.rational_ball(F(1, j + 1))) / omega for j, value in enumerate(f)]
    weighted_f = sum((vm.rational_ball(c) * value for c, value in zip(coefficients, f)), acb(0))
    weighted_d = sum((vm.rational_ball(c) * value for c, value in zip(coefficients, d)), acb(0))
    eps = vm.rational_ball(epsilon)
    return {
        "u": initial_u + drift * initial_w - eps * weighted_d,
        "w": phase * initial_w - eps * weighted_f,
        "M0": vm.rational_ball(source0),
        "Mexp": weighted_f,
        "Mdrift": weighted_d,
    }


def fabricated_response_checks():
    """All sources below are manufactured quadratics on eta in [0,1]."""
    ctx.prec = 256
    family = vm.EntireKernelFamily(degree=2)
    epsilon = F(7, 11)
    h = F(1, 128)
    frequencies = (F(0), F(1, 10**30), F(1, 10000), F(17, 3), F(256))
    count = 0
    for remainder in (F(0), F(1, 10000)):
        panels = []
        for j in range(64):
            c = F(2 * j + 1, 128)
            coefficients = (2 + 3 * c + 5 * c * c, h * (3 + 10 * c), 5 * h * h)
            # The coefficient radii test coverage of uncertain fabricated jets.
            balls = tuple(vm.inflate_real(vm.rational_ball(a), F(1, 10**80)) for a in coefficients)
            panels.append(vm.PolynomialPanel(c, h, balls, remainder))
        for k in frequencies:
            ctx.prec = 256
            incoming_u = acb(vm.rational_ball(F(1)), vm.rational_ball(F(1, 3)))
            incoming_w = acb(vm.rational_ball(F(2)), vm.rational_ball(F(-1, 5)))
            actual = vm.forced_response(panels, k, incoming_u, incoming_w, epsilon, family)
            assert actual["source_w_remainder"] == abs(epsilon) * remainder
            assert actual["source_u_remainder"] == abs(epsilon) * remainder / 2
            reference = closed_polynomial_reference(k, epsilon, source_offset=remainder)
            for key in ("u", "w", "M0", "Mexp", "Mdrift"):
                expect_contains(f"{key}, k={k}, remainder={remainder}", actual[key], reference[key])
            for key in ("u", "w", "Mexp", "Mdrift"):
                exported = vm.complex_endpoints(actual[key])
                for component in ("real", "imag"):
                    lower, upper = exported[component]["lower"], exported[component]["upper"]
                    assert F(lower) <= F(upper)
                    expected_component = getattr(reference[key], component)
                    assert vm.rational_ball(F(lower)) <= expected_component.lower()
                    assert vm.rational_ball(F(upper)) >= expected_component.upper()
            count += 1
    print(f"FABRICATED_GLOBAL_DUHAMEL_REFERENCE_CASES_PASS count={count}")


if __name__ == "__main__":
    fabricated_jet_cap_check()
    fabricated_response_checks()
    print("FABRICATED_ONLY_AUDIT_CHECKS_PASS; physical source callback never invoked")
