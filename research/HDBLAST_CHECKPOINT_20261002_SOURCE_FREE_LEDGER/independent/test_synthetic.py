"""Only fabricated modes and histories. Never opens any retained physical input."""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import numpy as np
from exact_binary import exact_binary, exact_complex, require
from ledger_core import (context, convert, decimal, calculate_profiles,
                         simpson_prefix, simpson_diagnostics, EPS_RATIO, PI_RATIO)

LD, CD = np.longdouble, np.clongdouble


def rejects(function, message):
    try:
        function()
    except (RuntimeError, ValueError, TypeError, OverflowError):
        return
    raise RuntimeError('Mutation accepted: ' + message)


def run():
    checks = []
    def record(name):
        checks.append(name)
    for value in (LD(0), LD('-0'), LD(1) + np.ldexp(LD(1), -60),
                  LD(-1) - np.ldexp(LD(1), -60), np.nextafter(LD(0), LD(1)),
                  np.finfo(LD).max, -np.finfo(LD).max):
        exact_binary(value)
    record('exact_binary_signed_zero_extended_increment_extreme_roundtrips')
    value = CD(LD(1) + np.ldexp(LD(1), -60)) + CD(1j) * (LD(2) - np.ldexp(LD(1), -60))
    real, imag = exact_complex(value)
    require(real.longdouble() == value.real and imag.longdouble() == value.imag,
            'Complex precision witness failed')
    record('complex_extended_precision_roundtrip')
    rejects(lambda: exact_binary(float(1)), 'Python binary64 float')
    rejects(lambda: exact_binary(np.float64(1)), 'NumPy binary64 scalar')
    rejects(lambda: exact_complex(np.complex128(1 + 2j)), 'NumPy complex binary64')
    record('precoercion_dtype_downcast_mutations_rejected')
    for value in (LD('nan'), LD('inf'), LD('-inf')):
        rejects(lambda value=value: exact_binary(value), 'Nonfinite LD')
    rejects(lambda: exact_complex(CD(LD('inf'))), 'Nonfinite complex component')
    record('nonfinite_components_rejected')
    require(LD('0.0001').as_integer_ratio() == EPS_RATIO, 'Pinned producer epsilon ratio differs')
    require(np.arccos(LD(-1)).as_integer_ratio() == PI_RATIO, 'Pinned producer pi ratio differs')
    record('exact_retained_producer_constant_ratios')
    nodes = np.array(['1.25', '7.5', '255.75'], dtype=LD)
    weights = np.array(['0.125', '0.25', '0.0625'], dtype=LD)
    ua = np.array([CD('0.00001') + CD(1j) * LD('0.00002'),
                   CD('-0.00002') + CD(1j) * LD('0.00003'),
                   CD('0.00004') - CD(1j) * LD('0.00001')], dtype=CD)
    wa = np.array([CD('0.00003') - CD(1j) * LD('0.00001'),
                   CD('-0.00001') + CD(1j) * LD('0.00002'),
                   CD('0.00002') + CD(1j) * LD('0.00001')], dtype=CD)
    for dps in (80, 100):
        ctx = context(dps)
        eps, pi = ctx.mpf(EPS_RATIO[0]) / EPS_RATIO[1], ctx.mpf(PI_RATIO[0]) / PI_RATIO[1]
        # Fabricated data intentionally includes Re(c)!=0; both65/129 test refresh boundaries.
        for samples in (65, 129):
            times = LD('-2.5') + np.arange(samples, dtype=LD) / LD(samples - 1)
            ub, wb = [], []
            for k0, u0, w0 in zip(nodes, ua, wa):
                k, u, w = convert(ctx, k0), ctx.mpc(convert(ctx, u0.real), convert(ctx, u0.imag)), ctx.mpc(convert(ctx, w0.real), convert(ctx, w0.imag))
                phase = ctx.exp(ctx.j * 2 * k)
                d = w / (ctx.j * 2 * k)
                final_u, final_w = u + (phase - 1) * d, phase * w
                # Retained synthetic endpoint components round only to the genuine LD dtype.
                ub.append(CD(LD(ctx.nstr(final_u.real, 110))) + CD(1j) * LD(ctx.nstr(final_u.imag, 110)))
                wb.append(CD(LD(ctx.nstr(final_w.real, 110))) + CD(1j) * LD(ctx.nstr(final_w.imag, 110)))
            _, ts, result = calculate_profiles(nodes, weights, ua, wa, np.array(ub, dtype=CD),
                                               np.array(wb, dtype=CD), times, (1, 2, 3), dps)
            for cutoff, profile in enumerate(result, 1):
                for j, eta in enumerate(ts):
                    R = ctx.mpf(0)
                    P = ctx.mpf(0)
                    L = -1 / eta
                    for k0, wt, u0, w0 in zip(nodes[:cutoff], weights[:cutoff], ua[:cutoff], wa[:cutoff]):
                        k, weight = convert(ctx, k0), convert(ctx, wt)
                        u_init = ctx.mpc(convert(ctx, u0.real), convert(ctx, u0.imag))
                        w_init = ctx.mpc(convert(ctx, w0.real), convert(ctx, w0.imag))
                        E = ctx.exp(ctx.j * 2 * k * (eta - ts[0]))
                        u = u_init + (E - 1) * w_init / (ctx.j * 2 * k)
                        w = E * w_init
                        v = ctx.exp(-ctx.j * k * eta) / ctx.sqrt(2 * k)
                        D = (-ctx.j * k - L) * v
                        dv = v * u
                        dD = v * (w + (-ctx.j * k - L) * u)
                        measure = weight * k * k / (2 * pi * pi)
                        kinetic = ctx.re(ctx.conj(D) * dD)
                        field = ctx.re(ctx.conj(v) * dv)
                        R += measure * (kinetic + (k * k + 2 * L * L) * field) / eps
                        P += measure * (kinetic - (k * k / 3 + 2 * L * L) * field) / eps
                    tolerance = ctx.mpf(10) ** (-(dps - 12))
                    require(abs(R - profile['R'][j]) < tolerance and abs(P - profile['P'][j]) < tolerance,
                            'Recurrence profile disagrees with raw complex bilinears')
                    require(abs(profile['F'][j] - L * (profile['R'][j] - 3 * profile['P'][j])) < tolerance,
                            'Independent F inventory mismatch')
                require(abs(profile['I_ab'] - (profile['R'][-1] - profile['R'][0])) < tolerance,
                        'Independent primitive/free stress increment mismatch')
                require(abs(profile['E_flow']) <= profile['E_flow_triangle'] + tolerance,
                        'Signed flow exceeds its mathematical triangle')
            require(abs(result[-1]['R'][0]) > ctx.mpf('0.1'), 'Nonzero Re(c) precision witness ineffective')
        record('raw_bilinear_profiles_fixed_refresh_nonzero_Re_c_' + str(dps))
    # Independent high-precision quadrature checks F integration on fabricated low-frequency nodes.
    ctx = context(100)
    k0, wt = LD('0.75'), LD('0.125')
    lowtimes = np.array(['-2.5', '-2', '-1.5'], dtype=LD)
    _, ts, profiles = calculate_profiles(np.array([k0], dtype=LD), np.array([wt], dtype=LD),
        ua[:1], wa[:1], ua[:1], wa[:1], lowtimes, (1,), 100)
    k = convert(ctx, k0)
    u0 = ctx.mpc(convert(ctx, ua[0].real), convert(ctx, ua[0].imag))
    w0 = ctx.mpc(convert(ctx, wa[0].real), convert(ctx, wa[0].imag))
    d, c = w0 / (ctx.j * 2 * k), u0 - w0 / (ctx.j * 2 * k)
    eps, pi = ctx.mpf(EPS_RATIO[0]) / EPS_RATIO[1], ctx.mpf(PI_RATIO[0]) / PI_RATIO[1]
    coefficient = convert(ctx, wt) * k / (4 * pi * pi * eps)
    def forcing(t):
        L, E = -1 / t, ctx.exp(ctx.j * 2 * k * (t - ts[0]))
        return coefficient * (6 * L ** 3 * c.real + ctx.re(((4 * k * k * L + 6 * L ** 3) * d + 2 * L * L * w0) * E))
    quadrature = ctx.quad(forcing, [ts[0], ts[-1]])
    require(abs(quadrature - profiles[0]['I_ab']) < ctx.mpf('1e-85'),
            'Primitive disagrees with independent high-precision quadrature')
    record('independent_F_quadrature_matches_primitive')
    # Fabricated cancellation histories verify exact inherited3-column ordering and global reset.
    history = np.array([[LD(2) ** 50, LD('-0.125'), LD('0.25')],
                        [LD(1), LD('0.5'), LD('-0.75')],
                        [-LD(2) ** 50, LD('0.125'), LD('0.25')],
                        [LD('0.125'), LD('0.375'), LD('0.5')],
                        [LD('0.625'), LD('-0.25'), LD('0.125')],
                        [LD('-0.125'), LD('0.25'), LD('0.375')],
                        [LD('0.875'), LD('0.5'), LD('-0.125')],
                        [LD('0.75'), LD('0.125'), LD('0.25')],
                        [LD('0.375'), LD('-0.5'), LD('0.625')]], dtype=LD)
    dt = LD(1) / LD(128)
    def literal(end, data=history, delta=dt):
        return delta / 3 * (data[0] + data[end] + 4 * np.sum(data[1:end:2], axis=0, dtype=LD) + 2 * np.sum(data[2:end:2], axis=0, dtype=LD))
    actual = simpson_diagnostics(history, 4, 8, dt, True)
    require(np.array_equal(actual['S_ab'], literal(8) - literal(4)), 'Global prefix arithmetic changed')
    decimated = history[::2]
    require(np.array_equal(actual['S_ab_double_step'], literal(4, decimated, dt * LD(2)) - literal(2, decimated, dt * LD(2))), 'Double step global arithmetic changed')
    record('literal_global_LD_prefix_native_double_and_subtraction')
    for value in (ctx.mpf('0'), ctx.mpf('1e-77'), ctx.mpf('-987654.123456789')):
        text = decimal(ctx, value)
        require(isinstance(text, str) and abs(ctx.mpf(text) - value) <= ctx.mpf('1e-12'), 'Decimal serialization failed')
    record('decimal_string_serialization')
    return {'status': 'PASS_SYNTHETIC_ONLY', 'physical_values_loaded': False,
            'checks_passed': len(checks), 'checks': checks,
            'optimized': not __debug__, 'numpy_version': np.__version__,
            'longdouble_nmant': int(np.finfo(LD).nmant)}


if __name__ == '__main__':
    result = run()
    destination = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if destination:
        require(not destination.exists(), 'Synthetic receipt must be new')
        destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
