"""Unrestricted bare bilinears, independent Fourier moments, and saved LD Simpson.

This module performs no file reads and no source/model evolution. It is called
only after the separate CLI's registration and immutable input checks.
"""
from __future__ import annotations
import numpy as np
import mpmath
from decimal import Decimal
from fractions import Fraction
from exact_binary import exact_binary, exact_complex, require
from profile_kernel import accumulate_node, empty_sums

LD = np.longdouble
CD = np.clongdouble
EPS_RATIO = (3777893186295716171, 37778931862957161709568)
PI_RATIO = (14488038916154245685, 4611686018427387904)
CUTOFFS = (64, 128, 256)
NAMES = ('q', 'q_prime', 'q_second', 'rho', 'p', 'Q0', 'rho0', 'p0', 'current')
GEOMETRY_NAMES = ('a0', 'L', 'L_prime', 'L_second', 'L_third', 'M',
                  'M_prime', 'M_second', 'M_third', 'M_fourth')
SETTINGS = {
    'coarse': {'dt_denominator': 128, 'momentum_width_denominator': 2,
               'nodes': 8192, 'history_samples': 577, 'profile_samples': 129},
    'fine': {'dt_denominator': 256, 'momentum_width_denominator': 4,
             'nodes': 16384, 'history_samples': 1153, 'profile_samples': 257},
}
MEMBERS = ('k', 'momentum_weights', 'u_4', 'w_4', 'u_5', 'w_5',
           'history_eta', 'history_values', 'history_ledger_integrand',
           'history_baseline_contact', 'history_source_jet', 'history_forcing_jet',
           'history_geometry', 'history_geometry_names',
           'history_quantity_names', 'history_cutoffs', 'observation_eta')


def context(dps):
    require(dps in (80, 100), 'Only the registered 80/100 decimal precision is allowed')
    ctx = mpmath.mp.clone()
    ctx.dps = dps
    return ctx


def convert(ctx, value):
    """No binary64 intermediate and no decimal rendering of retained inputs."""
    return ctx.make_mpf(exact_binary(value).mpf_tuple())


def convert_complex(ctx, value):
    real, imag = exact_complex(value)
    return ctx.mpc(ctx.make_mpf(real.mpf_tuple()), ctx.make_mpf(imag.mpf_tuple()))


def decimal(ctx, value):
    """Serialize beyond the working precision and check the decimal roundtrip."""
    require(ctx.isfinite(value), 'Nonfinite diagnostic output')
    rendered = ctx.nstr(value, ctx.dps + 8, strip_zeros=False, min_fixed=-8, max_fixed=8)
    sign, mantissa, exponent, _ = value._mpf_
    numerator = -mantissa if sign else mantissa
    exact_value = Fraction(numerator << exponent, 1) if exponent >= 0 \
        else Fraction(numerator, 1 << -exponent)
    require(abs(Fraction(Decimal(rendered)) - exact_value) <= Fraction('1e-12'),
            'Decimal serialization gap exceeds the fixed diagnostic gate')
    return rendered


def simpson_prefix(history, end, dt):
    """Literal inherited arithmetic: all three cutoff columns summed together."""
    require(end >= 2 and end % 2 == 0, 'Simpson prefix needs an even positive end')
    return dt / LD(3) * (history[0] + history[end]
        + LD(4) * np.sum(history[1:end:2], axis=0, dtype=LD)
        + LD(2) * np.sum(history[2:end:2], axis=0, dtype=LD))


def simpson_diagnostics(history, anchor, end, dt, fine):
    """Global prefixes subtract in LD before conversion, plus separate reset."""
    native = simpson_prefix(history, end, dt) - simpson_prefix(history, anchor, dt)
    reset = simpson_prefix(history[anchor:end + 1], end - anchor, dt)
    result = {'S_ab': native, 'S_reset_direct_ld': reset}
    if fine:
        decimated = history[::2]
        double = simpson_prefix(decimated, end // 2, LD(2) * dt) \
            - simpson_prefix(decimated, anchor // 2, LD(2) * dt)
        result['S_ab_double_step'] = double
        result['S_reset_direct_ld_double_step'] = simpson_prefix(
            decimated[anchor // 2:end // 2 + 1], (end - anchor) // 2, LD(2) * dt)
    return result


def validate_arrays(arrays, setting):
    require(set(arrays) == set(MEMBERS), 'Input members must match the frozen schema exactly')
    require(setting in SETTINGS, 'Unknown setting')
    spec = SETTINGS[setting]
    n, nt = spec['nodes'], spec['history_samples']
    shapes = {'k': (n,), 'momentum_weights': (n,), 'u_4': (n,), 'w_4': (n,),
              'u_5': (n,), 'w_5': (n,), 'history_eta': (nt,),
              'history_values': (nt, 3, 9), 'history_ledger_integrand': (nt, 3),
              'history_baseline_contact': (nt, 3), 'history_source_jet': (nt, 6),
              'history_forcing_jet': (nt, 4), 'history_geometry': (nt, 10),
              'history_cutoffs': (3,), 'observation_eta': (6,)}
    for name, shape in shapes.items():
        value = arrays[name]
        require(isinstance(value, np.ndarray) and value.shape == shape,
                'Frozen input shape mismatch: ' + name)
        expected_dtype = np.dtype(np.int64 if name == 'history_cutoffs' else
                                  CD if name in ('u_4', 'w_4', 'u_5', 'w_5') else LD)
        require(value.dtype == expected_dtype, 'Frozen input dtype mismatch: ' + name)
        require(np.all(np.isfinite(value)), 'Nonfinite retained array: ' + name)
    require(tuple(arrays['history_quantity_names'].tolist()) == NAMES,
            'History quantity order differs from the inherited producer')
    require(arrays['history_quantity_names'].shape == (9,)
            and arrays['history_quantity_names'].dtype == np.dtype('<U8'),
            'History quantity name dtype/shape mismatch')
    require(arrays['history_geometry_names'].shape == (10,)
            and arrays['history_geometry_names'].dtype == np.dtype('<U8')
            and tuple(arrays['history_geometry_names'].tolist()) == GEOMETRY_NAMES,
            'History geometry name dtype/shape/order mismatch')
    require(np.array_equal(arrays['history_cutoffs'], np.array(CUTOFFS, dtype=np.int64)),
            'Cutoff order differs from the registration')
    expected_observations = np.array(['-5.5', '-4.5', '-4', '-3.5', '-2.5', '-1.5'], dtype=LD)
    require(np.array_equal(arrays['observation_eta'], expected_observations),
            'Observation anchors differ from the registration')
    times = arrays['history_eta']
    dt = LD(1) / LD(spec['dt_denominator'])
    require(np.array_equal(times, LD(-6) + np.arange(nt, dtype=LD) * dt),
            'Retained history is not the exact registered dyadic grid')
    require(np.all(arrays['k'] > 0) and np.all(np.diff(arrays['k']) > 0),
            'Momentum nodes must be positive and strictly increasing')
    require(np.all(arrays['momentum_weights'] > 0), 'Momentum weights must be positive')
    boundaries = tuple(K * spec['momentum_width_denominator'] * 16 for K in CUTOFFS)
    for K, count in zip(CUTOFFS, boundaries):
        require(arrays['k'][count - 1] < LD(K), 'Cutoff node prefix ends outside its panel')
        if count < n:
            require(arrays['k'][count] > LD(K), 'Cutoff node prefix starts inside prior panel')
    mask = (times >= LD('-2.5')) & (times <= LD('-1.5'))
    indices = np.flatnonzero(mask)
    require(len(indices) == spec['profile_samples'], 'Wrong inclusive profile mask')
    require(indices[0] % 2 == 0 and indices[-1] % 2 == 0,
            'Selected endpoints are not Simpson prefix endpoints')
    for name in ('history_source_jet', 'history_forcing_jet', 'history_baseline_contact'):
        require(np.count_nonzero(arrays[name][indices]) == 0,
                'Saved source, forcing or baseline contact is nonzero on selected interval: ' + name)
    return indices, boundaries, dt


def calculate_profiles(k_values, weights, u_anchor, w_anchor, u_end, w_end,
                       sample_times, boundaries, dps, flow_sink=None, progress=None):
    """Independent free-flow calculation. The synthetic tests also call this core.

    Weighted unrestricted stresses are grouped into Re(Z), Omega Im(Z), k²Re(Z).
    The separately expanded F primitive is G=3L²(C0+ReZ)+L OmegaImZ.
    G is not evaluated from a measured density or a Ward-defined pressure.
    """
    ctx = context(dps)
    eps = ctx.mpf(EPS_RATIO[0]) / EPS_RATIO[1]
    pi = ctx.mpf(PI_RATIO[0]) / PI_RATIO[1]
    normalization = 1 / (4 * pi * pi * eps)
    times = [convert(ctx, value) for value in sample_times]
    require(len(times) >= 2, 'At least two profile samples required')
    first, last = times[0], times[-1]
    dt = times[1] - times[0]
    require(all(times[j] == first + j * dt for j in range(len(times))),
            'Profile phase recurrence requires an exactly uniform grid')
    L_values = [-1 / value for value in times]
    La, Lb = L_values[0], L_values[-1]
    sums = empty_sums(len(times))
    C0 = ctx.mpf(0)
    Ck2 = ctx.mpf(0)
    signed_flow = ctx.mpf(0)
    direct_ledger = ctx.mpf(0)
    triangle_flow = ctx.mpf(0)
    max_du = ctx.mpf(0)
    max_dw = ctx.mpf(0)
    result = []
    cut_index = 0
    for index in range(len(k_values)):
        k = convert(ctx, k_values[index])
        weight = convert(ctx, weights[index])
        ua, wa = convert_complex(ctx, u_anchor[index]), convert_complex(ctx, w_anchor[index])
        ub, wb = convert_complex(ctx, u_end[index]), convert_complex(ctx, w_end[index])
        omega = 2 * k
        k2 = k * k
        coefficient = normalization * weight * k
        d = ctx.mpc(wa.imag / omega, -wa.real / omega)
        creal = ua.real - d.real
        C0 += coefficient * creal
        Ck2 += coefficient * creal * k2
        z = coefficient * d
        accumulate_node((z.real._mpf_, z.imag._mpf_), omega._mpf_, k2._mpf_,
                        ctx.mpf(0)._mpf_, dt._mpf_, len(times), sums, ctx.prec)
        phase = ctx.exp(ctx.j * omega * (last - first))
        # Canonical ledger uses fresh directMPexp, independently of the profile recurrence.
        weighted_ledger = coefficient * (3 * creal * (Lb * Lb - La * La)
            + ctx.re(d * (3 * (Lb * Lb * phase - La * La)
                           - ctx.j * omega * (Lb * phase - La))))
        direct_ledger += weighted_ledger
        du = ub - ua - (phase - 1) * d
        dw = wb - phase * wa
        flow = coefficient * ((2 * k2 + 3 * Lb * Lb) * du.real
                               - k * dw.imag - Lb * dw.real)
        triangle = coefficient * ((2 * k2 + 3 * Lb * Lb) * abs(du)
                                   + (k + Lb) * abs(dw))
        signed_flow += flow
        triangle_flow += triangle
        max_du = max(max_du, abs(du))
        max_dw = max(max_dw, abs(dw))
        if flow_sink:
            flow_sink({'node_index': index,
                       'k': decimal(ctx, k), 'delta_u_real': decimal(ctx, du.real),
                       'delta_u_imag': decimal(ctx, du.imag),
                       'delta_w_real': decimal(ctx, dw.real),
                       'delta_w_imag': decimal(ctx, dw.imag),
                       'weighted_signed_flow': decimal(ctx, flow),
                       'weighted_triangle': decimal(ctx, triangle),
                       'weighted_direct_ledger': decimal(ctx, weighted_ledger)})
        if index + 1 == boundaries[cut_index]:
            R, P, F, G = [], [], [], []
            for j, L in enumerate(L_values):
                zr = ctx.make_mpf(sums[0][j])
                omegazi = ctx.make_mpf(sums[1][j])
                k2zr = ctx.make_mpf(sums[2][j])
                R.append(2 * Ck2 + 3 * L * L * (C0 + zr) + L * omegazi)
                P.append(ctx.mpf(2) / 3 * Ck2 - L * L * (C0 + zr)
                         - ctx.mpf(4) / 3 * k2zr + L * omegazi)
                F.append(6 * L ** 3 * (C0 + zr) + 4 * L * k2zr - 2 * L * L * omegazi)
                G.append(3 * L * L * (C0 + zr) + L * omegazi)
            result.append({'R': R, 'P': P, 'F': F, 'I_ab': +direct_ledger,
                           'I_recurrence': G[-1] - G[0],
                           'I_recurrence_minus_direct': G[-1] - G[0] - direct_ledger,
                           'E_flow': +signed_flow, 'E_flow_triangle': +triangle_flow,
                           'max_abs_delta_u': +max_du, 'max_abs_delta_w': +max_dw})
            cut_index += 1
            if cut_index == len(boundaries):
                break
        if progress and (index + 1) % 256 == 0:
            progress()
    require(cut_index == len(boundaries), 'Not all cutoff prefixes were processed')
    return ctx, times, result


def evaluate_arrays(arrays, setting, dps, flow_sink=None, progress=None):
    indices, boundaries, dt = validate_arrays(arrays, setting)
    sample_times = arrays['history_eta'][indices]
    ctx, times, profiles = calculate_profiles(
        arrays['k'], arrays['momentum_weights'], arrays['u_4'], arrays['w_4'],
        arrays['u_5'], arrays['w_5'], sample_times, boundaries, dps, flow_sink, progress)
    simpsons = simpson_diagnostics(arrays['history_ledger_integrand'],
                                  int(indices[0]), int(indices[-1]), dt, setting == 'fine')
    rows = []
    for cutoff_index, K in enumerate(CUTOFFS):
        profile = profiles[cutoff_index]
        measured_R = [convert(ctx, value) for value in arrays['history_values'][indices, cutoff_index, 3]]
        measured_P = [convert(ctx, value) for value in arrays['history_values'][indices, cutoff_index, 4]]
        measured_F = [convert(ctx, value) for value in arrays['history_ledger_integrand'][indices, cutoff_index]]
        DeltaR = measured_R[-1] - measured_R[0]
        S_ab = convert(ctx, simpsons['S_ab'][cutoff_index])
        I_ab = profile['I_ab']
        D_S, D_cont, E_Q = DeltaR - S_ab, DeltaR - I_ab, I_ab - S_ab
        R_errors = [a - b for a, b in zip(measured_R, profile['R'])]
        P_errors = [a - b for a, b in zip(measured_P, profile['P'])]
        F_errors = [a - b for a, b in zip(measured_F, profile['F'])]
        scalars = {name: profile[name] for name in ('I_ab', 'I_recurrence',
                                                   'I_recurrence_minus_direct', 'E_flow',
                                                   'max_abs_delta_u', 'max_abs_delta_w')}
        scalars['triangle_bound'] = profile['E_flow_triangle']
        scalars.update({'S_ab': S_ab, 'DeltaR': DeltaR, 'D_S': D_S, 'D_cont': D_cont,
                        'E_Q': E_Q, 'decomposition_error': D_S - D_cont - E_Q,
                        'flow_operator_closure': D_cont - profile['E_flow'],
                        'R_profile_max_error': max(map(abs, R_errors)),
                        'P_profile_max_error': max(map(abs, P_errors)),
                        'F_profile_max_error': max(map(abs, F_errors))})
        scalars['S_direct_reset'] = convert(ctx, simpsons['S_reset_direct_ld'][cutoff_index])
        scalars['S_global_minus_reset'] = S_ab - scalars['S_direct_reset']
        if setting == 'fine':
            scalars['S_ab_double_global'] = convert(ctx, simpsons['S_ab_double_step'][cutoff_index])
            scalars['S_direct_reset_double'] = convert(ctx, simpsons['S_reset_direct_ld_double_step'][cutoff_index])
            scalars['D_S_double_step'] = DeltaR - scalars['S_ab_double_global']
            scalars['E_Q_double_step'] = I_ab - scalars['S_ab_double_global']
            scalars['S_native_minus_double'] = S_ab - scalars['S_ab_double_global']
        full_profile = {'times': times, 'R': profile['R'], 'P': profile['P'], 'F': profile['F'],
                        'stored_R': measured_R, 'stored_P': measured_P,
                        'stored_F': measured_F,
                        'signed_R_error': R_errors, 'signed_P_error': P_errors}
        rows.append({'K': K, **{name: decimal(ctx, value) for name, value in scalars.items()},
                     **{name: [decimal(ctx, value) for value in values]
                        for name, values in full_profile.items()}})
    return rows
