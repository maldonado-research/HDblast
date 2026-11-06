"""Exact arithmetic for known quadratic proper-cut sources at e=0.

This consumes prior analytic/PDE bounds. The finite checks do not prove
those premises, the continuous identities, full curvature or a radius.
"""
from fractions import Fraction as F
import hashlib
import itertools
import json
from math import factorial

SPHERE_DIMENSION = 4


def need(ok, why):
    if not ok:
        raise ValueError(why)


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def eye(d):
    return [[F(i == j) for j in range(d)] for i in range(d)]


def moment(powers):
    if any(p % 2 for p in powers):
        return F(0)
    numerator = 1
    for p in powers:
        for k in range(1, p, 2):
            numerator *= k
    denominator = 1
    for k in range(sum(powers)//2):
        denominator *= SPHERE_DIMENSION+2*k
    return F(numerator, denominator)


def integrate_y_power(power, extra=(0, 0, 0, 0)):
    total = F(0)
    for k in range(power+1):
        p = list(extra)
        p[0] += 2*k
        p[1] += 2*(power-k)
        choose = factorial(power)//(factorial(k)*factorial(power-k))
        total += (-1)**(power-k)*choose*moment(p)
    return total


def det_coefficient(matrix, order):
    """Leibniz determinant of degree-two entry polynomials, no jet formula."""
    value = F(0)
    for perm in itertools.permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        for deg in itertools.product(range(3), repeat=3):
            if sum(deg) == order:
                product = F((-1)**inversions)
                for i in range(3):
                    product *= matrix[i][perm[i]][deg[i]]
                value += product
    return value


def evaluate(root):
    pin_doc = json.loads((root/'SOURCE_PINS.json').read_bytes())
    need(pin_doc['schema'] == 'HDBLAST_QUADRATIC_PROPER_AREA_SOURCE_PINS_V1', 'Source pin schema differs')
    pins = pin_doc['sources']
    need(len(pins) == 4 and len({r['local'] for r in pins}) == 4, 'Source roster differs')
    for row in pins:
        raw = (root/row['local']).read_bytes()
        need(len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'], 'Retained source identity differs')
        need(hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == row['git_blob_sha1'], 'Retained blob identity differs')
        need(row['ref'] == 'a715e8ee2963467b8988b7984a1cd0d402468f3d', 'Source commit differs')
    higher = json.loads((root/'inputs/HIGHER_JETS_NORMAL.json').read_bytes())
    need(higher['status'] == 'PASS' and higher['schema'] == 'hdblast-interior-master-jets-v1', 'Retained report status/schema differs')
    need(higher['input_pins']['inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json'] ==
         'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8', 'Same-state input join differs')
    need(higher['source_bound_adapter']['ell'] == '1/1000000' and
         higher['source_bound_adapter']['domain'] == '0<=x<=3, r>=0, x+r<=8', 'Scaled trace domain differs')
    passed = []

    def check(name, ok):
        need(ok, name)
        passed.append(name)

    check('immutable_actual_input_chain', True)
    ell, rho, U = F(1, 10**6), F(129), F(23, 2*10**9)
    T, qD = 25*ell/16, ell/4
    D, fullD = T-ell, T-qD
    y, yp = (F(higher['exact']['boundary'][key]) for key in ['y', 'yprime'])
    check('actual_first_field_and_distinct_quiet_durations',
          0 < y < 189*ell and 190*ell/rho**2 < U and
          D == 9*ell/16 and fullD == 21*ell/16 and T < 3*ell)
    omega_upper = F(315, 4)*ell**4
    transfer = F(253, 4000)
    check('complete_energy_derivative_endpoint_transfer',
          omega_upper < 5500*ell**2 and F(4, 1000) == F(1, 250) < transfer**2)
    Imax, Jmax = U*D, U*D**2/2
    Bmax = U*(1+T*D+F(4, 3)*D**2)
    uq = (yp/ell+transfer+F(3, 2)*T*190*ell)/rho**2
    Fmax = 3*uq+3*T*U+5*Imax
    check('actual_past_first_time_derivative_bound', uq < F(293, 10000) and Fmax < F(11, 125))
    mY2, mY4 = integrate_y_power(2), integrate_y_power(4)
    mY2s = sum(integrate_y_power(2, extra) for extra in [(2, 0, 0, 0), (0, 2, 0, 0)])
    check('sphere Y2 moment mismatch', mY2 == F(1, 6))
    mK = 8*(1-moment((2, 0, 0, 0))-moment((0, 2, 0, 0)))+F(8, 3)*mY2
    mY2K = 8*(mY2-mY2s)+F(8, 3)*mY4
    check('exact_round_and_denominator_weighted_moments',
          mY4 == F(3, 40) and mY2s == F(1, 8) and mK == F(40, 9) and
          mY2K == F(8, 15) and mY4/mY2 == F(9, 20) and mY2K/mY2 == F(16, 5) and
          integrate_y_power(3) == integrate_y_power(5) == 0)
    fixtures = [(F(1), F(0), F(0), F(0)), (F(0), F(0), F(1), F(0)), (F(0), F(0), F(0), F(1))]
    for a, b, c in [(F(1, 2), F(1, 3), F(1, 5)), (F(2), F(-1), F(1, 2)), (F(0), F(1), F(2)), (F(1), F(1), F(1))]:
        den = 1+a*a+b*b+c*c
        fixtures.append((2*a/den, 2*b/den, 2*c/den, (1-a*a-b*b-c*c)/den))
    amat = [[F((1, -1, 0, 0)[i] if i == j else 0) for j in range(4)] for i in range(4)]
    id4 = eye(4)
    compressed = []
    angular_ok = True
    for n in fixtures:
        Y, s = n[0]**2-n[1]**2, n[0]**2+n[1]**2
        p = [[id4[i][j]-n[i]*n[j] for j in range(4)] for i in range(4)]
        at = mm(mm(p, amat), p)
        h = [[2*at[i][j]+F(2, 3)*Y*p[i][j] for j in range(4)] for i in range(4)]
        K = 8*(1-s)+F(8, 3)*Y**2
        angular_ok &= (sum(t*t for t in n) == 1 and tr(h) == 0 and tr(mm(h, h)) == K and 0 <= K <= 8)
        v = [F(i == 3)-n[i] for i in range(4)]
        norm2 = sum(t*t for t in v)
        house = id4 if norm2 == 0 else [[id4[i][j]-2*v[i]*v[j]/norm2 for j in range(4)] for i in range(4)]
        tangent = [[house[i][j] for j in range(3)] for i in range(4)]
        ht = mm(mm([list(x) for x in zip(*tangent)], h), tangent)
        angular_ok &= tr(mm(ht, ht)) == K
        compressed.append((Y, ht, K))
    check('rational_ambient_projection_and_tangent_compression_fixtures', angular_ok)
    det_ok = True
    g2 = [[F(2), F(3), F(-1)], [F(3), F(-4), F(2)], [F(-1), F(2), F(5)]]
    for Y, h, K in compressed:
        for B, J in [(F(1), F(0)), (F(0), F(1)), (F(2, 3), F(-5, 7))]:
            g1 = [[2*B*Y*F(i == j)-2*J*h[i][j] for j in range(3)] for i in range(3)]
            polynomial = [[[F(i == j), g1[i][j], g2[i][j]/2] for j in range(3)] for i in range(3)]
            d1, d2 = det_coefficient(polynomial, 1), det_coefficient(polynomial, 2)
            observed_j2 = d2-d1*d1/4
            expected_j2 = tr(g2)/2+3*B*B*Y*Y-2*J*J*K
            det_ok &= observed_j2 == expected_j2
    check('leibniz_determinant_second_derivative_with_noncommuting_G2', det_ok)
    det_mean_pos, det_mean_neg = Bmax**2/2, F(80, 9)*Jmax**2
    det_w_pos, det_w_neg = F(27, 20)*Bmax**2, F(32, 5)*Jmax**2
    check('known_determinant_and_weighted_denominator_enclosures',
          3*mY4/mY2 == F(27, 20) and 2*mY2K/mY2 == F(32, 5) and
          det_mean_pos < F(6613, 10**20) and det_mean_neg < F(295, 10**43) and
          det_w_pos < F(1786, 10**19) and det_w_neg < F(212, 10**43))
    Cosh = (1-T*T/2)**-2
    check('rational_cosh_majorant', T*T/2 < 1 and Cosh < 1+F(3, 10**12))
    # Formal second derivative of -theta^2/3-Sigma^2 at theta0=3tanh(q).
    theta0, theta1, theta2 = F(2, 5), F(7, 3), F(-11, 4)
    second_square = 2*theta1**2+2*theta0*theta2
    check('raychaudhuri_derivative_factors_and_weighted_sources',
          -second_square/3+2*theta0*theta2/3 == -F(2, 3)*theta1**2 and
          F(2, 3)*mY4/mY2 == F(3, 10) and 2*mY2K/mY2 == F(32, 5) and
          9*mY4/mY2 == F(81, 20))
    mean_exp = Cosh*Fmax**2*D**2/18
    mean_shear = F(20, 27)*U**2*D**4
    weighted_exp = F(3, 20)*Cosh*Fmax**2*D**2
    weighted_shear = F(8, 15)*U**2*D**4
    variance = F(81, 20)*Bmax**2
    known_neg = weighted_exp+weighted_shear
    curvature = Cosh*fullD**2/2
    check('known_quadratic_ray_and_variance_budgets',
          mean_exp < F(1352, 10**19) and mean_shear < F(981, 10**44) and
          weighted_exp+weighted_shear < F(365, 10**18) and weighted_shear < F(7062, 10**45) and
          variance < F(5357, 10**19) and 2*Imax/rho < F(101, 10**18))
    check('curvature_projection_transfer_uses_full_quiet_duration',
          curvature < F(8614, 10**16) and curvature > Cosh*D**2/2)
    exact = {
        'ell':ell, 'rho_min':rho, 'target_q':T, 'nonlinear_quiet_q':qD,
        'first_active_duration':D, 'full_active_duration':fullD,
        'U':U, 'future_y_bound':y, 'future_yprime_bound':yp,
        'endpoint_derivative_transfer':transfer, 'I_upper':Imax, 'J_upper':Jmax,
        'B_upper':Bmax, 'u_q_upper':uq, 'theta1_coefficient_upper':Fmax,
        'cosh_squared_upper':Cosh, 'mean_determinant_positive_upper':det_mean_pos,
        'mean_determinant_negative_magnitude':det_mean_neg,
        'weighted_determinant_positive_upper':det_w_pos,
        'weighted_determinant_negative_magnitude':det_w_neg,
        'mean_log_expansion_magnitude':mean_exp, 'mean_log_shear_magnitude':mean_shear,
        'weighted_log_expansion_magnitude':weighted_exp,
        'weighted_log_shear_magnitude':weighted_shear,
        'weighted_j1_square_upper':variance,
        'known_denominator_negative_magnitude':known_neg,
        'curvature_projection_multiplier_upper':curvature,
    }
    return {
        'schema':'HDBLAST_QUADRATIC_PROPER_AREA_MATHEMATICS_V1',
        'status':'PASS', 'named_checks':passed, 'named_check_count':len(passed),
        'angular_fixtures':len(fixtures), 'determinant_fixtures':len(fixtures)*3,
        'exact':{key:str(value) for key,value in exact.items()},
        'scope':'Finite exact arithmetic, retained input joins and polynomial fixtures. '
                'Continuous observer/Raychaudhuri and analytic/PDE bounds require separate review. '
                'Known e=0 quadratic sources only; actual R2 projection, full G2, Q(a), '
                'uniform nonlinear remainder, radius and cosmological origin are not established.',
    }
