"""Deferred interval budgets for derivatives of the two inherited compact pulses.

Import performs no source, root, mode, response, or interval evaluation. Call
only after the enclosing stress registration is frozen. Unit amplitude and
unit half-width are used: u=eta+4, f=B or uB, B=exp(1-1/(1-u^2)).
The public result encloses N_j=|f^(j+2)(u)|+integral_-1^u |f^(j+3)(v)|dv
for j=0,1,2. The function accepts rational observation times; if a critical
point cannot be ordered at the fixed root precision it fails explicitly.
"""
from fractions import Fraction
from functools import lru_cache
from math import inf, nextafter


ROOT_BISECTIONS = 256
INTERVAL_DECIMAL_PRECISION = 60
CENTER = Fraction(-4)
HALF_WIDTH = Fraction(1)
AMPLITUDE = Fraction(1)
INHERITED_OBSERVATIONS = tuple(Fraction(x) for x in
                               ("-5.5", "-4.5", "-4", "-3.5", "-2.5", "-1.5"))

# Descending coefficients of P_n(u), f^(n)=P_n(u) B/(1-u^2)^(2n).
DERIVATIVE_POLYNOMIALS = {
    "positive_B": {
        0: (1,),
        1: (-2, 0),
        2: (6, 0, 0, 0, -2),
        3: (-24, 0, -12, 0, 40, 0, -12, 0),
        4: (120, 0, 180, 0, -528, 0, 232, 0, 24, 0, -12),
        5: (-720, 0, -2160, 0, 6120, 0, -2400, 0, -2112, 0, 1360, 0, -120, 0),
    },
    "signed_uB": {
        0: (1, 0),
        1: (1, 0, -4, 0, 1),
        2: (2, 0, 8, 0, -6, 0),
        3: (-6, 0, -48, 0, 52, 0, 0, 0, -6),
        4: (24, 0, 324, 0, -368, 0, -184, 0, 280, 0, -60, 0),
        5: (-120, 0, -2460, 0, 2280, 0, 4940, 0, -6952, 0, 2220, 0, 120, 0, -60),
    },
}

# Roots in z=u^2 of f^(m+1), for total variation of f^m, m=2,3,4.
# Constant nonzero factors are removed. All brackets have width 1/256.
CRITICAL_POLYNOMIALS = {
    "positive_B": {
        2: (6, 3, -10, 3),
        3: (30, 45, -132, 58, 6, -3),
        4: (90, 270, -765, 300, 264, -170, 15),
    },
    "signed_uB": {
        2: (3, 24, -26, 0, 3),
        3: (6, 81, -92, -46, 70, -15),
        4: (30, 615, -570, -1235, 1738, -555, -30, 15),
    },
}
CRITICAL_BRACKET_NUMERATORS = {
    "positive_B": {2: (95, 205), 3: (59, 184, 224), 4: (27, 166, 213, 233)},
    "signed_uB": {2: (116, 207), 3: (76, 188, 224), 4: (47, 170, 214, 233)},
}
ZERO_IS_CRITICAL = {"positive_B": {2: True, 3: False, 4: True},
                    "signed_uB": {2: False, 3: True, 4: False}}


def _polynomial(coefficients, value):
    result = 0
    for coefficient in coefficients:
        result = result*value+coefficient
    return result


@lru_cache(maxsize=None)
def _root_bracket(source, order, index):
    numerator = CRITICAL_BRACKET_NUMERATORS[source][order][index]
    left, right = Fraction(numerator,256), Fraction(numerator+1,256)
    coefficients = CRITICAL_POLYNOMIALS[source][order]
    fl = _polynomial(coefficients,left)
    if fl*_polynomial(coefficients,right) >= 0:
        raise RuntimeError("Frozen critical-root bracket lost its sign change")
    for _ in range(ROOT_BISECTIONS):
        middle = (left+right)/2
        fm = _polynomial(coefficients,middle)
        if fm == 0:
            return middle,middle
        if fl*fm < 0:
            right = middle
        else:
            left,fl = middle,fm
    return left,right


def _rational(iv, value):
    return iv.mpf(value.numerator)/iv.mpf(value.denominator)


def _root_before_endpoint(source, order, index, sign, endpoint):
    """Exact ordering only; no floating source or root approximation."""
    left,right = _root_bracket(source,order,index)
    if sign > 0:
        if endpoint <= 0:
            return False
        if right <= endpoint**2:
            return True
        if left >= endpoint**2:
            return False
    else:
        if endpoint >= 0:
            return True
        if left >= endpoint**2:
            return True
        if right <= endpoint**2:
            return False
    if _polynomial(CRITICAL_POLYNOMIALS[source][order],endpoint**2) == 0:
        return True
    raise RuntimeError("Endpoint and critical root cannot be ordered at frozen precision")


def _derivative(iv, source, order, u):
    denominator = 1-u*u
    return (_polynomial(DERIVATIVE_POLYNOMIALS[source][order],u)
            *iv.exp(1-1/denominator)/denominator**(2*order))


def _upper(value):
    if value == 0:
        return 0.0
    return nextafter(float(value.b),inf)


def pulse_interval_data(iv, source, eta):
    """Return normalized jets f0..f5 and N0..N2 in the supplied interval context.

    Runtime interval/source evaluations belong to the separately frozen run.
    The caller must supply an MPIntervalContext with dps exactly 60.
    """
    if source not in DERIVATIVE_POLYNOMIALS:
        raise ValueError("Source must be positive_B or signed_uB")
    time = Fraction(str(eta))
    endpoint = (time-CENTER)/HALF_WIDTH
    if iv.dps != INTERVAL_DECIMAL_PRECISION:
        raise ValueError("Pulse budgets require the frozen 60-digit interval context")
    zero = iv.mpf(0)
    details = {"source":source,"eta":str(eta),"amplitude":"1","half_width":"1",
               "root_bisections":ROOT_BISECTIONS,"interval_dps":INTERVAL_DECIMAL_PRECISION,
               "method":"Exact critical-root total variation with directed interval primitive evaluation",
               "positive_critical_root_counts":{
                   str(m):len(CRITICAL_BRACKET_NUMERATORS[source][m]) for m in (2,3,4)}}
    if endpoint <= -1:
        return {"jets":[zero for _ in range(6)],"budgets":[zero for _ in range(3)],"details":details}
    local = ({n:zero for n in range(6)} if endpoint >= 1 else
             {n:_derivative(iv,source,n,_rational(iv,endpoint)) for n in range(6)})
    stop = min(endpoint,Fraction(1))
    norms = {}
    for j in range(3):
        order = j+2
        count = len(CRITICAL_BRACKET_NUMERATORS[source][order])
        roots = [(-1,i) for i in reversed(range(count))]
        if ZERO_IS_CRITICAL[source][order]:
            roots.append((0,None))
        roots.extend((1,i) for i in range(count))
        previous = zero
        variation = zero
        for sign,index in roots:
            if sign == 0:
                if stop < 0:
                    continue
                value = _derivative(iv,source,order,zero)
            else:
                if not _root_before_endpoint(source,order,index,sign,stop):
                    continue
                left,right = _root_bracket(source,order,index)
                z = iv.mpf([_rational(iv,left).a,_rational(iv,right).b])
                value = _derivative(iv,source,order,sign*iv.sqrt(z))
            variation += abs(value-previous)
            previous = value
        variation += abs(local[order]-previous)
        norms[j] = abs(local[order])+variation
    return {"jets":[local[n] for n in range(6)],"budgets":[norms[j] for j in range(3)],"details":details}


def normalized_budget_intervals(source, eta):
    """Convenience wrapper constructing the same frozen interval context."""
    from mpmath.ctx_iv import MPIntervalContext
    iv = MPIntervalContext()
    iv.dps = INTERVAL_DECIMAL_PRECISION
    result = pulse_interval_data(iv,source,eta)
    return iv,dict(enumerate(result["budgets"])),dict(enumerate(result["jets"]))


def normalized_derivative_budgets(source, eta):
    """Outward binary64 upper bounds for all three N_j and local derivatives."""
    _,norms,local = normalized_budget_intervals(source,eta)
    return {"source":source,"eta":str(eta),"amplitude":"1","half_width":"1",
            "normalized_N_upper":{str(j):_upper(norms[j]) for j in range(3)},
            "local_abs_derivative_upper":{str(n):_upper(abs(local[n])) for n in range(6)},
            "root_bisections":ROOT_BISECTIONS,"interval_dps":INTERVAL_DECIMAL_PRECISION,
            "scope":"Analytic derivative-norm enclosures only; no response/mode/finite-quadrature error estimate"}
