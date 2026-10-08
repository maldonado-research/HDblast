"""Arb incoming-target adapter; no real callback is evaluated at import.

Adapted from the frozen 3 October verified_moments machinery. This candidate
uses exact geometric parent geometry and exact affine child recentering.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial, comb, ceil
from flint import arb, acb, arb_series, acb_poly, fmpq, ctx

DEGREE = 112
PHASE_DEGREE = 128
PRECISION_BITS = 512
DELTA = F(1, 128)
A, B = F(-5) + DELTA, F(-9, 2)

def ball(x):
    x = F(x)
    return arb(fmpq(x.numerator, x.denominator))

def inflate_real(x, radius):
    radius = F(radius)
    if radius < 0:
        raise ValueError('Negative remainder')
    return x + arb(0, ball(radius).upper())

def inflate_complex(x, radius):
    return acb(inflate_real(x.real, radius), inflate_real(x.imag, radius))

def upper(x):
    return F(str(x.abs_upper().fmpq()))

def qstring(x):
    x = F(str(x))
    return f'{x.numerator}/{x.denominator}'

def rectangle(x):
    if isinstance(x, arb):
        x = acb(x, arb(0))
    if not x.is_finite():
        raise ValueError('Nonfinite enclosure')
    return {name: {'lo': qstring(part.lower().fmpq()), 'hi': qstring(part.upper().fmpq())}
            for name, part in [('real', x.real), ('imag', x.imag)]}

def radius(encoded):
    return sum(((F(encoded[name]['hi']) - F(encoded[name]['lo'])) / 2
                for name in ('real', 'imag')), F(0))

def parent_geometry():
    left = DELTA
    output = []
    while left < F(1, 2):
        right = min(F(3, 2) * left, F(1, 2))
        center = (left + right) / 2
        half_width = (right - left) / 2
        disk_radius = left / 2
        q = 1 - center + disk_radius
        t = 5 - center - disk_radius
        d = 1 - q * q
        mb = 2 * (4/t**2 + 4*q/(t*d**2) + 2/d**2 + 8*q*q/d**3 + 4*q*q/d**4)
        mz = q*mb + 2*(2/t + 4*q/d**2)
        if not (0 < half_width <= disk_radius/2 and q < 1 and t > 0 and d > 0):
            raise ValueError('Unproved parent geometry')
        output.append({'center': -5 + center, 'half_width': half_width,
                       'disk_radius': disk_radius, 'positive_majorant': mb,
                       'signed_majorant': mz})
        left = right
    if len(output) != 11:
        raise ValueError('Incorrect geometric panel count')
    return tuple(output)

PARENTS = parent_geometry()

@dataclass(frozen=True)
class Panel:
    center: F
    half_width: F
    coefficients: tuple
    remainder: F
    parent_index: int

    def __post_init__(self):
        if self.half_width <= 0 or self.remainder < 0 or len(self.coefficients) != DEGREE + 1:
            raise ValueError('Invalid source polynomial')
        if any(not x.is_finite() for x in self.coefficients):
            raise ValueError('Nonfinite source coefficient')

def forcing_coefficients(eta, source, half_width):
    if ctx.cap < DEGREE + 3 or source.prec < DEGREE + 3 or eta.prec < DEGREE + 3:
        raise ValueError('Insufficient source-jet capacity')
    l = -eta.inv()
    hp = source.derivative() / ball(half_width)
    hpp = hp.derivative() / ball(half_width)
    g = 4*l*l*source - 2*l*hp - hpp
    if g.prec < DEGREE + 1:
        raise ValueError('Missing differentiated coefficients')
    return tuple(g[j] for j in range(DEGREE + 1))

def registered_parent(index, source, authorize=None):
    """Root-owned authorization must precede the first real Taylor coefficient."""
    if not callable(authorize):
        raise RuntimeError('No physical source evaluation before public freeze')
    if type(index) is not int or not 0 <= index < len(PARENTS) or source not in ('positive_B', 'signed_uB'):
        raise ValueError('Unregistered source/parent')
    if ctx.prec < PRECISION_BITS:
        raise ValueError('512-bit source arithmetic required')
    p = PARENTS[index]
    authorize('real_source_taylor', {'source': source, 'center': qstring(p['center']),
              'scope': 'UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY'})
    previous_cap = ctx.cap
    try:
        ctx.cap = DEGREE + 3
        eta = arb_series([ball(p['center']), ball(p['half_width'])], prec=DEGREE + 3)
        z = eta + 4
        bump = (1 - (1 - z*z).inv()).exp()
        h = bump if source == 'positive_B' else z*bump
        coefficients = forcing_coefficients(eta, h, p['half_width'])
    finally:
        ctx.cap = previous_cap
    majorant = p['positive_majorant'] if source == 'positive_B' else p['signed_majorant']
    return Panel(p['center'], p['half_width'], coefficients, majorant / 2**DEGREE, index)

def child_panels(parent):
    """Exact subdivision/recentering: source tail is inherited, never sampled."""
    count = ceil(2 * parent.half_width * 64)
    width = 2 * parent.half_width / count
    child_half = width / 2
    previous_cap = ctx.cap
    output = []
    try:
        ctx.cap = DEGREE + 1
        for j in range(count):
            center = parent.center - parent.half_width + (F(j) + F(1, 2)) * width
            alpha = (center - parent.center) / parent.half_width
            beta = child_half / parent.half_width
            affine = arb_series([ball(alpha), ball(beta)], prec=DEGREE + 1)
            poly = arb_series([parent.coefficients[-1]], prec=DEGREE + 1)
            for coefficient in reversed(parent.coefficients[:-1]):
                poly = poly * affine + coefficient
            output.append(Panel(center, child_half, tuple(poly[j] for j in range(DEGREE + 1)),
                                parent.remainder, parent.parent_index))
    finally:
        ctx.cap = previous_cap
    if any(2*p.half_width > F(1, 64) for p in output):
        raise ValueError('Child phase domain exceeds original proof')
    return tuple(output)

def mixed_exact(j, m):
    return 2**(m+1) * sum((F(comb(j, ell)*(-2)**ell, m+ell+1)
                           for ell in range(j+1)), F(0))

def mixed_table():
    """I[j,m]=integral x^j(1-x)^m = I[j-1,m]-I[j-1,m+1]."""
    top = PHASE_DEGREE + DEGREE + 1
    rows = [[F(2**(m+1), m+1) for m in range(top+1)]]
    for j in range(1, DEGREE+1):
        previous = rows[-1]
        rows.append([previous[m]-previous[m+1] for m in range(len(previous)-1)])
    return rows

class Kernels:
    def __init__(self):
        if ctx.prec < PRECISION_BITS:
            raise ValueError('512-bit kernel arithmetic required')
        rows = mixed_table()
        self.exp = [acb_poly([ball(rows[j][m]/factorial(m)) for m in range(PHASE_DEGREE+1)])
                    for j in range(DEGREE+1)]
        self.drift = [acb_poly([ball(rows[j][m+1]/factorial(m+1)) for m in range(PHASE_DEGREE+1)])
                      for j in range(DEGREE+1)]
        self.e_tail = F(2*4096*8**(PHASE_DEGREE+1), factorial(PHASE_DEGREE+1))
        self.q_tail = F(4*4096*8**(PHASE_DEGREE+1), factorial(PHASE_DEGREE+2))
        self.cache = {}

    def at(self, k, h, include_tail=True):
        if not (0 <= k <= 256 and h > 0 and 2*k*h <= 4):
            raise ValueError('Outside proved entire-kernel phase domain')
        key = (k, h, include_tail)
        if key not in self.cache:
            z = acb(0, ball(2*k*h))
            e = tuple(poly(z) for poly in self.exp)
            q = tuple(poly(z) for poly in self.drift)
            if include_tail and k:
                e = tuple(inflate_complex(value, self.e_tail) for value in e)
                q = tuple(inflate_complex(value, self.q_tail) for value in q)
            self.cache[key] = e, q
        return self.cache[key]

def stable_drift(k, distance, include_tail=True):
    if k < 0 or distance < 0:
        raise ValueError('Invalid drift argument')
    z = acb(0, ball(2*k*distance))
    if 2*k*distance <= 1:
        poly = acb_poly([ball(F(1, factorial(m+1))) for m in range(PHASE_DEGREE+1)])
        tail = distance * F(3, factorial(PHASE_DEGREE+2)) if k and include_tail else F(0)
        return inflate_complex(ball(distance)*poly(z), tail)
    return (z.exp()-1) / acb(0, ball(2*k))

def moments(panel, kernels, k, model=True):
    e, q = kernels.at(k, panel.half_width, model)
    h = ball(panel.half_width)
    m0 = h * sum((a*ball(F(2,j+1) if j%2==0 else 0)
                  for j,a in enumerate(panel.coefficients)), arb(0))
    me = h * sum((a*b for a,b in zip(panel.coefficients,e)), acb(0))
    mu = h*h * sum((a*b for a,b in zip(panel.coefficients,q)), acb(0))
    r = panel.remainder if model else F(0)
    me = inflate_complex(me, 2*panel.half_width*r)
    mu = inflate_complex(mu, 2*panel.half_width**2*r)
    m0 = inflate_real(m0, 2*panel.half_width*r)
    if not k:
        me, mu = acb(me.real, arb(0)), acb(mu.real, arb(0))
    return m0, me, mu

def incoming(panels, k, family, model=True):
    panels = tuple(panels)
    if not panels or panels[0].center-panels[0].half_width != A or panels[-1].center+panels[-1].half_width != B:
        raise ValueError('Incomplete fixed prehistory interior')
    for left,right in zip(panels,panels[1:]):
        if left.center+left.half_width != right.center-right.half_width:
            raise ValueError('Noncontiguous source panels')
    mw,mu = acb(0),acb(0)
    for p in panels:
        distance = B-p.center-p.half_width
        rotation = acb(0,ball(2*k*distance)).exp()
        m0,me,mdrift = moments(p,family,k,model)
        mw += rotation*me
        mu += stable_drift(k,distance,model)*m0+rotation*mdrift
    return {'U': -mu, 'W': -mw}

def cap_bounds(source):
    d = 2*DELTA-DELTA*DELTA
    lam = 1/(5-DELTA)
    h = F(3,8)**63
    if source == 'positive_B':
        g = DELTA*h*(4*lam*lam+4*lam/d**2+4/d**4)
    elif source == 'signed_uB':
        g = DELTA*h*(4*lam*lam+2*lam+(4*lam+4)/d**2+4/d**4)
    else:
        raise ValueError('Unknown source')
    return {'U':g/2,'W':g}
