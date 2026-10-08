"""Validated later endpoint flow. Caller supplies authenticated source polynomials.

No I/O, physical-source construction, numerical-array reader, or implicit float.
The finite polynomial source response is integrated by exact moment identities;
Arb contains arithmetic rounding, with an explicit analytic kernel tail.
"""
from fractions import Fraction as Q
from math import comb, factorial
from flint import arb, acb, acb_poly, fmpq, ctx

BITS = 1024
DEGREE = 2048
SOURCE_DEGREE = 24
EXPORT_BITS = 96
SCALE = 1 << EXPORT_BITS
EXPORT_GATE = Q(1, 10**18)
ANCHOR = Q(-9, 2)
ENDPOINTS = (Q(-4), Q(-7, 2))
MAX_K = Q(256)
HALF = Q(1, 128)

def require(ok, text):
    if not ok:
        raise ValueError(text)

def av(q):
    require(type(q) is Q, 'exact Fraction required')
    return arb(fmpq(q.numerator, q.denominator))

def cq(z):
    require(type(z) is tuple and len(z) == 2 and all(type(v) is Q for v in z),
            'exact complex pair required')
    return acb(av(z[0]), av(z[1]))

def exact_q(a):
    q = a.fmpq()
    return Q(int(q.numerator), int(q.denominator))

def floor_grid(q):
    return q.numerator*SCALE//q.denominator

def ceil_grid(q):
    return -floor_grid(-q)

def upward_tail(q):
    # A fixed outward dyadic bound keeps enormous factorial denominators out of
    # exported JSON; this arithmetic enlargement is included once in the tail.
    scale=1<<512
    return Q((q.numerator*scale+q.denominator-1)//q.denominator,scale)

def rectangle(z, error=Q(0)):
    require(type(error) is Q and error >= 0, 'nonnegative exact error required')
    result = {}
    for axis, part in (('real', z.real), ('imag', z.imag)):
        require(part.is_finite(), 'nonfinite endpoint ball')
        result[axis] = [floor_grid(exact_q(part.lower())-error),
                        ceil_grid(exact_q(part.upper())+error)]
    radius = Q(sum(result[a][1]-result[a][0] for a in ('real', 'imag')), 2*SCALE)
    require(radius <= EXPORT_GATE, 'complete exported endpoint radius exceeds fixed gate')
    return result

def phase_coefficient(value, n):
    return (acb(value, 0), acb(0, value), acb(-value, 0), acb(0, -value))[n % 4]

def exponential_tail(x, n):
    require(type(x) is Q and x >= 0 and type(n) is int and n >= 0 and x < n+2,
            'invalid factorial-tail premise')
    return x**(n+1)/factorial(n+1)/(1-x/Q(n+2))

def geometry():
    return tuple({'left':ANCHOR+Q(j,64), 'right':ANCHOR+Q(j+1,64),
                  'center':ANCHOR+Q(2*j+1,128), 'half':HALF} for j in range(64))

def polynomial_in_lag(coefficients, shift):
    """g(s)=p(s-c), s-c=shift-y/2; all transformation arithmetic exact."""
    out = [Q(0)]*len(coefficients)
    for m, coefficient in enumerate(coefficients):
        for j in range(m+1):
            out[j] += coefficient*comb(m,j)*shift**(m-j)*Q(-1,2)**j
    return tuple(out)

class EndpointFlow:
    """Exact finite-polynomial Duhamel response plus certified uniform errors.

    There are exactly64 source rows, ordinary powers in s-center, each25 exact
    coefficients. 'uniform_error' encloses actual-g minus chosen-p on the cell.
    Caller must establish that premise independently and authenticate the rows.
    """
    def __init__(self, rows, endpoint, check=lambda:None):
        ctx.prec = BITS
        require(endpoint in ENDPOINTS, 'unregistered endpoint')
        require(type(rows) is list and len(rows) == 64, 'complete64cell source required')
        self.endpoint = endpoint
        self.duration = endpoint-ANCHOR
        self.source_l1 = Q(0)
        self.source_error_W = Q(0)
        self.source_error_U = Q(0)
        moments = [arb(0) for _ in range(DEGREE+2)]
        for index, (row, g) in enumerate(zip(rows, geometry())):
            require(row['center'] == g['center'] and row['half'] == g['half'],
                    'source geometry changed')
            co = row['coefficients']; error = row['uniform_error']
            require(type(co) in (tuple,list) and len(co) == SOURCE_DEGREE+1 and
                    all(type(v) is Q for v in co), '25 exact source coefficients required')
            require(type(error) is Q and error >= 0, 'nonnegative source error required')
            if g['right'] > endpoint:
                continue
            check()
            self.source_l1 += 2*HALF*sum((abs(v)*HALF**j for j,v in enumerate(co)), Q(0))
            self.source_error_W += 2*HALF*error
            self.source_error_U += 2*HALF*(endpoint-g['center'])*error
            d = polynomial_in_lag(co, endpoint-g['center'])
            yl, yr = av(2*(endpoint-g['left'])), av(2*(endpoint-g['right']))
            lp, rp = [arb(1)], [arb(1)]
            for n in range(DEGREE+SOURCE_DEGREE+3):
                lp.append(lp[-1]*yl); rp.append(rp[-1]*yr)
            terms = [(j,av(v)/2) for j,v in enumerate(d) if v]
            for n in range(DEGREE+2):
                total = arb(0)
                for j, v in terms:
                    total += v*(lp[n+j+1]-rp[n+j+1])/(n+j+1)
                moments[n] += total
        w, u, fact = [], [], 1
        for n in range(DEGREE+1):
            if n:
                fact *= n
            w.append(phase_coefficient(-moments[n]/fact, n))
            u.append(phase_coefficient(-moments[n+1]/(2*(n+1)*fact), n))
        self.W = acb_poly(w)
        self.U = acb_poly(u)
        x = 2*MAX_K*self.duration
        tail = exponential_tail(x, DEGREE)
        self.kernel_error_W = upward_tail(self.source_l1*tail)
        self.kernel_error_U = upward_tail(self.source_l1*self.duration*tail)

    def evaluate(self, k, incoming_u, incoming_w):
        require(ctx.prec == BITS, 'arithmetic precision changed')
        require(type(k) is Q and 0 <= k <= MAX_K, 'exact k outside registered band')
        x = av(k)
        t = av(self.duration)
        E = acb(0, 2*x*t).exp()
        # Entire form at k=0, with no reciprocal momentum or exp cancellation.
        Phi = t*acb(0, x*t).exp()*(x*t).sinc()
        W = E*cq(incoming_w)+self.W(x)
        U = cq(incoming_u)+Phi*cq(incoming_w)+self.U(x)
        return {'U':rectangle(U, self.source_error_U+self.kernel_error_U),
                'W':rectangle(W, self.source_error_W+self.kernel_error_W)}

    def budget(self):
        return {'endpoint':str(self.endpoint), 'source_polynomial_L1_upper':str(self.source_l1),
                'source_error_U':str(self.source_error_U), 'source_error_W':str(self.source_error_W),
                'kernel_tail_U':str(self.kernel_error_U), 'kernel_tail_W':str(self.kernel_error_W),
                'arithmetic':'Arb1024bits contained in exported endpoint rectangles',
                'integration':'finite polynomial moments integrated exactly; degree2048 exponential tail rounded outward to dyadic512'}
