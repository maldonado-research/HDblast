#!/usr/bin/env python3
"""Exact independent second-Taylor checker from diagonal metric Christoffels.

Only the Python standard library is required.  Every operation uses Fraction;
the fields are arbitrary two-jets and no background equation is imposed.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def q(value):
    return value if isinstance(value, F) else F(value)


@dataclass(frozen=True)
class P:
    """Truncated epsilon polynomial through epsilon squared."""
    a: tuple

    def __init__(self, a=0, b=0, c=0):
        object.__setattr__(self, "a", (q(a), q(b), q(c)))

    def __add__(self, other):
        other = promote(other)
        return P(*(x+y for x, y in zip(self.a, other.a)))

    __radd__ = __add__

    def __neg__(self):
        return P(*(-x for x in self.a))

    def __sub__(self, other):
        return self + (-promote(other))

    def __rsub__(self, other):
        return promote(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, Jet):
            return NotImplemented
        other = promote(other)
        return P(*(sum(self.a[j]*other.a[i-j] for j in range(i+1)) for i in range(3)))

    __rmul__ = __mul__

    def inverse(self):
        a, b, c = self.a
        if not a:
            raise ValueError("Cannot invert an epsilon polynomial with zero constant term")
        return P(1/a, -b/a**2, b**2/a**3-c/a**2)

    def __truediv__(self, other):
        return self * promote(other).inverse()

    def __rtruediv__(self, other):
        return promote(other)*self.inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse()**(-n)
        out = P(1)
        for _ in range(n):
            out *= self
        return out


def promote(x):
    return x if isinstance(x, P) else P(x)


ZERO = P()
EPS = P(0, 1)
DIM = 5


class Jet:
    """Coordinate two-jet with epsilon-polynomial entries."""
    def __init__(self, value, gradient=None, hessian=None):
        self.value = promote(value)
        self.d = [ZERO for _ in range(DIM)] if gradient is None else list(map(promote, gradient))
        self.dd = ([[ZERO for _ in range(DIM)] for _ in range(DIM)]
                   if hessian is None else [list(map(promote, row)) for row in hessian])

    def __add__(self, other):
        other = asjet(other)
        return Jet(self.value+other.value,
                   [self.d[i]+other.d[i] for i in range(DIM)],
                   [[self.dd[i][j]+other.dd[i][j] for j in range(DIM)] for i in range(DIM)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, [-d for d in self.d], [[-d for d in row] for row in self.dd])

    def __sub__(self, other):
        return self + (-asjet(other))

    def __rsub__(self, other):
        return asjet(other) + (-self)

    def __mul__(self, other):
        other = asjet(other)
        return Jet(self.value*other.value,
                   [self.d[i]*other.value+self.value*other.d[i] for i in range(DIM)],
                   [[self.dd[i][j]*other.value+self.d[i]*other.d[j]
                     +self.d[j]*other.d[i]+self.value*other.dd[i][j]
                     for j in range(DIM)] for i in range(DIM)])

    __rmul__ = __mul__

    def inverse(self):
        iv = self.value.inverse()
        return Jet(iv, [-iv**2*d for d in self.d],
                   [[2*iv**3*self.d[i]*self.d[j]-iv**2*self.dd[i][j]
                     for j in range(DIM)] for i in range(DIM)])

    def __truediv__(self, other):
        return self*asjet(other).inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse()**(-n)
        out = Jet(1)
        for _ in range(n):
            out *= self
        return out


def asjet(value):
    return value if isinstance(value, Jet) else Jet(value)


def field_jet(value, dy=0, dt=0, dyy=0, dyt=0, dtt=0):
    d = [dy, dt, 0, 0, 0]
    dd = [[0 for _ in range(DIM)] for _ in range(DIM)]
    dd[0][0], dd[0][1], dd[1][0], dd[1][1] = dyy, dyt, dyt, dtt
    return Jet(value, d, dd)


def angular_data(st, ct, sf, cf):
    """Hopf metric, Y and all ordinary coordinate derivatives, exactly."""
    st, ct, sf, cf = map(q, (st, ct, sf, cf))
    if st*st+ct*ct != 1 or sf*sf+cf*cf != 1 or not st or not ct:
        raise ValueError("Invalid rational Hopf point")
    x, w, C, S = st*st, st*ct, cf*cf-sf*sf, 2*sf*cf
    xd = [0, 0, 2*w, 0, 0]
    xdd = [[0 for _ in range(DIM)] for _ in range(DIM)]
    xdd[2][2] = 2*(1-2*x)
    xj = Jet(x, xd, xdd)
    yd = [0, 0, 2*w*C, -2*x*S, 0]
    ydd = [[0 for _ in range(DIM)] for _ in range(DIM)]
    ydd[2][2], ydd[2][3], ydd[3][2], ydd[3][3] = 2*(1-2*x)*C, -4*w*S, -4*w*S, -4*x*C
    yj = Jet(x*C, yd, ydd)
    k = [F(1), x, 1-x]
    Y = x*C
    H = x-F(1, 2)
    Q = Y*Y-x/2+F(1, 12)
    dY = [2*w*C, -2*x*S, F(0)]
    dH = [2*w, F(0), F(0)]
    dQ = [2*Y*dY[i]-dH[i]/2 for i in range(3)]
    hY = [[2*(1-2*x)*C, -2*w*S, F(0)],
          [-2*w*S, -2*x*(1+x)*C, F(0)],
          [F(0), F(0), -2*x*(1-x)*C]]
    gradnorm = sum(dY[i]**2/k[i] for i in range(3))
    A = [[Y*(hY[i][j]+(F(8, 3)*Y*k[i] if i == j else 0))
          for j in range(3)] for i in range(3)]
    B = [[dY[i]*dY[j]-(gradnorm*k[i]/3 if i == j else 0)
          for j in range(3)] for i in range(3)]
    hx = [2*(1-2*x), 2*x*(1-x), -2*x*(1-x)]
    hQTF = [[2*Y*hY[i][j]+2*dY[i]*dY[j]
             +(-hx[i]/2+8*Q*k[i] if i == j else 0)
             for j in range(3)] for i in range(3)]
    return {"x": xj, "Yjet": yj, "Qjet": yj*yj-xj/2+F(1, 12),
            "k": k, "H": H, "Q": Q, "hQTF": hQTF,
            "dH": dH, "dQ": dQ, "A": A, "B": B}


def direct_residual(case, point, control="none"):
    """Build Gamma, dGamma, Ricci and Box directly in coordinates."""
    d = angular_data(*point)
    u = field_jet(*(case[key] for key in ("u", "uy", "ut", "uyy", "uyt", "utt")))
    v = field_jet(*(case[key] for key in ("v", "vy", "vt", "vyy", "vyt", "vtt")))
    rho = field_jet(case["rho"], case["rhoy"], dyy=case["rhoyy"])
    phi0 = field_jet(case["phi0"], case["p"], dyy=case["py"])
    c, s = case["c"], case["s"]
    c2 = field_jet(c*c, dt=2*c*s, dtt=2*(c*c+s*s))
    m = 1+2*EPS*u*d["Yjet"]
    g = [1-4*EPS*u*d["Yjet"], -rho**2*m, rho**2*c2*m,
         rho**2*c2*m*d["x"], rho**2*c2*m*(1-d["x"])]
    iv = [g_i.inverse() for g_i in g]
    Gamma = [[[ZERO for _ in range(DIM)] for _ in range(DIM)] for _ in range(DIM)]
    dGamma = [[[[ZERO for _ in range(DIM)] for _ in range(DIM)]
               for _ in range(DIM)] for _ in range(DIM)]
    for a in range(DIM):
        for b in range(DIM):
            for cidx in range(DIM):
                numerator = ((g[a].d[b] if a == cidx else ZERO)
                             +(g[a].d[cidx] if a == b else ZERO)
                             -(g[b].d[a] if b == cidx else ZERO))
                Gamma[a][b][cidx] = iv[a].value*numerator/2
                for j in range(DIM):
                    derivative = ((g[a].dd[b][j] if a == cidx else ZERO)
                                  +(g[a].dd[cidx][j] if a == b else ZERO)
                                  -(g[b].dd[a][j] if b == cidx else ZERO))
                    dGamma[j][a][b][cidx] = (iv[a].d[j]*numerator+iv[a].value*derivative)/2
    Ricci = [[ZERO for _ in range(DIM)] for _ in range(DIM)]
    for b in range(DIM):
        for cidx in range(DIM):
            R = sum((dGamma[a][a][b][cidx]-dGamma[cidx][a][b][a]
                     for a in range(DIM)), ZERO)
            if control != "drop-ricci-products":
                R += sum((Gamma[a][b][cidx]*Gamma[j][a][j]
                          -Gamma[a][b][j]*Gamma[j][cidx][a]
                          for a in range(DIM) for j in range(DIM)), ZERO)
            Ricci[b][cidx] = -R if control == "flip-ricci-sign" else R
    phi = phi0+EPS*v*d["Yjet"]
    delta_phi = EPS*v.value*d["Yjet"].value
    U = case["U0"]+case["U1"]*delta_phi+case["U2"]*delta_phi**2/2
    Uprime = case["U1"]+case["U2"]*delta_phi+case["U3"]*delta_phi**2/2
    potential_factor = F(1, 3) if control == "wrong-potential-normalization" else F(2, 3)
    E = [[Ricci[i][j]-(ZERO if control == "drop-gradient-square" else phi.d[i]*phi.d[j])
          -(potential_factor*U*g[i].value if i == j else ZERO)
          for j in range(DIM)] for i in range(DIM)]
    box = sum((iv[i].value*(phi.dd[i][i]-sum((Gamma[j][i][i]*phi.d[j]
               for j in range(DIM)), ZERO)) for i in range(DIM)), ZERO)
    factor = 1 if control == "omit-d2-factor" else 2
    E2 = [[factor*E[i][j].a[2] for j in range(DIM)] for i in range(DIM)]
    S2 = factor*(box-Uprime).a[2]
    tr = sum(E2[i+2][i+2]/d["k"][i] for i in range(3))
    TF = [[E2[i+2][j+2]-(tr*d["k"][i]/3 if i == j else 0)
           for j in range(3)] for i in range(3)]
    return d, {"yy": E2[0][0], "yt": E2[0][1], "tt": E2[1][1],
               "angular_trace": tr/(3*case["c"]**2), "scalar": S2,
               "ya": E2[0][2:], "ta": E2[1][2:], "TF": TF}


def solve(rows, values):
    """Exact overdetermined Gaussian solve; inconsistency raises explicitly."""
    if not rows or len(rows) != len(values):
        raise ValueError("Invalid projection system")
    width = len(rows[0])
    M = [[q(x) for x in row]+[q(y)] for row, y in zip(rows, values)]
    pivots, outrow = [], 0
    for col in range(width):
        pivot = next((i for i in range(outrow, len(M)) if M[i][col]), None)
        if pivot is None:
            continue
        M[outrow], M[pivot] = M[pivot], M[outrow]
        scale = M[outrow][col]
        M[outrow] = [x/scale for x in M[outrow]]
        for i in range(len(M)):
            if i != outrow and M[i][col]:
                scale = M[i][col]
                M[i] = [a-scale*b for a, b in zip(M[i], M[outrow])]
        pivots.append(col)
        outrow += 1
    if any(not any(row[:width]) and row[-1] for row in M):
        raise ArithmeticError("Direct angular samples do not lie in the specified harmonic span")
    if len(pivots) != width:
        raise ArithmeticError("Angular projection system is rank deficient")
    result = [F(0) for _ in range(width)]
    for row, col in enumerate(pivots):
        result[col] = M[row][-1]
    return result


POINTS = [
    (F(3, 5), F(4, 5), F(0), F(1)),
    (F(3, 5), F(4, 5), F(3, 5), F(4, 5)),
    (F(5, 13), F(12, 13), F(5, 13), F(12, 13)),
    (F(8, 17), F(15, 17), F(8, 17), F(15, 17)),
    (F(4, 5), F(3, 5), F(12, 13), F(5, 13)),
]


def angular_projection_checks():
    """Obtain the projection factors from Laplace eigenvalues and integration.

    The normalized S3 measure for these psi-independent functions is uniform
    in x=sin(theta)^2 and in phi.  Integral Q^2 uses <cos^2>=1/2 and
    <cos^4>=3/8.  Two integrations by parts give the B:T contraction below.
    """
    normQ = F(3, 40)+F(1, 12)+F(1, 144)-F(1, 8)+F(1, 36)-F(1, 24)
    lambdaY, lambdaQ = F(8), F(24)
    gradient_projection = lambdaY-lambdaQ/2
    hessian_norm_factor = lambdaQ*(lambdaQ-2)-lambdaQ**2/3
    B_contraction_factor = lambdaY*lambdaQ/2-lambdaQ*gradient_projection/2+lambdaQ*gradient_projection/3
    cB = B_contraction_factor/hessian_norm_factor
    cA = F(1, 2)-cB
    if normQ != F(19, 720) or cA != F(1, 6) or cB != F(1, 3):
        raise ArithmeticError("Independent angular integration has failed")
    for point in POINTS:
        data = angular_data(*point)
        traceTF = sum(data["hQTF"][i][i]/data["k"][i] for i in range(3))
        if traceTF:
            raise ArithmeticError("Direct angular Hessian is inconsistent with Delta Q = -24 Q")
    return {"Q_norm_squared_normalized_measure": normQ,
            "HessQ_TF_norm_squared_normalized_measure": hessian_norm_factor*normQ,
            "Y_HessY_TF_projection": cA, "DY_DY_TF_projection": cB,
            "gradient_square_Q_projection": gradient_projection}


def project(case, control="none"):
    angular, residuals = zip(*(direct_residual(case, point, control) for point in POINTS))
    rows = [[1, d["H"], d["Q"]] for d in angular]
    result = {key: solve(rows, [r[key] for r in residuals])[2]
              for key in ("yy", "yt", "tt", "angular_trace", "scalar")}
    for key in ("ya", "ta"):
        rows, values = [], []
        for d, r in zip(angular, residuals):
            for i in range(3):
                rows.append([d["dH"][i], d["dQ"][i]])
                values.append(r[key][i])
        result[key] = solve(rows, values)[1]
    rows, values = [], []
    for d, r in zip(angular, residuals):
        for i in range(3):
            for j in range(i, 3):
                rows.append([d["A"][i][j], d["B"][i][j]])
                values.append(r["TF"][i][j])
    result["TF_A"], result["TF_B"] = solve(rows, values)
    angular_factors = angular_projection_checks()
    cA, cB = angular_factors["Y_HessY_TF_projection"], angular_factors["DY_DY_TF_projection"]
    if control == "wrong-angular-projection":
        cA *= 2
    result["angular_TF"] = cA*result["TF_A"]+cB*result["TF_B"]
    return result


def expected(case):
    """Candidate closed-form coefficients, compared with direct geometry."""
    u, uy, ut, uyy, uyt, utt = (case[key] for key in ("u", "uy", "ut", "uyy", "uyt", "utt"))
    v, vy, vt = (case[key] for key in ("v", "vy", "vt"))
    rho, rhoy, rhoyy, c, s = (case[key] for key in ("rho", "rhoy", "rhoyy", "c", "s"))
    a, ay, b = rhoy/rho, rhoyy/rho-(rhoy/rho)**2, s/c
    U1, U2 = case["U1"], case["U2"]
    K = 8*uy**2+8*u*uyy+104*a*u*uy+48*(ay+4*a*a)*u*u+F(2, 3)*U2*v*v+F(8, 3)*U1*u*v
    A = -4*u*utt-12*b*u*ut-4*ut**2-48*u*u/(c*c)-rho*rho*K
    return {
        "scalar": 12*uy*vy+8*u*case["vyy"]+32*a*u*vy+48*case["p"]*u*uy
                  +32*u*u*(case["py"]+4*a*case["p"])
                  +4*u*(case["vtt"]+3*b*vt+8*v/(c*c))/(rho*rho)-case["U3"]*v*v,
        "yy": -8*uy**2+16*u*uyy-32*a*u*uy-2*vy**2-F(2, 3)*U2*v*v+F(16, 3)*U1*u*v
              +(8*u*utt+24*b*u*ut-16*ut**2)/(rho*rho),
        "yt": 12*u*uyt-48*a*u*ut-2*vy*vt,
        "tt": 24*u*utt+12*ut**2-2*vt**2-A,
        "ya": 6*u*uy-24*a*u*u-v*vy,
        "ta": 18*u*ut-12*b*u*u-v*vt,
        "angular_trace": A-24*b*u*ut+(-80*u*u+F(8, 3)*v*v)/(c*c),
        "angular_TF": 8*u*u-F(2, 3)*v*v,
        "TF_A": 24*u*u,
        "TF_B": 12*u*u-2*v*v,
    }


def derivative_jet(jet, coordinate):
    """A first derivative with sufficient gradient data for this operator."""
    return Jet(jet.d[coordinate], jet.dd[coordinate])


def base_geometry(case):
    rho = field_jet(case["rho"], case["rhoy"], dyy=case["rhoyy"])
    c = field_jet(case["c"], dt=case["s"], dtt=case["c"])
    R = rho*c
    qb = [[Jet(1), Jet(0)], [Jet(0), -rho*rho]]
    qi = [[Jet(1), Jet(0)], [Jet(0), -(rho**-2)]]
    gamma = [[[Jet(0) for _ in range(2)] for _ in range(2)] for _ in range(2)]
    for a in range(2):
        for b in range(2):
            for cidx in range(2):
                gamma[a][b][cidx] = sum((qi[a][j]*(derivative_jet(qb[j][cidx], b)
                    +derivative_jet(qb[j][b], cidx)-derivative_jet(qb[b][cidx], j))/2
                    for j in range(2)), Jet(0))
    r = [derivative_jet(R, a)/R for a in range(2)]
    return rho, R, qb, qi, gamma, r


def linear_amplitudes(index):
    """Five independently assigned arbitrary amplitude two-jets."""
    result = {}
    for n, name in enumerate(("fyy", "fyt", "ftt", "H", "w")):
        values = [F((-1)**(index+n+j)*(11+(7*index+5*n+3*j)%19),
                    17+(index+2*n+5*j)%13) for j in range(6)]
        result[name] = field_jet(*values)
    return result


def expected_linear(case, amp, control="none"):
    """Evaluate the supplied 2D covariant formulas, separately from 5D Ricci."""
    rho, R, qb, qi, gamma, r = base_geometry(case)
    f = [[amp["fyy"], amp["fyt"]], [amp["fyt"], amp["ftt"]]]
    H, w = amp["H"], amp["w"]
    trace = sum((qi[a][b]*f[a][b] for a in range(2) for b in range(2)), Jet(0))
    T = trace/2+3*H
    mixed = [[sum((qi[a][c]*f[c][b] for c in range(2)), Jet(0))
              for b in range(2)] for a in range(2)]
    upper = [[sum((qi[a][c]*qi[b][d]*f[c][d] for c in range(2) for d in range(2)), Jet(0))
              for b in range(2)] for a in range(2)]
    rupper = [sum((qi[a][b]*r[b] for b in range(2)), Jet(0)) for a in range(2)]
    Df = [[[derivative_jet(f[b][c], a)-sum((gamma[d][a][b]*f[d][c]
                  +gamma[d][a][c]*f[b][d] for d in range(2)), Jet(0))
            for c in range(2)] for b in range(2)] for a in range(2)]
    J = [[[sum((qi[a][d]*(Df[b][c][d]+Df[c][b][d]-Df[d][b][c])/2
                for d in range(2)), Jet(0)) for c in range(2)] for b in range(2)] for a in range(2)]
    divmixed = [sum((mixed[b][a].d[b]+sum((gamma[b][b][d].value*mixed[d][a].value
                          -gamma[d][b][a].value*mixed[b][d].value for d in range(2)), ZERO)
                         for b in range(2)), ZERO) for a in range(2)]
    P = [case["p"], F(0)]
    Pupper = [sum(qi[a][b].value.a[0]*P[b] for b in range(2)) for a in range(2)]
    pderivative = [[(case["py"] if a == b == 0 else F(0))
                   -sum(gamma[c][a][b].value.a[0]*P[c] for c in range(2))
                   for b in range(2)] for a in range(2)]
    Ebase = [[ZERO for _ in range(2)] for _ in range(2)]
    for a in range(2):
        for b in range(2):
            divJ = sum((J[c][a][b].d[c]
                +sum((gamma[c][c][d].value*J[d][a][b].value
                      -gamma[d][c][a].value*J[c][d][b].value
                      -gamma[d][c][b].value*J[c][a][d].value for d in range(2)), ZERO)
                for c in range(2)), ZERO)
            warpJ = 3*sum((r[c].value*J[c][a][b].value for c in range(2)), ZERO)
            hessT = T.dd[b][a]-sum((gamma[c][b][a].value*T.d[c] for c in range(2)), ZERO)
            Ebase[a][b] = (divJ+(ZERO if control == "linear-omit-base-warp" else warpJ)
                +12*f[a][b].value/R.value**2-3*(r[a].value*H.d[b]+r[b].value*H.d[a])-hessT
                -P[a]*w.d[b]-P[b]*w.d[a]
                -F(2, 3)*(case["U0"]*f[a][b].value+case["U1"]*w.value*qb[a][b].value))
    Emixed = []
    for a in range(2):
        warp = (sum((rupper[b].value*f[a][b].value for b in range(2)), ZERO)+r[a].value*trace.value)
        Emixed.append((divmixed[a]-trace.d[a]+(ZERO if control == "linear-omit-mixed-warp" else warp))/2
                      -2*H.d[a]-P[a]*w.value)
    Z = [sum((r[b]*upper[a][b] for b in range(2)), Jet(0))-2*rupper[a]*H
         -sum((qi[a][b]*derivative_jet(H, b) for b in range(2)), Jet(0)) for a in range(2)]
    divZ = sum((Z[a].d[a]+sum((gamma[a][a][b].value*Z[b].value for b in range(2)), ZERO)
                for a in range(2)), ZERO)
    B = -trace.value/2-H.value
    L = (divZ+3*sum((r[a].value*Z[a].value for a in range(2)), ZERO)
         +2*sum((rupper[a].value*H.d[a] for a in range(2)), ZERO)
         -sum((rupper[a].value*T.d[a] for a in range(2)), ZERO)
         +24*H.value/R.value**2-F(4, 3)*case["U0"]*H.value
         -(ZERO if control == "linear-omit-angular-potential" else F(2, 3)*case["U1"]*w.value))
    scalar = (sum((qi[a][b].value*(w.dd[a][b]-sum((gamma[c][a][b].value*w.d[c]
                         for c in range(2)), ZERO)) for a in range(2) for b in range(2)), ZERO)
              +3*sum((rupper[a].value*w.d[a] for a in range(2)), ZERO)
              -(24/R.value**2+case["U2"])*w.value
              -sum((upper[a][b].value*pderivative[a][b] for a in range(2) for b in range(2)), ZERO)
              -sum((Pupper[b]*(divmixed[b]+3*sum((r[a].value*mixed[a][b].value for a in range(2)), ZERO)
                                 -trace.d[b]/2-3*H.d[b]) for b in range(2)), ZERO))
    result = {"yy": Ebase[0][0], "yt": Ebase[0][1], "tt": Ebase[1][1],
              "ya": Emixed[0], "ta": Emixed[1], "angular_trace": L-8*B/R.value**2,
              "angular_TF": B, "scalar": scalar}
    return {key: value.a[0] for key, value in result.items()}


def direct_linear(case, amp, point):
    """Full 5D Christoffel calculation with an off-diagonal y,t block."""
    d = angular_data(*point)
    rho, R, _, _, _, _ = base_geometry(case)
    Qj = d["Qjet"]
    g = [[Jet(0) for _ in range(DIM)] for _ in range(DIM)]
    g[0][0], g[0][1], g[1][0], g[1][1] = (1+EPS*amp["fyy"]*Qj,
        EPS*amp["fyt"]*Qj, EPS*amp["fyt"]*Qj, -rho*rho+EPS*amp["ftt"]*Qj)
    angular_factor = R*R*(1+2*EPS*amp["H"]*Qj)
    g[2][2], g[3][3], g[4][4] = angular_factor, angular_factor*d["x"], angular_factor*(1-d["x"])
    iv = [[Jet(0) for _ in range(DIM)] for _ in range(DIM)]
    determinant = g[0][0]*g[1][1]-g[0][1]*g[1][0]
    iv[0][0], iv[0][1], iv[1][0], iv[1][1] = (g[1][1]/determinant,
        -g[0][1]/determinant, -g[1][0]/determinant, g[0][0]/determinant)
    for a in range(2, DIM):
        iv[a][a] = g[a][a].inverse()
    support = [(a,b) for a in range(DIM) for b in range(DIM) if a == b or (a < 2 and b < 2)]
    Gamma = [[[ZERO for _ in range(DIM)] for _ in range(DIM)] for _ in range(DIM)]
    dGamma = [[[[ZERO for _ in range(DIM)] for _ in range(DIM)]
               for _ in range(DIM)] for _ in range(DIM)]
    for a in range(DIM):
        inverse_columns = [j for i,j in support if i == a]
        for b in range(DIM):
            for cidx in range(DIM):
                for j in inverse_columns:
                    n = g[j][cidx].d[b]+g[j][b].d[cidx]-g[b][cidx].d[j]
                    Gamma[a][b][cidx] += iv[a][j].value*n/2
                    for k in range(DIM):
                        dn = g[j][cidx].dd[b][k]+g[j][b].dd[cidx][k]-g[b][cidx].dd[j][k]
                        dGamma[k][a][b][cidx] += (iv[a][j].d[k]*n+iv[a][j].value*dn)/2
    phi = field_jet(case["phi0"], case["p"], dyy=case["py"])+EPS*amp["w"]*Qj
    deltaphi = EPS*amp["w"].value*Qj.value
    U = case["U0"]+case["U1"]*deltaphi+case["U2"]*deltaphi**2/2
    Uprime = case["U1"]+case["U2"]*deltaphi+case["U3"]*deltaphi**2/2
    E = [[ZERO for _ in range(DIM)] for _ in range(DIM)]
    for b in range(DIM):
        for cidx in range(DIM):
            Ricci = sum((dGamma[a][a][b][cidx]-dGamma[cidx][a][b][a]
                +sum((Gamma[a][b][cidx]*Gamma[j][a][j]-Gamma[a][b][j]*Gamma[j][cidx][a]
                     for j in range(DIM)), ZERO) for a in range(DIM)), ZERO)
            E[b][cidx] = Ricci-phi.d[b]*phi.d[cidx]-F(2, 3)*U*g[b][cidx].value
    box = sum((iv[a][b].value*(phi.dd[a][b]-sum((Gamma[c][a][b]*phi.d[c]
                    for c in range(DIM)), ZERO)) for a,b in support), ZERO)
    E1 = [[entry.a[1] for entry in row] for row in E]
    tr = sum(E1[i+2][i+2]/d["k"][i] for i in range(3))
    TF = [[E1[i+2][j+2]-(tr*d["k"][i]/3 if i == j else 0) for j in range(3)] for i in range(3)]
    return d, {"yy": E1[0][0], "yt": E1[0][1], "tt": E1[1][1],
               "ya": E1[0][2:], "ta": E1[1][2:], "angular_trace": tr/(3*R.value.a[0]**2),
               "TF": TF, "scalar": (box-Uprime).a[1]}


def project_linear(case, amp):
    angular, residuals = zip(*(direct_linear(case, amp, point) for point in POINTS[:3]))
    result = {key: solve([[d["Q"]] for d in angular], [r[key] for r in residuals])[0]
              for key in ("yy", "yt", "tt", "angular_trace", "scalar")}
    for key in ("ya", "ta"):
        result[key] = solve([[d["dQ"][i]] for d in angular for i in range(3)],
                            [r[key][i] for r in residuals for i in range(3)])[0]
    result["angular_TF"] = solve([[d["hQTF"][i][j]] for d in angular for i in range(3) for j in range(i, 3)],
                                  [r["TF"][i][j] for r in residuals for i in range(3) for j in range(i, 3)])[0]
    return result


LINEAR_CONTROLS = ["linear-omit-base-warp", "linear-omit-mixed-warp", "linear-omit-angular-potential"]


def verify_linear(control="none"):
    checks = []
    for index, case in enumerate(make_cases()):
        amp = linear_amplitudes(index)
        actual, candidate = project_linear(case, amp), expected_linear(case, amp, control)
        differences = {key: str(actual[key]-value) for key, value in candidate.items() if actual[key] != value}
        if differences:
            raise VerificationFailure(json.dumps({"case": case["name"], "differences": differences}, sort_keys=True))
        checks.append({"case": case["name"], "projected": {k: str(v) for k,v in actual.items()},
                       "differences": differences, "component_relations": component_relation_checks(case, amp, actual)})
    return checks


def component_relation_checks(case, amp, verified_operator):
    """Check note (16), (17), (18) against the verified covariant operator."""
    rho = field_jet(case["rho"], case["rhoy"], dyy=case["rhoyy"])
    a, b, p = case["rhoy"]/case["rho"], case["s"]/case["c"], case["p"]
    N, B, T, K = amp["fyy"], amp["fyt"], amp["ftt"]/(rho*rho), 2*amp["H"]
    value = lambda z: z.value.a[0]
    dy = lambda z: z.d[0].a[0]
    dt = lambda z: z.d[1].a[0]
    component_y = dy(T)/2+3*a*value(N)/2-(dt(B)+b*value(B))/(2*case["rho"]**2)-dy(K)-p*value(amp["w"])
    component_t = dy(B)/2-dt(N)/2+a*value(B)+b*value(N)/2-b*value(T)-dt(K)
    differences = {"equation_16_y": component_y-verified_operator["ya"],
                   "equation_16_t": component_t-verified_operator["ta"]}
    u = field_jet(*(case[key] for key in ("u", "uy", "ut", "uyy", "uyt", "utt")))
    v = field_jet(*(case[key] for key in ("v", "vy", "vt", "vyy", "vyt", "vtt")))
    d = 16*u*u-F(4, 3)*v*v
    G = amp["ftt"]
    forced_amplitudes = {"fyy": d+G-2*K, "fyt": B, "ftt": rho*rho*(G-K), "H": K/2,
                         "w": amp["w"]}
    forced_operator = expected_linear(case, forced_amplitudes)
    forcing = expected(case)
    lhs17 = dy(B)+2*a*value(B)-dt(G)-b*value(G)
    rhs17 = -4*case["u"]*case["ut"]-F(2, 3)*case["v"]*case["vt"]+8*b*case["u"]**2+F(4, 3)*b*case["v"]**2
    differences["equation_17_identity"] = 2*(forced_operator["ta"]+forcing["ta"])-(lhs17-rhs17)
    rhs18 = ((dy(G)+3*a*value(G))/2-F(3, 2)*(dy(K)+2*a*value(K))
             +F(3, 2)*a*value(d)-(dt(B)+b*value(B))/(2*case["rho"]**2)+forcing["ya"])
    differences["equation_18_identity"] = forced_operator["ya"]+forcing["ya"]-(rhs18-p*value(amp["w"]))
    if not p:
        raise ArithmeticError("Reconstruction check requires the stated collar assumption p != 0")
    forced_amplitudes["w"] = Jet(rhs18/p, amp["w"].d, amp["w"].dd)
    reconstructed_operator = expected_linear(case, forced_amplitudes)
    differences["equation_18_forced_y_residual"] = reconstructed_operator["ya"]+forcing["ya"]
    if any(differences.values()):
        raise VerificationFailure(json.dumps({"case": case["name"], "component_differences": {
            key: str(value) for key,value in differences.items() if value}}, sort_keys=True))
    return {key: str(value) for key,value in differences.items()}


def make_cases():
    fields = ["rho", "rhoy", "rhoyy", "phi0", "p", "py", "u", "uy", "ut", "uyy", "uyt", "utt",
              "v", "vy", "vt", "vyy", "vyt", "vtt", "U0", "U1", "U2", "U3"]
    raw = [
        [2, F(3, 7), F(-4, 9), F(1, 5), F(5, 11), F(-2, 3), F(2, 7), F(-3, 5), F(4, 9), F(5, 8), F(-7, 10), F(2, 11), F(-4, 13), F(5, 9), F(-2, 7), F(3, 10), F(7, 12), F(-4, 11), F(-2, 5), F(3, 8), F(-5, 7), F(8, 9)],
        [F(5, 3), F(-7, 8), F(9, 10), F(-2, 9), F(-3, 4), F(7, 13), F(-5, 11), F(7, 9), F(-8, 13), F(11, 15), F(4, 7), F(-9, 8), F(5, 12), F(-7, 10), F(6, 11), F(-8, 9), F(3, 8), F(5, 7), F(4, 9), F(-5, 12), F(7, 11), F(-2, 3)],
        [F(7, 4), F(2, 5), F(-3, 7), F(4, 11), F(6, 13), F(-9, 14), F(7, 15), F(3, 8), F(-4, 9), F(-2, 11), F(8, 13), F(5, 12), F(-3, 7), F(2, 9), F(7, 10), F(4, 13), F(-5, 11), F(8, 15), F(-6, 7), F(9, 10), F(-4, 5), F(3, 11)],
    ]
    times = [(F(5, 4), F(3, 4)), (F(13, 12), F(-5, 12)), (F(17, 15), F(8, 15))]
    cases = []
    for index, (values, (c, s)) in enumerate(zip(raw, times)):
        case = dict(zip(fields, map(q, values)))
        case.update(c=c, s=s, name=f"arbitrary_jet_{index+1}")
        if c*c-s*s != 1:
            raise ValueError("Hyperbolic point is invalid")
        cases.append(case)
    return cases


CONTROLS = ["flip-ricci-sign", "drop-ricci-products", "drop-gradient-square",
            "wrong-potential-normalization", "omit-d2-factor", "wrong-angular-projection"]


class VerificationFailure(RuntimeError):
    pass


def verify(control="none"):
    checks = []
    for case in make_cases():
        actual = project(case, control)
        candidate = expected(case)
        differences = {key: str(actual[key]-value) for key, value in candidate.items() if actual[key] != value}
        checks.append({"case": case["name"], "projected": {k: str(v) for k, v in actual.items()},
                       "differences": differences})
        if differences:
            raise VerificationFailure(json.dumps({"case": case["name"], "differences": differences}, sort_keys=True))
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write JSON evidence to this path")
    parser.add_argument("--control", choices=["none"]+CONTROLS+LINEAR_CONTROLS, default="none",
                        help="Deliberately corrupt one formula; successful rejection exits nonzero")
    parser.add_argument("--include-controls", action="store_true", help="Check that every deliberate corruption is rejected")
    args = parser.parse_args(argv)
    report = {"method": "exact Fraction epsilon and coordinate two-jets; direct Christoffels and Ricci",
              "epsilon_convention": "D2 residual = 2 times epsilon^2 coefficient",
              "python_optimized": not __debug__, "control": args.control, "angular_points": len(POINTS)}
    try:
        report["angular_projection"] = {key: str(value) for key, value in angular_projection_checks().items()}
        report["cases"] = verify(args.control if args.control in CONTROLS else "none")
        report["linear_cases"] = verify_linear(args.control if args.control in LINEAR_CONTROLS else "none")
        report["exact_coefficient_comparisons"] = {"quadratic": 30, "linear": 24, "component_and_forced_relations": 15}
        report["status"] = "PASS"
        if args.include_controls:
            report["negative_controls"] = {}
            for control in CONTROLS:
                try:
                    verify(control)
                except (VerificationFailure, ArithmeticError) as error:
                    report["negative_controls"][control] = {"status": "REJECTED", "reason": str(error)}
                else:
                    raise VerificationFailure(f"Negative control {control!r} unexpectedly passed")
            for control in LINEAR_CONTROLS:
                try:
                    verify_linear(control)
                except (VerificationFailure, ArithmeticError) as error:
                    report["negative_controls"][control] = {"status": "REJECTED", "reason": str(error)}
                else:
                    raise VerificationFailure(f"Negative control {control!r} unexpectedly passed")
    except (VerificationFailure, ArithmeticError) as error:
        report.update(status="FAIL", error=str(error))
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
