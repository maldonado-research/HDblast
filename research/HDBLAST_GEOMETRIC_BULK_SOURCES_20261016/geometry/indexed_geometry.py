"""Coordinate-independent quadratic geometry from covariant metric two-jets.

The scalar ring may be Fraction, a polynomial ring, or an interval ring.
All tensors are covariant. dq[d,a,b] = nabla_d q_ab and
ddq[e,d,a,b] = nabla_e nabla_d q_ab; the derivative slot d IS covariant.
The inverse metric is diagonal in the supplied frame, with entries signature.
Redundant ambient tangent slots are supported when all inputs are tangent.
No field equations, harmonic condition, or angular projection is imposed here.
"""
from fractions import Fraction as F
from itertools import product


def indexed_geometry(q, dq, ddq, signature, mutation=None):
    """Return exact P, J, Ricci Hessian, gauge Hessians and reduction terms.

    QC_actual is d²[g(eps) C(eps)]/deps². QC_background is g0 times
    d²C(eps)/deps². These coincide only when q.C1 = 0.
    """
    n = len(signature)
    I = range(n)
    z = q[0, 0] * 0
    half = F(1, 2)
    P = {(a,b,c): half*(dq[b,c,a]+dq[c,b,a]-dq[a,b,c])
         for a,b,c in product(I, repeat=3)}
    dP = {(e,a,b,c): half*(ddq[e,b,c,a]+ddq[e,c,b,a]-ddq[e,a,b,c])
          for e,a,b,c in product(I, repeat=4)}
    J = {(a,b,c): -2*sum((signature[d]*q[a,d]*P[d,b,c] for d in I), z)
         for a,b,c in product(I, repeat=3)}
    required_dJ = {(c,c,a,b) for a,b,c in product(I,repeat=3)}
    required_dJ.update((b,c,a,c) for a,b,c in product(I,repeat=3))
    dJ = {(e,a,b,c): -2*sum((signature[d]*(dq[e,a,d]*P[d,b,c]
                                +q[a,d]*dP[e,d,b,c]) for d in I), z)
          for e,a,b,c in sorted(required_dJ)}
    C1 = {(a,): sum((signature[b]*P[a,b,b] for b in I), z) for a in I}
    QC_actual = {(a,): -2*sum((signature[b]*signature[c]*q[b,c]*P[a,b,c]
                               for b,c in product(I, repeat=2)), z) for a in I}
    QC_background = {(a,): QC_actual[a,]-2*sum(
        (q[a,d]*signature[d]*C1[d,] for d in I), z) for a in I}
    dQC = {(e,a): -2*sum((signature[b]*signature[c]*(
        dq[e,b,c]*P[a,b,c]+q[b,c]*dP[e,a,b,c])
        for b,c in product(I, repeat=2)), z) for e,a in product(I, repeat=2)}
    QRic = {}
    harmonic_reduction = {}
    general_reduction = {}
    for a,b in product(I, repeat=2):
        div = sum((signature[c]*(dJ[c,c,a,b]-dJ[b,c,a,c]) for c in I), z)
        square = sum((signature[c]*signature[d]*(P[c,c,d]*P[d,a,b]
                       -P[c,b,d]*P[d,a,c]) for c,d in product(I, repeat=2)), z)
        QRic[a,b] = div + (1 if mutation == 'ricci_square_half' else 2)*square
        h = half*(dQC[a,b]+dQC[b,a])
        harmonic_reduction[a,b] = h
        general_reduction[a,b] = h-2*sum(
            (signature[c]*P[c,a,b]*C1[c,] for c in I), z)
    return dict(P=P, J=J, C1=C1, QC=QC_actual, QC_actual=QC_actual,
                QC_background=QC_background, dQC=dQC, QRic=QRic,
                harmonic_reduction=harmonic_reduction,
                general_reduction=general_reduction)


def natural_sources(q, eta, deta, ddeta, background_hessian, background_gradient,
                    g0, signature, geometry, Vp, Vpp, Vppp,
                    harmonic=True, mutation=None):
    """Natural Einstein/scalar and reduced Einstein second derivatives.

    harmonic=True applies the registered C1=0 simplification. False retains
    all C1 terms and gives the general affine-path reduced Hessian.
    """
    I = range(len(signature)); z = eta*0; P = geometry['P']; C1 = geometry['C1']
    two = 1 if mutation == 'matter_half' else 2
    matter = {(a,b): -two*deta[a,]*deta[b,]-F(2,3)*(
        Vpp*eta*eta*g0[a,b]+2*Vp*eta*q[a,b]) for a,b in product(I, repeat=2)}
    QE = {k: geometry['QRic'][k]+matter[k] for k in matter}
    reduction = geometry['harmonic_reduction' if harmonic else 'general_reduction']
    QEred = {k: QE[k]-reduction[k] for k in QE}
    qsq = sum((signature[a]*signature[b]*signature[d]*q[a,d]*q[d,b]
               *background_hessian[a,b] for a,b,d in product(I, repeat=3)), z)
    qh = sum((signature[a]*signature[b]*q[a,b]*ddeta[a,b]
              for a,b in product(I, repeat=2)), z)
    qp = sum((signature[a]*signature[b]*signature[c]*q[a,b]*P[c,a,b]
              *background_gradient[c,] for a,b,c in product(I, repeat=3)), z)
    Qphi = 2*(qsq-qh+qp)-Vppp*eta*eta
    if not harmonic:
        Qphi += 2*sum((signature[c]*signature[d]*q[c,d]*C1[d,]
                      *background_gradient[c,] for c,d in product(I, repeat=2)), z)
        Qphi -= 2*sum((signature[c]*C1[c,]*deta[c,] for c in I), z)
    return dict(QE_matter=matter, QE=QE, QEred=QEred, Qphi=Qphi)
