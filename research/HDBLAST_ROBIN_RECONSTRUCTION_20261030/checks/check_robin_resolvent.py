#!/usr/bin/env python3
"""Independent exact polynomial controls for Robin Volterra reconstruction.

No numerical HDBLAST evaluator or upstream scientific code is imported.
The fixture uses mixed projectors, dense noncommuting bulk couplings,
distinct frozen/actual Robin matrices, and nonzero affine shell forcing.
"""
from fractions import Fraction as Q
from math import factorial
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def poly(values):
    out = list(map(Q, values))
    while out and out[-1] == 0:
        out.pop()
    return out


def add(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, z in enumerate(a): out[i] += z
    for i, z in enumerate(b): out[i] += z
    return poly(out)


def scale(a, c): return poly(c*z for z in a)
def sub(a, b): return add(a, scale(b, -1))
def deriv(a): return poly(i*a[i] for i in range(1, len(a)))
def integ(a): return poly([Q(0)] + [z/Q(i+1) for i, z in enumerate(a)])
def times_x(a): return poly([Q(0)] + a)
def reflect(a): return poly((-1)**i*z for i, z in enumerate(a))
def evaluate(a, x): return sum((z*x**i for i, z in enumerate(a)), Q(0))


def primitive(a, k):
    for _ in range(k): a = integ(a)
    return a


def va(a, b): return [add(x, y) for x, y in zip(a, b)]
def vs(a, b): return [sub(x, y) for x, y in zip(a, b)]
def vc(a, c): return [scale(x, c) for x in a]
def J(a, k=1): return [primitive(x, k) for x in a]
def dx(a): return [deriv(x) for x in a]
def vv(a, x): return [evaluate(p, x) for p in a]


def mv(matrix, vector):
    out = []
    for row in matrix:
        p = []
        for c, z in zip(row, vector): p = add(p, scale(z, c))
        out.append(p)
    return out


def matrix(rows): return [list(map(Q, row)) for row in rows]
def ma(a, b): return [[u+v for u, v in zip(ar, br)] for ar, br in zip(a, b)]
def ms(a, b): return [[u-v for u, v in zip(ar, br)] for ar, br in zip(a, b)]
def PN(a): return [p if parity > 0 else [] for p, parity in zip(a, P)]
def PD(a): return [p if parity < 0 else [] for p, parity in zip(a, P)]
def parity(a): return [scale(p, sign) for p, sign in zip(a, P)]


P = (1, -1, 1)
A = matrix(((1, 2, -1), (3, -2, 4), (2, 1, 3)))
B = matrix(((-2, 3, 1), (1, 4, -3), (5, -1, 2)))
D = matrix(((2, -1, 4), (-3, 2, 1), (1, 5, -2)))
M = matrix(((Q(1,10), 0, Q(1,20)), (0, 0, 0), (Q(-1,30), 0, Q(1,12))))
M0 = matrix(((Q(1,11), 0, Q(1,21)), (0, 0, 0), (Q(-1,31), 0, Q(1,13))))
need(any(A[i][j] and P[i] != P[j] for i in range(3) for j in range(3)),
     'noncommuting parity fixture')
need(any(sum(M[i][k]*A[k][j] for k in range(3)) !=
         sum(A[i][k]*M[k][j] for k in range(3))
         for i in range(3) for j in range(3)), 'Robin/bulk noncommuting fixture')
checks = ['noncommuting_parity_and_bulk_fixture']

# Check the endpoint formulas for every temporal moment of the physical
# reflected source trace; direct univariate integration is the reference.
u = [poly([0,0,0,0,1]), poly([0,0,0,0,0,2]), poly([0,0,0,0,0,0,-1])]
v = [poly([0,0,0,0,0,0,-2]), poly([0,0,0,0,0,0,0,1]), poly([0,0,0,0,3])]
f0 = vc(va(mv(ma(A,B), va(dx(u),v)), mv(D, va(u,J(v)))), Q(1,2))
g0 = vc(va(mv(ms(B,A), vs(dx(u),v)), mv(D, vs(u,J(v)))), Q(1,2))
H = va(vs([times_x(p) for p in J(f0)], J(f0,2)),
       vc(va(J(g0,2), [reflect(p) for p in J(g0,2)]), Q(1,2)))
h0 = PN(va(va(u,J(v)),H))
for j in range(6):
    closedH = va(vs([times_x(p) for p in J(f0,j+1)], vc(J(f0,j+2),j+1)),
                 vc(va(J(g0,j+2),vc([reflect(p) for p in J(g0,j+2)],(-1)**j)),Q(1,2)))
    need(closedH == J(H,j), 'source trace primitive endpoint formula')
    need(PN(va(va(J(u,j),J(v,j+1)),closedH)) == J(h0,j),
         'complete h primitive endpoint formula')
    checks.append('endpoint_moment_' + str(j))
    k=j+1
    direct_minus = [reflect(p) for p in J(g0,k)]
    closed_minus = vc(vs(mv(ms(B,A), parity(va(J(u,k-1),J(v,k)))),
                         mv(D,parity(va(J(u,k),J(v,k+1))))),Q((-1)**(k-1),2))
    need(direct_minus == closed_minus, 'negative ridge primitive parity after coupling')
    checks.append('negative_primitive_' + str(k))

# Manufacture an exact polynomial trace Z, nonzero G, and h from the actual
# Volterra equation. This gives a separate exact solution for testing the
# resolvent, rather than comparing two implementations of its same series.
Z = [poly([0,0,0,0,Q(1,3),Q(-1,2)]), [], poly([0,0,0,Q(1,5),0,0,Q(1,7)])]
G = [poly([0,0,0,Q(1,7),0,Q(-1,11)]), [], poly([0,0,0,0,Q(-2,9)])]
q = va(mv(M,Z),G)
h = va(Z,J(q))
need(vs(h,J(q)) == Z, 'exact trace integral equation')
need(vs(dx(h),q) == dx(Z), 'exact trace derivative equation')
need(PD(Z) == [[],[],[]] and PD(G) == [[],[],[]], 'mixed boundary support')
checks.extend(('exact_Volterra_solution','exact_trace_derivative','mixed_boundary_support'))

# Coefficient dictionary in (x,r). Construct an independently chosen base S
# with trace h and zero PN normal derivative, then add the characteristic lift.
def bi_add(a,b):
    out=dict(a)
    for k,z in b.items():out[k]=out.get(k,Q(0))+z
    return {k:z for k,z in out.items() if z}
def bi_scale(a,c):return {k:c*z for k,z in a.items() if c*z}
def bi_deriv(a,axis):
    out={}
    for ij,z in a.items():
        if ij[axis]:
            kk=list(ij);kk[axis]-=1;out[tuple(kk)]=z*ij[axis]
    return out
def bi_shift(a,k):return {(i,j+k):z for (i,j),z in a.items()}
def xpoly(p):return {(i,0):z for i,z in enumerate(p) if z}
def characteristic(p):
    from math import comb
    out={}
    for n,z in enumerate(p):
        for j in range(n+1):out[(n-j,j)]=out.get((n-j,j),Q(0))+z*comb(n,j)*(-1)**j
    return {k:z for k,z in out.items() if z}
def shell(p):
    degree=max((i for (i,j) in p if j==0),default=-1)
    out=[Q(0)]*(degree+1)
    for (i,j),z in p.items():
        if j==0:out[i]+=z
    return poly(out)
base=[bi_add(xpoly(h[0]),bi_shift(xpoly(poly([0,0,1])),2)),
      bi_add(bi_shift(xpoly(poly([0,0,0,2])),1),bi_shift(xpoly(poly([0,0,3])),3)),
      bi_add(xpoly(h[2]),bi_shift(xpoly(poly([0,0,-1])),2))]
lift=[bi_scale(characteristic(p),-1) for p in J(q)]
W=[bi_add(s,l) for s,l in zip(base,lift)]
Wt=[shell(p) for p in W]
Wr=[shell(bi_deriv(p,1)) for p in W]
need(PN(Wt)==Z and PD(Wt)==[[],[],[]], 'reconstruction trace and Dirichlet law')
need(vs(vs(PN(Wr),mv(M,PN(Wt))),G)==[[],[],[]], 'exact actual affine Robin law')
checks.extend(('reconstruction_trace','actual_affine_Robin_law'))
for i in range(3):
    waveW=bi_add(bi_deriv(bi_deriv(W[i],0),0),bi_scale(bi_deriv(bi_deriv(W[i],1),1),-1))
    waveS=bi_add(bi_deriv(bi_deriv(base[i],0),0),bi_scale(bi_deriv(bi_deriv(base[i],1),1),-1))
    need(waveW==waveS,'characteristic lift leaves principal bulk source unchanged')
    checks.append('principal_bulk_identity_' + str(i))
need(any(vs(vs(PN(Wr),mv(M0,PN(Wt))),G)), 'substituting frozen M0 must fail')
need(any(vs(PN(Wr),mv(M,PN(Wt)))), 'omitting affine G must fail')
wrongW=[bi_add(s,bi_scale(l,-1)) for s,l in zip(base,lift)]
wrongtrace=[shell(p) for p in wrongW];wrongr=[shell(bi_deriv(p,1)) for p in wrongW]
need(any(vs(vs(PN(wrongr),mv(M,PN(wrongtrace))),G)), 'wrong lift sign must fail')
checks.extend(('wrong_frozen_M_rejected','omitted_G_rejected','wrong_lift_sign_rejected'))

# Exact finite-series residual and independent factorial tail inequalities.
def power_apply(m,z,n):
    for _ in range(n):z=mv(m,z)
    return z
def cap(polys,d):return max(sum(abs(c)*d**j for j,c in enumerate(p)) for p in polys)
def vmax(values):return max(map(abs,values))
duration=Q(1,4)
m=max(sum(map(abs,row)) for row in M)
a=m*duration
Hcap=cap(h,duration);Gcap=cap(G,duration)
need(a<1,'fixture stable rational tail denominator')
source_bound=duration*Gcap/(1-a/2)
for N in range(6):
    zh=[[],[],[]];zg=[[],[],[]]
    for j in range(N+1):
        zh=va(zh,vc(power_apply(M,J(h,j),j),(-1)**j))
        zg=va(zg,vc(power_apply(M,J(G,j+1),j),(-1)**(j+1)))
    zn=va(zh,zg)
    rhs=vs(h,J(G))
    residual=vs(va(zn,mv(M,J(zn))),rhs)
    expected=vc(power_apply(M,J(rhs,N+1),N+1),(-1)**N)
    need(residual==expected,'finite Neumann residual telescopes')
    need(any(residual),'finite series cannot be labeled exact Robin trace')
    tailh=Hcap*a**(N+1)/Q(factorial(N+1))/(1-a/Q(N+2))
    tailg=duration*Gcap*a**(N+1)/Q(factorial(N+2))/(1-a/Q(N+3))
    for k in range(17):
        t=duration*k/16
        need(vmax(vv(vs(Z,zn),t))<=tailh+tailg,'factorial full resolvent tail')
        need(vmax(vv(vs(Z,zh),t))<=tailh+source_bound,'unresolved affine source bound')
    checks.extend(('finite_series_residual_'+str(N),'factorial_tail_'+str(N),'affine_source_bound_'+str(N)))

# Correction forced by an old candidate boundary residual: d=Z-V,
# d=J R-M Jd; lift forcing Q=M d-R. This checks the transfer sign.
Z0=PN(u)
V=vs(h,mv(M0,J(Z0)))
R=vs(vs(mv(M0,Z0),mv(M,V)),G)
correction=vs(Z,V)
need(correction==vs(J(R),mv(M,J(correction))),'boundary residual correction integral equation')
Qforcing=vs(mv(M,correction),R)
need(vs(correction,vc(J(Qforcing),-1))==[[],[],[]],'correction characteristic-lift sign')
checks.extend(('boundary_residual_correction','correction_lift_sign'))

print(json.dumps({'status':'PASS','named_check_groups':len(checks),'checks':checks,
 'method':'exact engineered polynomial Volterra solution and bivariate characteristic lift; independent direct primitives for trace moments; rational factorial tail controls',
 'scope':'finite algebra controls for analytic exact Robin reconstruction and its enclosure; not full coupled HDBLAST evolution, numerical C2 error, external specialist review or physical observable'},sort_keys=True))
