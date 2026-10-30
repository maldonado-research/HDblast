#!/usr/bin/env python3
"""Independent exact fixtures for folded causal shell-source quadrature.

Only Python's standard library is used. Sources are arbitrary physical-side
polynomials on s+y>=1, reflected after a noncommuting matrix acts. Neither
the project's evaluator nor its quadrature implementation is imported.
"""
from fractions import Fraction as Q
from math import comb, factorial
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def clean(p): return {k: Q(v) for k, v in p.items() if v}
def add(p, q):
    z = dict(p)
    for k, v in q.items(): z[k] = z.get(k, Q(0)) + v
    return clean(z)
def scale(p, a): return clean({k: a*v for k, v in p.items()})
def sub(p, q): return add(p, scale(q, -1))
def at(p, t): return sum((v*t**k for k, v in p.items()), Q(0))
def J(p):
    z = {k+1: v/Q(k+1) for k, v in p.items()}
    z[0] = -sum(z.values(), Q(0))
    return clean(z)
def pd(p): return clean({k-1: k*v for k, v in p.items() if k})


def triangle(f, a):
    """Physical monomial integral over s,y>=0, 1<=s+y<=a."""
    return sum((v*Q(factorial(i)*factorial(j), factorial(i+j+2))
                *(a**(i+j+2)-1) for (i,j),v in f.items()), Q(0))


def triangle_poly(f):
    out = {}
    for (i,j), v in f.items():
        c = v*Q(factorial(i)*factorial(j), factorial(i+j+2))
        out[i+j+2] = out.get(i+j+2, Q(0))+c
        out[0] = out.get(0, Q(0))-c
    return clean(out)


def physical_moment(f, a, n):
    """Expand (a-s-y)^n independently in physical coordinates."""
    out = Q(0)
    for k in range(n+1):
        for l in range(n-k+1):
            c = Q((-1)**(k+l)*a**(n-k-l),
                  factorial(k)*factorial(l)*factorial(n-k-l))
            out += c*triangle({(i+k,j+l): v for (i,j),v in f.items()}, a)
    return out


def char_moment(f, a, n, alo, ahi, tlo, thi):
    """Exact polynomial integration of alpha/theta rectangular panels."""
    out = Q(0)
    for (i,j), v in f.items():
        for k in range(i+1):
            for l in range(j+1):
                angular = ((thi**(k+l+1)-tlo**(k+l+1))/Q(k+l+1)
                           *(-1)**k*comb(i,k)*comb(j,l))
                for h in range(n+1):
                    power = i+j+h+2
                    radial = ((ahi**power-alo**power)/power
                              *(-1)**h*comb(n,h)*a**(n-h))
                    out += v*angular*radial/Q(2**(i+j+1)*factorial(n))
    return out


def mass(a, n, lo, hi, theta_width=Q(2)):
    return char_moment({(0,0): Q(1)}, a, n, lo, hi, 0, theta_width)


def edge(f, a):
    return sum((v*Q(factorial(i)*factorial(j),factorial(i+j+1))
                *a**(i+j+1) for (i,j),v in f.items()), Q(0))


def line(f, c, d, lo, hi):
    out = Q(0)
    for (i,j), v in f.items():
        for k in range(j+1):
            power = i+k+1
            out += v*comb(j,k)*c**(j-k)*d**k*(hi**power-lo**power)/power
    return out


def strip(f, lower, upper, lo, hi):
    ant = {(i,j+1): v/Q(j+1) for (i,j),v in f.items()}
    return line(ant,*upper,lo,hi)-line(ant,*lower,lo,hi)


def physical_K(f, p, t, r):
    """Direct physical regions, valid t-r>=1 and r>=0."""
    need(t-r >= 1 and r >= 0, 'physical branch')
    direct = (strip(f,(1,-1),(t+r,-1),Q(0),Q(1))
              +strip(f,(0,0),(t+r,-1),Q(1),t-r)
              +strip(f,(-t+r,1),(t+r,-1),t-r,t))
    reflected = (strip(f,(1,-1),(t-r,-1),Q(0),Q(1))
                 +strip(f,(0,0),(t-r,-1),Q(1),t-r))
    return (direct+p*reflected)/2


def derivative_by_interpolation(fn, t, step, degree):
    """Derivative at t, independent exact Lagrange polynomial samples."""
    result = Q(0)
    for k in range(degree+1):
        numerator = {0: Q(1)}
        denominator = Q(1)
        for l in range(degree+1):
            if l == k: continue
            numerator = add(scale(numerator,-l),
                            {n+1:c for n,c in numerator.items()})
            denominator *= k-l
        coefficient = numerator.get(1,Q(0))/denominator/step
        result += coefficient*fn(t+k*step)
    return result


def matvec(M, v):
    return [sum((a*b for a,b in zip(row,v)),Q(0)) for row in M]
def mp(M, v):
    out = []
    for row in M:
        p = {}
        for a,b in zip(row,v): p=add(p,scale(b,a))
        out.append(p)
    return out
def matmul(A,B): return [matvec(list(zip(*B)),row) for row in A]
def diag(P): return [[Q(P[i]) if i==j else Q(0) for j in range(len(P))] for i in range(len(P))]
def iadd(a,b): return (a[0]+b[0],a[1]+b[1])
def imul(a,b):
    z=[x*y for x in a for y in b]; return min(z),max(z)
def ipow(a,n):
    z=(Q(1),Q(1))
    for _ in range(n): z=imul(z,a)
    return z
def irange(f,alo,ahi,tlo,thi):
    s=imul((alo/2,ahi/2),(1-thi,1-tlo))
    y=imul((alo/2,ahi/2),(1+tlo,1+thi))
    z=(Q(0),Q(0))
    for (i,j),v in f.items(): z=iadd(z,imul((v,v),imul(ipow(s,i),ipow(y,j))))
    return z


def run():
    groups=[]; controls=[]; a=Q(9,8); d=a-1
    P=(1,-1,1,1)
    PN=diag(tuple((p+1)//2 for p in P)); PD=diag(tuple((1-p)//2 for p in P))
    A=((2,3,-1,5),(7,-2,4,-1),(3,2,1,6),(-2,5,3,4))
    M=((0,0,Q(1,7),0),(0,0,0,0),(0,0,0,Q(1,11)),(0,0,0,0))
    Zmat=[[Q(0)]*4 for _ in range(4)]
    need(matmul(M,matmul(M,M))==Zmat,'nilpotent order-three matrix')
    need(matmul(M,M)!=Zmat,'second Robin moment must contribute')
    need(matmul(A,diag(P))!=matmul(diag(P),A),'noncommuting coupling/parity fixture')
    need(matmul(M,PN)==matmul(PN,M)==[list(r) for r in M], 'actual mixed Robin structure')
    groups+=['noncommuting_source_coupling','nonzero_second_Robin_moment','mixed_Robin_structure']
    raw=[{(0,0):Q(1),(2,1):Q(2,3)}, {(1,0):Q(-2),(0,2):Q(3,5)},
         {(1,1):Q(4),(3,0):Q(-1,7)}, {(0,1):Q(5,3),(2,2):Q(2)}]
    source=[]
    for row in A:
        z={}
        for c,f in zip(row,raw):
            for k,v in f.items(): z[k]=z.get(k,Q(0))+c*v
        source.append(clean(z))
    moments=[]
    for n in range(5):
        m=[]
        for f in source:
            exact=physical_moment(f,a,n)
            total=Q(0); enclosed=(Q(0),Q(0)); mass_sum=Q(0)
            for k in range(8):
                lo=1+d*k/8; hi=1+d*(k+1)/8
                for l in range(8):
                    tlo=-1+Q(l,4); thi=tlo+Q(1,4)
                    total+=char_moment(f,a,n,lo,hi,tlo,thi)
                    w=mass(a,n,lo,hi,thi-tlo)
                    mass_sum+=w
                    enclosed=iadd(enclosed,imul((w,w),irange(f,lo,hi,tlo,thi)))
            area=d**(n+1)/factorial(n)*(a/Q(n+1)-d/Q(n+2))
            need(total==exact,'physical and characteristic weighted integral')
            need(mass_sum==area,'partition exact positive masses')
            need(enclosed[0]<=exact<=enclosed[1],'whole-panel range encloses weighted integral')
            m.append(exact)
        moments.append(m)
        groups+=['physical_characteristic_moment_'+str(n),'positive_panel_mass_'+str(n),
                 'continuous_interval_fixture_'+str(n)]
    need(moments[0]!=[2*x for x in moments[0]],'Jacobian omission detectable')
    controls.append('omitted_Jacobian_rejected')
    h=[triangle_poly(f) if p==1 else {} for f,p in zip(source,P)]
    Z=[dict(v) for v in h]
    term=h
    for _ in range(2):
        term=[scale(J(v),-1) for v in mp(M,term)]
        Z=[add(u,v) for u,v in zip(Z,term)]
    need([add(u,v) for u,v in zip(Z,[J(v) for v in mp(M,Z)])]==h,
         'polynomial trace solves Volterra equation')
    folded=matvec(PN,moments[0]); power=matvec(M,moments[1])
    folded=[x-y for x,y in zip(folded,power)]
    folded=[x+y for x,y in zip(folded,matvec(matmul(M,M),moments[2]))]
    zv=[at(z,a) for z in Z]
    need(zv==folded,'folded moments equal independently integrated polynomial Volterra solution')
    groups+=['polynomial_Volterra_equation','folded_source_moments_equal_trace']
    load=mp(M,Z); load_primitive=[J(v) for v in load]
    def corrected(t,r):
        return [physical_K(f,p,t,r)-at(b,t-r) for f,p,b in zip(source,P,load_primitive)]
    need(corrected(a,Q(0))==zv,'physical reflected response plus lift has trace Z')
    e=[edge(f,a) for f in source]; mz=matvec(M,zv)
    dx=[p-u for p,u in zip(matvec(PN,e),mz)]
    dr=[p+u for p,u in zip(matvec(PD,e),mz)]
    deg=max(i+j for f in source for i,j in f)+6
    step=d/Q(2*deg)
    for i in range(4):
        numx=derivative_by_interpolation(lambda t:corrected(t,Q(0))[i],a,step,deg)
        numr=derivative_by_interpolation(lambda r:corrected(a,r)[i],Q(0),step,deg)
        need(numx==dx[i],'time derivative of independent physical response')
        need(numr==dr[i],'normal derivative of independent physical response')
    need(matvec(PD,zv)==[0]*4 and matvec(PD,dx)==[0]*4,'Dirichlet value/time law')
    need(matvec(PN,dr)==mz,'homogeneous actual Robin law')
    groups+=['physical_response_trace','physical_time_derivative','physical_normal_derivative',
             'Dirichlet_shell_value_time','Robin_shell_normal']
    wrong_reflection=matvec(A,matvec(diag(P),[triangle(f,a) for f in raw]))
    correct_reflection=matvec(diag(P),[triangle(f,a) for f in source])
    need(wrong_reflection!=correct_reflection,'source matrix must precede output reflection')
    controls.append('commuting_coupling_through_parity_rejected')
    need(dx!=[p+u for p,u in zip(matvec(PN,e),mz)],'lift sign corruption')
    controls.append('wrong_lift_sign_rejected')
    need(folded!=[x-y for x,y in zip(matvec(PN,moments[0]),matvec(M,moments[1]))],
         'omitted second Robin moment')
    controls.append('omitted_Robin_moment_rejected')
    need(dr!=matvec(M,zv),'Dirichlet edge derivative must be retained')
    controls.append('dropped_Dirichlet_normal_rejected')
    need(dx!=[d/a*p-u for p,u in zip(matvec(PN,e),mz)],'edge length is a not d')
    controls.append('wrong_edge_length_rejected')
    # General factorial tails, independent of the finite fixture.
    m=Q(7,100); N=3
    areas=[d**(j+1)/factorial(j)*(a/Q(j+1)-d/Q(j+2)) for j in range(32)]
    for j in range(31):
        need(areas[j+1]<=d/Q(j+1)*areas[j],'positive weighted-area ratio')
    head=m**(N+1)*areas[N+1]
    tail_majorant=head/(1-m*d/Q(N+2))
    finite_tail=sum((m**j*areas[j] for j in range(N+1,32)),Q(0))
    need(finite_tail<=tail_majorant,'factorial-tail geometric majorant')
    need(tail_majorant>head,'a retained first omitted term is not the full tail')
    groups+=['weighted_area_ratio','factorial_tail_majorant']
    controls.append('omitted_tail_rejected')
    return {'status':'PASS_EXACT_FOLDED_CAUSAL_GEOMETRY','named_groups':groups,
            'group_count':len(groups),'corruption_controls':controls,
            'control_count':len(controls),'source_components':4,'moments_tested':5,
            'panels_per_moment_component':64,'point':['9/8','0'],
            'scope':'Exact manufactured geometry and algebra controls; not an actual HDBLAST source evaluation.'}


if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True,indent=2))
