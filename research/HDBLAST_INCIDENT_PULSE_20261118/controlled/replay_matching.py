"""Independent standard-library fixed-step replay of the two-ended B1 initial data."""
from bisect import bisect_left
import json
from math import exp, log, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = json.loads((ROOT/'MATCHING_INPUT.json').read_text())


def check(condition, message):
    if not condition:
        raise ValueError(message)


def profile(z, center, width):
    s = (z-center)/width
    if abs(s) >= 1:
        return 0., 0.
    d = 1-s*s
    h = exp(1-1/d)
    return h, -2*s*h/(d*d*width)


def potential(phi):
    w = 1-phi+phi**3/3
    dw = phi*phi-1
    return dw*dw/2-2*w*w/3, 2*phi*dw-4*w*dw/3


def wall():
    a = SOURCE['params']
    phi = SOURCE['phi_b']
    W = 1-phi+phi**3/3
    s = 2*W+a['delta']*(1+a['c']*phi+a['d']*phi**2/2)
    sp = 2*(phi**2-1)+a['delta']*(a['c']+a['d']*phi)
    U, Up = potential(phi)
    # Assign the exactly compatible radius algebraically before floating evaluation.
    rho = 1/sqrt(s*s/36-sp*sp/48+U/6)
    q = rho*s/6
    r = -rho*sp/2
    b = log(rho)
    determinant = q*r*r-2*rho*rho*r*Up/3
    return (b,q,phi,r,b,q,1.), dict(compatible_radius=rho, archived_radius=SOURCE['rho_b'],
                                     bracket_determinant=determinant)


def sum_step(state, change, step):
    return [a+step*b for a,b in zip(state,change,strict=True)]


def stage(z, state, eta, coefficients, mode):
    b,q,phi,r,B,Q,u = state[:7]
    f, df = profile(z,SOURCE['scalar_center'],SOURCE['scalar_width'])
    h = [profile(z,center,mode['width'])[0] for center in mode['centers']]
    v = sum(c*x for c,x in zip(coefficients,h,strict=True))
    base_U, dU = potential(phi)
    base_E = exp(2*b)
    # Use the reference's static evolution equation, independent of the
    # adaptive generator's Hamiltonian equation for its reference solution.
    bzz = 3-3*q*q-2*base_E*base_U/3
    rz = base_E*dU-3*q*r
    p = -eta*df
    radial = r+eta*df
    trial_U = potential(phi+eta*f)[0]
    EU = exp(2*B)*trial_U
    Qz = u*u+u*v-Q*Q-(p*p+radial*radial)/6-EU/3
    uz = Q*v-p*radial/3
    out = [q,bzz,r,rz,Q,Qz,uz]
    for j, basis in enumerate(h):
        i = 7+3*j
        dB,dQ,du = state[i:i+3]
        # Tangent equation about this particular controlled state.
        out.extend((dQ, -2*EU*dB/3-2*Q*dQ+(2*u+v)*du+u*basis, v*dQ+Q*basis))
    return out


def determinant(m):
    a,b,c = m
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def condition_number(m):
    # Inverse-row sums from transposed cofactors.
    d = abs(determinant(m))
    inverse_row_sums = []
    for i in range(3):
        total = 0.
        for j in range(3):
            sub = [[m[r][c] for c in range(3) if c!=i] for r in range(3) if r!=j]
            total += abs(sub[0][0]*sub[1][1]-sub[1][0]*sub[0][1])/d
        inverse_row_sums.append(total)
    return max(sum(abs(x) for x in row) for row in m)*max(inverse_row_sums)


def diagnostics(rows, eta, coef, mode, step):
    # Samples were integrated from the wall towards negative Z. Use increasing
    # Z for the independent centered-difference residuals and field interpolation.
    fields = list(reversed(rows))
    dz = -step
    cols = list(zip(*fields,strict=True))
    B = cols[5]
    u = cols[7]
    phi = []
    p = []
    v = []
    for z,_,_,phi0,*_ in fields:
        f, df = profile(z,SOURCE['scalar_center'],SOURCE['scalar_width'])
        phi.append(phi0+eta*f)
        p.append(-eta*df)
        v.append(sum(c*profile(z,x,mode['width'])[0] for c,x in zip(coef,mode['centers'],strict=True)))

    def first(y,i):
        return (y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*dz)

    def second(y,i):
        return (-y[i+2]+16*y[i+1]-30*y[i]+16*y[i-1]-y[i-2])/(12*dz*dz)

    maxH = 0.
    maxM = 0.
    for i in range(2,len(fields)-2):
        Q = first(B,i)
        r = first(phi,i)
        H_terms = [-2*exp(2*B[i])*potential(phi[i])[0],6*u[i]**2,6*u[i]*v[i],
                   -12*Q*Q,6*Q*Q,-6*second(B,i),-p[i]*p[i],-r*r]
        M_terms = [-3*first(u,i),-3*u[i]*Q,3*u[i]*Q,3*Q*v[i],-p[i]*r]
        maxH = max(maxH,abs(sum(H_terms))/sum(map(abs,H_terms)))
        maxM = max(maxM,abs(sum(M_terms))/sum(map(abs,M_terms)))

    # Fixed scalar level, comparing coordinate scalars. This probe lies inside
    # the original pulse and outside the delayed control support.
    index = min(range(len(fields)), key=lambda i:abs(fields[i][0]+0.11))
    z,b,q,basePhi,r,Btrial,Qtrial,utrial = fields[index]
    _,df = profile(z,SOURCE['scalar_center'],SOURCE['scalar_width'])
    X = exp(-2*Btrial)*((r+eta*df)**2-(eta*df)**2)
    level = phi[index]
    quiet_phi = cols[3]
    j = bisect_left(quiet_phi,level)
    check(0<j<len(fields) and all(a<b for a,b in zip(quiet_phi,quiet_phi[1:])),
          'Reference scalar has to be monotone and bracket the comparison level')
    weight = (level-quiet_phi[j-1])/(quiet_phi[j]-quiet_phi[j-1])
    X0 = [(exp(-2*cols[1][k])*cols[4][k]**2) for k in (j-1,j)]
    matched_X0 = X0[0]+weight*(X0[1]-X0[0])
    return dict(max_normalized_H=maxH,max_normalized_M=maxM,z_probe=z,
                comparison_phi=level,X=X,quiet_X_at_equal_phi=matched_X0,
                invariant_difference=X-matched_X0)


def run(eta, coef, mode, intervals, include_checks=True):
    initial,_ = wall()
    y = [*initial,*([0.]*9)]
    h = mode['cut']/intervals
    rows = [(0.,*y[:7])]
    for i in range(intervals):
        z = i*h
        a = stage(z,y,eta,coef,mode)
        b = stage(z+h/2,sum_step(y,a,h/2),eta,coef,mode)
        c = stage(z+h/2,sum_step(y,b,h/2),eta,coef,mode)
        d = stage(z+h,sum_step(y,c,h),eta,coef,mode)
        y = [yy+(h/6)*(aa+2*bb+2*cc+dd)
             for yy,aa,bb,cc,dd in zip(y,a,b,c,d,strict=True)]
        rows.append(((i+1)*h,*y[:7]))
    b,q,_,_,B,Q,u = y[:7]
    jacobian = [[y[7+i+3*j] for j in range(3)] for i in range(3)]
    out = dict(intervals=intervals,raw_outer_difference=[B-b,Q-q,u-1],
               jacobian=jacobian,determinant=determinant(jacobian),
               condition_infinity=condition_number(jacobian))
    if include_checks:
        out.update(diagnostics(rows,eta,coef,mode,h))
    return out


def main():
    _, metadata = wall()
    results = {}
    for name, mode in SOURCE['profiles'].items():
        quiet = run(0.,[0.]*3,mode,1600,False)
        trials = {}
        signed_runs = [(1,'coefficients'),(-1,'negative_coefficients')]
        if name == 'delayed':
            signed_runs.extend([(0.5,'half_positive'),(-0.5,'half_negative')])
        for sign, key in signed_runs:
            eta = sign*SOURCE['amplitude']
            coefficients = mode[key]
            coarse = run(eta,coefficients,mode,800)
            fine = run(eta,coefficients,mode,1600)
            omitted = run(eta,[0.]*3,mode,800,False)
            check(max(map(abs,fine['raw_outer_difference']))<3e-9,'Matching residual')
            check(max(map(abs,omitted['raw_outer_difference']))>1e-5,'Omitted controls')
            check(max(fine['max_normalized_H'],fine['max_normalized_M'])<2e-4,'Constraint residual')
            check(abs(fine['invariant_difference'])>1e-4,'Invariant must differ at fixed scalar')
            trials[str(eta)] = dict(coarse=coarse,fine=fine,
                                    control_omission_endpoint=omitted['raw_outer_difference'])
        results[name] = dict(quiet=quiet,trials=trials)
    print(json.dumps(dict(background=metadata,profiles=results),indent=2,allow_nan=False))


if __name__=='__main__':
    main()
