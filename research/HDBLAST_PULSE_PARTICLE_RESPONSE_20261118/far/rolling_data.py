"""Constraint-matched Cauchy perturbation on a rolling B1 background."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


def add_pulse(run, T, state, m, amplitude, options):
    import static_w as S
    fields = run.fields(T, state, m)
    z = run.z
    center, width = -.12, .02
    control, spread = options.get('control', -.85), options.get('spread', .25)
    left, right = control - spread, center + width
    active = (z >= left) & (z <= right)
    local = (z >= left - .02) & (z <= right + .02)
    names = ('A_Z', 'A_T', 'B', 'B_T', 'B_Z', 'phi', 'phi_Z', 'phi_T')
    interp = CubicSpline(z[local], np.stack([fields[k][local] for k in names], axis=1))
    def bump(x, loc, wid):
        t = np.asarray((x-loc)/wid)
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            s = np.maximum(1.-t*t, 1e-18)
            val = np.exp(1.-1./s)
            grad = -2.*t*val/(wid*s*s)
        return val, grad
    def coeff(x, b):
        Q, t0, B0, v0, Bz, phi, r, p = interp(x).T
        shape, grad = bump(x, center, width)
        h, hp = bump(x, control, spread)
        kind=options.get('kind','incident')
        if kind in ('velocity','metric','metric_field'):
            F, Fp = amplitude*shape, amplitude*grad
            P = -Fp+(b*h if kind=='velocity' else 0.)
        else:
            F, Fp = amplitude*shape+b*h, amplitude*grad+b*hp
            P = -Fp
        w, wp = (b*h,b*hp) if kind=='metric' else (0.*h,0.*hp)
        t = t0+w
        dotQ = np.zeros_like(Q)
        ddQ = np.zeros_like(Q)
        if kind=='metric_field':
            u=(x-control)/spread
            inner=np.maximum(1-u*u,1e-18)
            dotQ=b*hp
            ddQ=b*h*(4*u*u/inner**4-2/inner**2-8*u*u/inner**3)/spread**2
        newQ=Q+dotQ
        if np.any(newQ*Q[0]<=0) or np.min(abs(newQ))<1e-7:
            raise ValueError('Degenerate new radial scale gradient')
        D = newQ-t*t/newQ
        K = -(p*Fp+(r+Fp)*P)
        L = P*(2*p+P)+Fp*(2*r+Fp)
        s = (-w*(2*t0+w+v0)
             - t/newQ*(wp-w*(Bz-Q)+(t-v0)*dotQ)
             + 2*dotQ*(2*Q+dotQ)-dotQ*Bz+ddQ
             + L/6+t/newQ*K/3)
        W = 1-phi+phi**3/3
        dW = F*((phi*phi-1)+F*(phi+F/3))
        dWp = F*(2*phi+F)
        dU = .5*dWp*(2*(phi*phi-1)+dWp) - (2/3)*dW*(2*W+dW)
        E = np.exp(2*B0)
        return 2*(E*S.U(phi)/3-s)/D, -2*(s+E*dU/3)/D, t, Q, K, w, wp, Bz, F, P, D, dotQ
    # Solve y'=ay+b, y=exp(-2 delta B)-1, by integrating the two
    # integrating factors with grid quadrature. Check full-state constraints separately.
    from scipy.integrate import cumulative_simpson
    count = int(np.ceil((right-left)/.000025))+1
    x = np.linspace(left, right, count)
    def solve(b, complete=False):
        aa, ab, t, Q, K, w, wp, Bz, F, P, D, dotQ = coeff(x, b)
        if np.any(D*D[0]<=0) or np.min(abs(D))<1.e-6:
            raise ValueError('Characteristic degeneracy in compensation interval')
        I = cumulative_simpson(aa,x=x,initial=0.)
        J = cumulative_simpson(ab*np.exp(-I),x=x,initial=0.)
        d = np.exp(I)*(J-J[-1])
        return (d, aa*d+ab) if complete else d[0]
    fzero = solve(0.)
    if options.get('probe'):
        from scipy.optimize import minimize_scalar
        # This is a numerical profile check on a declared finite interval.
        # Local bounded refinements follow every sampled local minimum.
        bs=np.linspace(-.15,.15,151)
        vals=[]
        for b in bs:
            try:
                vals.append(float(solve(b)))
            except ValueError:
                vals.append(float('nan'))
        minima=[]
        for j in range(1,len(vals)-1):
            if vals[j]<=vals[j-1] and vals[j]<=vals[j+1]:
                opt=minimize_scalar(solve,bounds=(bs[j-1],bs[j+1]),method='bounded',
                                    options={'xatol':1e-11})
                minima.append((float(opt.x),float(opt.fun)))
        finite=[(float(b),f) for b,f in zip(bs,vals) if np.isfinite(f)]
        bmin,fmin=min(minima+finite,key=lambda bf:bf[1])
        return state, dict(amplitude=amplitude,kind=options.get('kind','incident'),
                           lower=-.15,upper=.15,total_samples=151,admissible_samples=len(finite),
                           zero_value=fzero,minimum_at=bmin,minimum=fmin)

    old = {1.: (0., fzero), -1.: (0.,fzero)}
    root = None
    for radius in np.geomspace(1e-8, options.get('radius', .30), 32):
        for sign in old:
            b0, residual0 = old[sign]
            b1 = sign*radius
            try:
                residual1 = solve(b1)
            except ValueError:
                continue
            if residual0*residual1 <= 0:
                root = brentq(solve,min(b0,b1),max(b0,b1),xtol=1e-14)
                break
            old[sign] = b1, residual1
        if root is not None:
            break
    if root is None:
        raise ValueError(f'No compact compensator, uncorrected endpoint {fzero}: {old}')
    d, dp = solve(root, True)
    if d.min()<=-1: raise ValueError('Singular correction to metric')
    y = CubicSpline(x,d)(z[active])
    a, ab, t, Q, K, w, wp, Bz, F, P, D, dotQ = coeff(z[active],root)
    g = -.5*np.log1p(y)
    gp = -.5*(a*y+ab)/(1+y)
    dv = (wp-w*(Bz-Q)+(t-interp(z[active])[:,3])*dotQ-t*gp-K/3)/(Q+dotQ)
    if options.get('kind')=='metric_field':
        changed_A=root*bump(z[active],control,spread)[0]
    else:
        changed_A=0.
    changed = state.copy()
    for k, value in ((0,changed_A),(1,g),(2,F),(3,w),(4,dv),(5,P)):
        changed[k,active] += value
    def residual(fields, state):
        dm = run.matter_rhs(fields,m,T)
        phitt = run.accel(fields)[2][-1]
        return run.constraints(fields,state[3],m,dm,phitt,T)
    before = residual(fields,state)
    after = residual(run.fields(T,changed,m), changed)
    info = dict(T=float(T), mode=run.mode, amplitude=amplitude,center=center,halfwidth=width,
                control_center=control, control_halfwidth=spread, control_kind=options.get("kind","incident"), compensator_coefficient=float(root),
                uncompensated_delta_y_left=float(fzero), matched_delta_y_left=float(d[0]),
                min_abs_D=float(np.min(abs(D))), min_abs_Q=float(np.min(abs(Q+dotQ))),
                max_abs_delta_B=float(max(abs(g))), max_abs_delta_BT=float(max(abs(dv))))
    for j,k in enumerate(('H','M')):
        delta = after[j]-before[j]
        scale = np.maximum(before[j+2],after[j+2])
        info[k+'_relative_delta_max']=float(np.max(abs(delta[local])/scale[local]))
        info[k+'_relative_background_max']=float(np.max(abs(before[j][local]/before[j+2][local])))
    return changed, info
