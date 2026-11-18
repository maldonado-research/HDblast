"""Evolve the initial B1 data in an explicit finite difference domain of dependence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
SOURCE_HASH = 'bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f'
LEFT = -.46
MONITOR_LEFT = -.05
END = .18
CFL = .5
WALL_CROSSING = .5


def potential(phi):
    w = 1 - phi + phi**3 / 3
    wp = phi**2 - 1
    return .5*wp**2 - 2*w**2/3, 2*phi*wp - 4*w*wp/3


def tension(phi, params):
    delta, c, d = (params[key] for key in ('delta', 'c', 'd'))
    w = 1 - phi + phi**3/3
    return (2*w + delta*(1+c*phi+d*phi**2/2),
            2*(phi**2-1) + delta*(c+d*phi))


def bump_sum(z, bumps):
    z = np.asarray(z)
    field = np.zeros_like(z, dtype=float)
    derivative = np.zeros_like(z, dtype=float)
    for center, width, amplitude in bumps:
        s = (z-center)/width
        mask = abs(s) < 1
        d = 1 - s[mask]**2
        f = np.exp(1-1/d)
        field[mask] += amplitude*f
        derivative[mask] -= 2*amplitude*s[mask]*f/(width*d**2)
    return field, derivative


def prepare(summary, cases, grid):
    rho, phi_wall = summary['rho_b'], summary['phi_b_static']
    sigma, sigma1 = tension(phi_wall, summary['params'])
    at_wall = np.array([np.log(rho), rho*sigma/6, phi_wall, -rho*sigma1/2])

    def static(z, y):
        b, q, phi, r = y
        U, Up = potential(phi)
        return q, 1-q*q-r*r/6-np.exp(2*b)*U/3, r, np.exp(2*b)*Up-3*q*r

    def integrate(eq, right, endpoint):
        sol = solve_ivp(eq, (right, LEFT), endpoint, method='DOP853',
                        rtol=2e-13, atol=2e-15, max_step=.0005, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.sol

    static_solution = integrate(static, 0., at_wall)
    b0, q0, phi0, r0 = static_solution(grid)
    C = 1-q0**2+r0**2/12-np.exp(2*b0)*potential(phi0)[0]/6
    preparations = {}
    for name, profile in cases.items():
        scalar, metric = profile['scalar'], profile['metric_velocity']
        F, Fprime = bump_sum(grid, scalar)
        vB, _ = bump_sum(grid, metric)
        B, Q, u = b0.copy(), q0.copy(), np.ones_like(b0)
        if scalar or metric:
            right = max(center+width for center,width,amp in scalar+metric)
            def constraints(z, state):
                b, q, a_t = state
                phi, r = static_solution(z)[2:]
                f, df = bump_sum(z, scalar)
                b_t = bump_sum(z, metric)[0]
                p = -df
                r = r+df
                return (q, a_t*a_t+a_t*b_t-q*q-(p*p+r*r)/6
                        - np.exp(2*b)*potential(phi+f)[0]/3, q*b_t-p*r/3)
            endpoint = [*static_solution(right)[:2], 1.]
            evaluator = integrate(constraints, right, endpoint)
            mask = grid < right
            B[mask], Q[mask], u[mask] = evaluator(grid[mask])
        cut = profile['matching_cut']
        nearest = np.argmin(abs(grid-cut))
        preparations[name] = (np.array([B, B, phi0+F]),
                              np.array([u, vB, -Fprime]),
                              np.array([B[nearest]-b0[nearest], Q[nearest]-q0[nearest], u[nearest]-1]))
    return preparations, float(np.max(abs(C)))


def spatial(fields, params, h):
    # Result on indices 1..last. The farthest-left node is discarded at each step.
    dz = np.empty_like(fields[:, 1:])
    dzz = np.empty_like(dz)
    dz[:, :-1] = (fields[:, 2:]-fields[:, :-2])/(2*h)
    dzz[:, :-1] = ((fields[:, 2:]-fields[:, 1:-1])-(fields[:, 1:-1]-fields[:, :-2]))/h**2
    sigma, sigma1 = tension(fields[2,-1], params)
    normal = np.exp(fields[1,-1])*np.array([sigma/6, sigma/6, -.5*sigma1])
    dz[:, -1] = normal
    # Cubic Hermite closure using two interior nodes and the actual wall Neumann slope.
    dzz[:, -1] = (4*(fields[:, -2]-fields[:, -1]) - .5*(fields[:, -3]-fields[:, -1])
                   + 3*h*normal)/h**2
    return dz, dzz


def source(dz, dzz, geometry, phi, velocity):
    az, bz, rz = dz
    at, bt, pt = velocity
    U, Up = potential(phi)
    e2b = np.exp(2*geometry)
    kinetic = at*at-az*az
    return dzz + np.array([-3*kinetic + 2*e2b*U/3,
                            3*kinetic - .5*pt*pt + .5*rz*rz - e2b*U/3,
                            -3*at*pt + 3*az*rz - e2b*Up])


def constraints(fields, velocity, grid, h):
    mask = (grid[1:-1] >= MONITOR_LEFT-1e-12)
    az, bz, rz = (fields[:,2:]-fields[:,:-2])/(2*h)
    azz = ((fields[0,2:]-fields[0,1:-1])-(fields[0,1:-1]-fields[0,:-2]))/h**2
    atz = (velocity[0,2:]-velocity[0,:-2])/(2*h)
    at,bt,pt = velocity[:,1:-1]
    potential_term = np.exp(2*fields[1,1:-1])*potential(fields[2,1:-1])[0]
    hterms = np.array([-2*potential_term, 6*at**2, 6*at*bt,-12*az**2,6*az*bz,-6*azz,-pt**2,-rz**2])
    mterms = np.array([-3*atz,-3*at*az,3*at*bz,3*az*bt,-pt*rz])
    return [np.max((abs(terms.sum(axis=0))/(abs(terms).sum(axis=0)+1e-300))[mask])
            for terms in (hterms,mterms)]


HEADER = ['T', 'A_wall','B_wall','phi_wall','A_T_wall','B_T_wall','phi_T_wall',
          *[f'phi_minus_{i}' for i in range(1,7)], 'H_constraint_rel','M_constraint_rel']


def record(t, fields, velocity, grid, h):
    return [t, *fields[:,-1], *velocity[:,-1], *fields[2,-2:-8:-1],
            *constraints(fields,velocity,grid,h)]


def evolve(fields, vel, grid, params, h):
    dt = CFL*h
    steps = round(END/dt)
    if not np.isclose(steps*dt, END) or grid[steps+4] > MONITOR_LEFT:
        raise ValueError('Final slice does not contain sufficient diagnostic domain')
    # Every older state is sliced to the same still-valid coordinates before the next update.
    rows = [record(0.,fields,vel,grid,h)]
    dz,dzz = spatial(fields,params,h)
    a0 = source(dz,dzz,fields[1,1:],fields[2,1:],vel[:,1:])
    prev = fields[:,1:]
    now = prev+dt*vel[:,1:]+.5*dt**2*a0
    grid = grid[1:]
    for n in range(1,steps+1):
        dz, dzz = spatial(now,params,h)
        recent = now[:,1:]
        past = prev[:,1:]
        speed_half = (recent-past)/dt
        speed = speed_half
        # Center the velocity-dependent nonlinear source by local fixed-point corrections.
        # Only same-node velocities are updated, preserving the numerical characteristic cone.
        for correction in range(3):
            acc = source(dz,dzz,recent[1],recent[2],speed)
            speed = speed_half + .5*dt*acc
        future = recent + dt*speed_half + dt**2*acc
        velocity = (future-past)/(2*dt)
        grid = grid[1:]
        if not np.isfinite(future).all():
            raise ArithmeticError('Nonfinite evolution')
        if max(recent[2,-1],future[2,-1]) >= WALL_CROSSING:
            raise ArithmeticError('The matter-free system leaves the B1 preproduction branch')
        rows.append(record(n*dt,recent,velocity,grid,h))
        prev, now = recent, future
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('profiles', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--spacing', type=float, required=True)
    args = parser.parse_args()
    raw = args.summary.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_HASH:
        raise ValueError('Incorrect registered B1 summary')
    summary = json.loads(raw)
    cases = json.loads(args.profiles.read_text())
    h = args.spacing
    grid = np.arange(-round(-LEFT/h),1)*h
    if h <= 0 or not np.isclose(grid[0],LEFT):
        raise ValueError('Expected positive spacing dividing patch length')
    args.output.mkdir(exist_ok=True, parents=True)
    initial, residual = prepare(summary,cases,grid)
    metadata = dict(spacing=h,dt=CFL*h,T_end=END,left=LEFT,static_residual_initial_max=residual,
                    constraint_monitor_left=MONITOR_LEFT,profile_sha256=hashlib.sha256(args.profiles.read_bytes()).hexdigest(),
                    summary_sha256=SOURCE_HASH,cases={})
    for name,(fields,velocity,mismatch) in initial.items():
        rows = evolve(fields,velocity,grid,summary['params'],h)
        target = args.output/(name+'.csv')
        with target.open('w',newline='') as f:
            out = csv.writer(f,lineterminator='\n')
            out.writerow(HEADER)
            out.writerows([[format(v,'.17g') for v in row] for row in rows])
        metadata['cases'][name] = {'initial_cut_mismatch_B_Q_u':mismatch.tolist(),
                                  'rows':len(rows),'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':
    main()
