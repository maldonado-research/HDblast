"""Mechanism screen, part 4: numerical calibration of a dissipative scalar junction (toy, flat 1D).

Doubled half-line z<0 (two copies), shell at z=0, kappa_5 = 1:
   bulk:     phi_tt = phi_zz - m^2 phi                      (m = 0: massless known limit;
             m = 3: bulk 'spring' that stores static gradient energy, shell settles at 1-phi^2 = m phi)
   junction: phi_z(0) = -(sigma'(phi_b) + Y phi_t(0))/2,    sigma = 2W(phi) = 2(1 - phi + phi^3/3)
   matter:   dE_m/dt = Y phi_t(0)^2    (j = Y v; radiation energy per unit shell volume, no expansion)
Energy E = 2*int (phi_t^2+phi_z^2+m^2(phi-phi_inf)^2)/2 dz + sigma(phi_b) + E_m is conserved.
Tests: (i) energy conservation (numerical control), (ii) massless case: matter fraction of the tension
drop = Y/(2+Y) (exact short-wave result E6), (iii) massive case: fraction differs, i.e.
the partition is regime-dependent and must be computed in the real warped bulk (motivates the
proposed next calculation). This is NOT the registered geometry: no warp, no expansion, flat bulk.
"""
import time, json
import numpy as np
from common import dump

def run(Y, m, Lz=40.0, N=16000, tf=30.0, phi0=0.0, phi_inf=0.0):
    z = np.linspace(-Lz, 0.0, N+1); h = z[1]-z[0]; dt = 0.4*h/max(1.0, Y/4)
    phi = np.full(N+1, phi0); pi = np.zeros(N+1)
    sig = lambda p: 2*(1 - p + p**3/3); sig1 = lambda p: 2*(p*p - 1)
    Em = 0.0; E0 = None; t = 0.0; hist = []
    def rhs(phi, pi):
        lap = np.empty_like(phi)
        lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2])/h**2
        # shell ghost: (g - phi[-2])/(2h) = phi_z(0)
        gz = -(sig1(phi[-1]) + Y*pi[-1])/2
        lap[-1] = (2*phi[-2] + 2*h*gz - 2*phi[-1])/h**2
        # far boundary: outgoing (left-moving) wave phi_t = phi_z  -> use one-sided
        lap[0] = 0.0
        acc = lap - m*m*(phi - phi_inf)
        dphi = pi.copy()
        dphi[0] = (phi[1] - phi[0])/h              # phi_t = phi_z at far end (first order)
        acc[0] = 0.0
        return dphi, acc
    def energy(phi, pi):
        grad = np.diff(phi)/h
        dens_k = pi**2; dens_m = m*m*(phi-phi_inf)**2
        eb = 2*(0.5*np.trapezoid(dens_k + dens_m, z) + 0.5*np.sum(grad**2)*h)
        return eb, sig(phi[-1])
    nsteps = int(tf/dt)
    eb, es = energy(phi, pi); E0 = eb + es
    out_flux = 0.0
    for n in range(nsteps):
        v0 = pi[-1]
        k1 = rhs(phi, pi); k2 = rhs(phi+dt/2*k1[0], pi+dt/2*k1[1]); k3 = rhs(phi+dt/2*k2[0], pi+dt/2*k2[1]); k4 = rhs(phi+dt*k3[0], pi+dt*k3[1])
        # energy leaving through the far boundary (both copies): 2 * (phi_t^2) at z=-Lz for left-moving waves
        out_flux += dt*2*k1[0][0]**2
        phi = phi + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); pi = pi + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        pi[0] = k4[0][0]
        Em += dt*Y*0.5*(v0**2 + pi[-1]**2)
        t += dt
    eb, es = energy(phi, pi)
    drop = sig(phi0) - es
    return dict(Y=Y, m=m, t_final=t, phi_b_final=float(phi[-1]), tension_drop=float(drop), matter_energy=float(Em),
                bulk_energy_in_box=float(eb), energy_left_box=float(out_flux),
                matter_fraction=float(Em/drop) if drop > 0 else None, short_wave_prediction=Y/(2+Y),
                energy_conservation_rel_error=float(abs(eb + es + Em + out_flux - E0)/max(drop, 1e-300)))

def main():
    t0 = time.time(); rows = []
    for m in [0.0, 3.0]:
        for Y in [0.5, 2.0, 8.0]:
            r = run(Y, m, phi_inf=0.0)
            rows.append(r); print(json.dumps(r), flush=True)
    # resolution control (massless, Y=2): halve N
    rc = run(2.0, 0.0, N=8000); rows.append(dict(resolution_control=True, **rc)); print(json.dumps(rc))
    checks = dict(
        massless_fraction_matches_Y_over_2_plus_Y=all(abs(r['matter_fraction']-r['short_wave_prediction']) < 5e-3
                                                     for r in rows if r['m'] == 0 and not r.get('resolution_control')),
        energy_conserved=all(r['energy_conservation_rel_error'] < 5e-3 for r in rows),
    )
    dump('M4_DISSIPATIVE_JUNCTION_TOY.json', dict(status='numerical toy (flat 1D, known-limit calibration)', rows=rows,
                                                  checks=checks, runtime_s=time.time()-t0))
    print(checks)

if __name__ == '__main__':
    main()
