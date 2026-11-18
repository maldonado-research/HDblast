"""Energy, scalar source and phase rotation for a finite comoving band on a rolling FRW trajectory."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

START = 0.80
CUT = 1.55
END = 2.10
BAND = 8.0
G = 100.0
ARMS = ('quiet', 'near_positive', 'near_negative', 'far_positive')


def histories(source, grid, arm):
    parent = source / 'classical' / 'runs' / grid
    branch = {'quiet': 'quiet', 'near_positive': '0.001', 'near_negative': '-0.001'}.get(arm)
    folder = parent / branch if branch else source / 'far' / 'runs' / grid / '0.001'
    with np.load(parent / 'checkpoint.npz', allow_pickle=False) as checkpoint:
        early = checkpoint['hist']
    summary = json.loads((parent / 'quiet/rolling_summary.json').read_text())
    normalization = summary['matter']['events'][0]['A_star']
    with np.load(folder / 'rolling_timeseries.npz', allow_pickle=False) as saved:
        late = dict(zip(saved['cols'], saved['rec'].T, strict=True))
    t0 = late['H0tau'][0]
    select = early[:, 0] < t0 - 1e-10
    times = np.r_[early[select, 0], late['H0tau']]
    phi = CubicSpline(times, np.r_[early[select, 1], late['phi_b']])
    hubble = CubicSpline(times, np.r_[early[select, 2], late['H_over_H0']])
    ln_a = hubble.antiderivative()
    offset = late['A_b'][0] - normalization - ln_a(t0)
    return phi, hubble, lambda t: np.exp(offset + ln_a(t))


def basis(t, k, phi, h, a):
    """Common order-two adiabatic state for X=a^(3/2) chi in the preceding study."""
    eps = .001
    def omega(time):
        return np.sqrt((k / a(time))**2 + G**2 * (phi(time) - .5)**2
                       - 2.25 * h(time)**2 - 1.5 * h(time, 1))
    def w(time):
        low, mid, high = (omega(time + s) for s in (-eps, 0, eps))
        first = (high - low) / (2*eps)
        second = (high - 2*mid + low) / eps**2
        return mid - second / (4*mid**2) + 3*first**2 / (8*mid**3)
    freq = w(t)
    slope = (w(t+eps) - w(t-eps)) / (2*eps)
    x = (2*freq)**(-.5)
    xp = (-1j*freq - slope/(2*freq)) * x
    return x/a(t)**1.5, a(t)**1.5 * (xp - 1.5*h(t)*x)


class Experiment:
    def __init__(self, source, grid, arm, nodes, rtol):
        self.phi, self.hubble, self.scale = histories(source, grid, arm)
        x, weights = np.polynomial.legendre.leggauss(nodes)
        self.k = (x+1)*BAND/2
        self.weight = weights*BAND/2 * self.k**2/(2*np.pi**2)
        self.phi0 = float(self.phi(START))
        self.rtol = rtol

    def background(self, s):
        a = self.scale(s)
        phi = self.phi(s)
        h = self.hubble(s)
        u = phi - .5
        q2 = (self.k/a)**2
        m2 = G**2*u**2
        return a, phi, h, u, q2, m2, np.sqrt(q2+m2)

    def powers(self, s, q, p):
        a, phi, h, u, q2, m2, f = self.background(s)
        v = self.phi(s, 1)
        volume = a**3
        qq, pp = abs(q)**2, abs(p)**2
        energy = .5 * (pp/volume + volume*(q2+m2)*qq)
        pressure = .5 * (pp/volume - volume*(q2/3+m2)*qq)
        force = G**2*u*volume*qq
        pressure0 = q2/(6*f)
        force0 = G**2*u/(2*f)
        # delta E=a^3 C(phi), with C even about phi=0.5, C(phi0)=C'(phi0)=0.
        z = u**2 - (self.phi0-.5)**2
        pot = z*z
        dpot = 4*u*z
        total = self.weight @ energy
        bare_ground = self.weight @ (f/2)
        rates = np.array([v*(self.weight @ force), -3*h*(self.weight @ pressure),
                          v*(self.weight @ force0), -3*h*(self.weight @ pressure0),
                          volume*dpot*v, 3*h*volume*pot])
        return total, bare_ground, volume*pot, rates

    def propagate(self, start, end, q, p):
        n = len(q)
        def rhs(s, state):
            q, p = state[:n], state[n:2*n]
            volume = self.scale(s)**3
            freq2 = (self.k/self.scale(s))**2 + G**2*(self.phi(s)-.5)**2
            return np.r_[p/volume, -volume*freq2*q, self.powers(s, q, p)[-1]]
        state = np.r_[q, p, np.zeros(6)]
        sol = solve_ivp(rhs, (start,end), state, method='DOP853', rtol=self.rtol,
                        atol=self.rtol*1e-2, max_step=.005)
        if not sol.success:
            raise ValueError(sol.message)
        qf, pf = sol.y[:n,-1], sol.y[n:2*n,-1]
        e, e0, ct, _ = self.powers(end, qf, pf)
        ei, e0i, cti, _ = self.powers(start, q, p)
        w, g, w0, g0, wc, gc = sol.y[2*n:,-1].real
        exc = e-e0
        wr = float(np.max(abs(1j*(qf.conj()*pf-pf.conj()*qf)-1)))
        if wr > 1e-7 or exc < -1e-8:
            raise ValueError('Mode normalization or excitation energy')
        ledger = dict(energy_change=float(e-ei), mass_work=w, expansion_work=g,
                      ground_energy_change=e0-e0i, ground_mass_work=w0, ground_expansion_work=g0,
                      exc_energy_change=exc-(ei-e0i), exc_mass_work=w-w0, exc_expansion_work=g-g0,
                      raw_defect=e-ei-w-g, exc_defect=exc-(ei-e0i)-(w-w0)-(g-g0),
                      end_exc_energy=exc, end_raw_energy=e,
                      potential_cnumber_work_unit=wc, potential_cnumber_expansion_unit=gc,
                      potential_cnumber_energy_change_unit=ct-cti,
                      potential_cnumber_defect=ct-cti-wc-gc, wronskian_defect=wr)
        return qf, pf, ledger

    def future_envelopes(self, n, coherence, cross, actual, rotated):
        """Transport all quadratic forms for arbitrary phases at fixed occupation."""
        a = self.scale(CUT)
        f = np.hypot(self.k/a, G*(self.phi(CUT)-.5))
        M = a**3*f
        zeros = np.zeros_like(self.k)
        # Rows 0..3 are both real columns of the transfer from normalized (x,y).
        # Rows 4..9 integrate the 11,22,12 elements for mass and pressure work.
        initial = np.stack((M**(-.5), zeros, zeros, M**.5, *([zeros]*6)))
        def rhs(s, state):
            X1, X2, P1, P2 = state.reshape(initial.shape)[:4]
            a, phi, h, u, q2, m2, f = self.background(s)
            V = a**3
            xx = np.stack((X1*X1, X2*X2, X1*X2))
            pp = np.stack((P1*P1, P2*P2, P1*P2))
            wm = G**2*u*self.phi(s,1)*V*xx
            wg = -1.5*h*(pp/V - V*(q2/3+m2)*xx)
            return np.vstack((P1/V,P2/V,-V*f*f*X1,-V*f*f*X2,wm,wg)).ravel()
        sol = solve_ivp(rhs, (CUT,END), initial.ravel(), method='DOP853',
                        t_eval=[END], rtol=self.rtol, atol=self.rtol*.01, max_step=.005)
        if not sol.success:
            raise ValueError(sol.message)
        out = sol.y[:,-1].reshape(initial.shape)
        X1,X2,P1,P2 = out[:4]
        symplectic_defect = float(np.max(abs(X1*P2-X2*P1-1)))
        if symplectic_defect > 1e-7:
            raise ValueError('Transfer determinant')
        a = self.scale(END)
        f_end = np.hypot(self.k/a, G*(self.phi(END)-.5))
        xx = np.stack((X1*X1,X2*X2,X1*X2))
        pp = np.stack((P1*P1,P2*P2,P1*P2))
        energy_form = (pp/a**3+a**3*f_end**2*xx)/2
        initial_form = np.stack((f/2,f/2,zeros))
        forms = {'mass_work':out[4:7], 'expansion_work':out[7:10],
                 'energy_change':energy_form-initial_form}
        report = {'transfer_symplectic_defect':symplectic_defect,
                  'ward_matrix_defect':float(np.max(abs(forms['mass_work']+forms['expansion_work']-forms['energy_change'])))}
        # A fixed n permits every point on a circle of radius sqrt(n(n+1)).
        for label, K in forms.items():
            diagonal = self.weight @ ((n+.5)*(K[0]+K[1]))
            vector = np.stack((K[0]-K[1],2*K[2]))
            length = np.hypot(vector[0],vector[1])
            radius = np.sqrt(np.maximum(0,n*(n+1)))
            half = self.weight @ (radius*length)
            offset = self.weight @ (coherence*vector[0]+cross*vector[1])
            observed = diagonal+offset
            twin = diagonal-offset
            largest_replay_error = max(abs(observed-actual[label]),abs(twin-rotated[label]))
            if abs(offset)>half+1e-7 or largest_replay_error>2e-4:
                raise ValueError(('Work support function disagrees with modes',label,largest_replay_error))
            # Construct a pure extremizing state independently and propagate its endpoint energy.
            angle = np.arctan2(vector[1],vector[0])
            fmax = radius*np.cos(angle)
            cmax = radius*np.sin(angle)
            covariance_xx = n+.5+fmax
            extremal_q = np.sqrt(covariance_xx)
            extremal_p = (cmax-.5j)/extremal_q
            kcheck = self.weight @ (K[0]*abs(extremal_q)**2 + K[1]*abs(extremal_p)**2
                                     + 2*K[2]*np.real(extremal_q.conj()*extremal_p))
            if abs(kcheck-(diagonal+half))>1e-7:
                raise ValueError('Unattained occupation-conditioned envelope')
            report[label] = dict(center=diagonal,half_range=half,minimum=diagonal-half,
                                 maximum=diagonal+half,actual=observed,quarter_turn=twin,
                                 max_direct_replay_error=largest_replay_error)
        return report

    def run(self):
        q0, p0 = basis(START, self.k, self.phi, self.hubble, self.scale)
        q, p, early = self.propagate(START, CUT, q0, p0)
        a, phi, h, u, q2, m2, f = self.background(CUT)
        # A quarter turn in (sqrt(M) q, p/sqrt(M)) where M=a^3 omega.
        M = a**3*f
        rq, rp = p/M, -M*q
        def numbers(x,y):
            return 0.5*(M*abs(x)**2+abs(y)**2/M-1)
        occ = numbers(q,p)
        qq = .5*(M*abs(q)**2-abs(p)**2/M)
        qp = np.real(q.conj()*p)
        radius = np.sqrt(np.maximum(0,occ*(1+occ)))
        rform = np.hypot(qq,qp)
        err = float(np.max(abs(radius**2-rform**2)))
        same_n = float(np.max(abs(occ-numbers(rq,rp))))
        coeff = G**2 * u / f
        force = self.weight @ (coeff*(occ+.5+qq))
        rotated = self.weight @ (coeff*(occ+.5-qq))
        diagonal = self.weight @ (coeff*(occ+.5))
        span = self.weight @ (abs(coeff)*radius)
        if err > 1e-7 or same_n > 2e-12 or abs(force-diagonal) > span+1e-8:
            raise ValueError('Coincident spectrum, covariance purity or force envelope')
        qf, pf, orig = self.propagate(CUT, END, q, p)
        rf, rpf, turn = self.propagate(CUT, END, rq, rp)
        a2 = self.scale(END)
        f2 = np.hypot(self.k/a2, G*(self.phi(END)-.5))
        def n_final(x,y):
            return self.weight @ (.5*(a2**3*f2*abs(x)**2+abs(y)**2/(a2**3*f2)-1))
        mix = {k:.5*(orig[k]+turn[k]) for k in ('mass_work','expansion_work','end_raw_energy','end_exc_energy')}
        contrast = {k:orig[k]-turn[k] for k in ('energy_change','mass_work','expansion_work','end_exc_energy')}
        envelope = self.future_envelopes(occ,qq,qp,orig,turn)
        return dict(early=early, future_envelopes=envelope, cut=dict(phi=float(phi), H=float(h), scale=float(a),
                    mass_over_H0=float(np.sqrt(m2)), band_number=float(self.weight@occ),
                    max_mode_n_gap=same_n, purity_gap=err, force_original=force,
                    force_quarter_turn=rotated, force_diagonal=diagonal,
                    force_ground=self.weight@(coeff/2), sharp_force_half_range=span,
                    energy_gap_at_rotation=float(self.powers(CUT,q,p)[0]-self.powers(CUT,rq,rp)[0])),
                    original=orig, quarter_turn=turn, dephased_mixture=mix,
                    original_minus_quarter=contrast,
                    original_final_inst_number=n_final(qf,pf), rotated_final_inst_number=n_final(rf,rpf))


def compare_arms(arms):
    """Original state only. The example source-potential convention shifts by lambda times the last columns."""
    def total(result, k):
        return result['early'][k]+result['original'][k]
    keys = ('mass_work','expansion_work','exc_mass_work','exc_expansion_work',
            'energy_change','exc_energy_change','potential_cnumber_work_unit',
            'potential_cnumber_energy_change_unit','potential_cnumber_expansion_unit')
    return {arm:{k:total(arms[arm],k)-total(arms['quiet'],k) for k in keys} for arm in ARMS[1:]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prior', type=Path)
    opt = parser.parse_args()
    reports = {}
    for grid, nodes, rtol in [('coarse',32,2e-10),('fine',32,2e-10),('fine',48,4e-11)]:
        result = {arm:Experiment(opt.prior,grid,arm,nodes,rtol).run() for arm in ARMS}
        reports[f'{grid}_{nodes}'] = dict(arms=result, pulse_minus_quiet=compare_arms(result))
    print(json.dumps(reports,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
