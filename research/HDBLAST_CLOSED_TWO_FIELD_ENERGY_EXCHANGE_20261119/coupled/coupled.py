"""Closed, radial Einstein equations with two classical fields and a curved field metric."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre


def load(prior):
    path = prior/"gravity/evolution/evolve.py"
    spec = importlib.util.spec_from_file_location("closed_gravity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TwoFields:
    def __init__(self, module, branch, degree, nodes, matter):
        self.mod = module
        self.branch = branch
        self.e = module.Evolution(branch, degree, nodes, 4*degree+1)
        self.matter = matter
        self.ftop = branch.fitted.spline(branch.top)[1]
        self.domain = branch.fitted.theta[[0, -1]]-branch.top
        self.qcurvature = self.e.degrees*(self.e.degrees+1)*(self.e.degrees+2)/3
        self.angular_error_factor = self.qcurvature*(self.e.monitor_chi[1]**2/8)

    def coupling(self, q):
        if not self.matter:
            return np.ones_like(q), np.zeros_like(q)
        path = self.branch.fitted
        theta = self.branch.top+q
        f = path.spline(theta)[..., 1]
        if np.any(f <= 0):
            raise ArithmeticError("Non-positive field metric")
        return self.ftop/f, -path.spline(theta, 1)[..., 1]/f

    def initial(self, parent, seed, phase):
        """Solve the Hamiltonian constraint for the prescribed physical initial stress."""
        e = self.e
        chi, y, U, weights = e.chi, e.y, e.U, e.weight
        harmon = 4*y*y-1
        g2, du = self.mod.constraint.source(self.branch, parent, chi)
        z, _ = self.coupling(parent*harmon)
        state = np.zeros((8, e.n))
        state[4, 1] = parent
        if phase == "velocity":
            du = du + (seed*harmon)**2/6
            state[7] = e.project @ (seed*self.branch.h*harmon/np.sqrt(z))
        elif phase == "displacement":
            # Integrate dchi_field/dU2=(b/sqrt(8))/sqrt[Z(A U2)] from U2=0.
            points, wg = roots_legendre(24)
            sample = harmon[:, None]*(1+points[None, :])/2
            zpath, _ = self.coupling(parent*sample)
            field = seed*harmon*(1/np.sqrt(zpath) @ wg)/(2*np.sqrt(8))
            state[6] = e.project @ field
            g2 = g2 + (seed*8*y)**2*(1-y*y)/8
        else:
            raise ValueError("Unknown daughter phase")
        k2 = e.degrees*(e.degrees+2)
        c = np.zeros(e.n)
        for _ in range(15):
            s = U @ c
            res = self.mod.constraint.residual(s, U@(k2*c), g2, du)
            delta = e.project @ res
            if max(abs(delta)) < 2e-16:
                break
            derivative = 6-30*(1+s)**4*(1+du)-g2
            jac = np.diag(8*k2)+e.project @ (derivative[:, None]*U)
            c -= np.linalg.solve(jac, delta)
        else:
            raise ArithmeticError("Initial constraint failed to converge")
        if 1+c[0] <= (e.degrees[1:]+1)@abs(c[1:]):
            raise ArithmeticError("Non-positive conformal factor")
        state[0] = e.project @ (2*np.log1p(U @ c))
        return state.ravel()

    def equation(self, t, state):
        e = self.e
        c = state.reshape(8, e.n)
        # Certify the entire angular interpolation is within the tabulated model
        # at every trial RK stage; the second derivative bound includes the axes.
        field_grid = e.monitor @ c[4]
        error = abs(c[4]) @ self.angular_error_factor
        if min(field_grid)-error <= self.domain[0] or max(field_grid)+error >= self.domain[1]:
            raise ArithmeticError("A trial stage reached the calibrated field limit")
        m, d, beta_dev, shear, q, pi, chi, velocity = (e.U @ c.T).T
        my, dy, qy, chiy = (e.Uy @ c[[0, 1, 4, 6]].T).T
        myy, qyy, chiyy = (e.Uyy @ c[[0, 4, 6]].T).T
        fields = ((m, d, beta_dev, shear, q, pi), (my, dy, None, qy), (myy, qyy))
        _, H, inv, ly, _, _, _, dric, _ = e.geometry(t, fields, e.y, e.w)
        u, force = self.branch.potential(q)
        z, log_derivative = self.coupling(q)
        theta = 3*(H+beta_dev)+e.w*shear
        beta_dot = (-dric - 6*H*beta_dev - 3*beta_dev**2
                    - e.w*shear*(H+beta_dev) + 3*(u-self.branch.u_top))
        nonlinear_metric = 4*d*d*self.mod.exp_remainder(2*e.w*d)
        shear_dot = inv*(myy-ly*my+2*d+e.y*dy+nonlinear_metric+qy*qy+z*chiy*chiy)-theta*shear
        first = -3*e.y+e.w*(2*my-ly)
        phi_lap = inv*(e.w*qyy+first*qy)
        chi_lap = inv*(e.w*chiyy+first*chiy)
        phi_dot = phi_lap-theta*pi+force+0.5*z*log_derivative*(velocity**2-inv*e.w*chiy**2)
        chi_dot = chi_lap-theta*velocity+log_derivative*(inv*e.w*qy*chiy-pi*velocity)
        out = np.empty_like(c)
        out[[0, 1, 4, 6]] = c[[2, 3, 5, 7]]
        out[[2, 3, 5, 7]] = (e.project @ np.stack([beta_dot, shear_dot, phi_dot, chi_dot], axis=1)).T
        return out.ravel()

    def diagnostics(self, t, state, constraints=False):
        e, branch = self.e, self.branch
        c = state.reshape(8, e.n)
        U, Uy, Uyy = e.monitor_all if constraints else e.check
        y = np.cos(e.monitor_chi) if constraints else e.check_y
        w = 1-y*y
        m, d, b, shear, q, pi, chi, vel = (U @ c.T).T
        my, dy, by, qy, cy = (Uy @ c[[0, 1, 2, 4, 6]].T).T
        myy, qyy = (Uyy @ c[[0, 4]].T).T
        fields = ((m,d,b,shear,q,pi),(my,dy,by,qy),(myy,qyy))
        a, H, inv, _, F, _, _, _, dr3 = e.geometry(t, fields, y, w)
        u, _ = branch.potential(q)
        z, log_derivative = self.coupling(q)
        gsq = inv*w*qy*qy
        csq = inv*w*cy*cy
        beta = H+b
        alpha = beta+w*shear
        if constraints:
            ham = (dr3+12*H*b+6*b*b+4*w*beta*shear-pi*pi-gsq
                   -z*(vel*vel+csq)-6*(u-branch.u_top))
            mom = np.sqrt(inv*w)*(2*(by+F*shear)+pi*qy+z*vel*cy)
            rhs = self.equation(t,state).reshape(8,e.n)
            expansion_dot = 3*(U @ rhs[2]+branch.h**2/np.cosh(branch.h*t)**2)+w*(U @ rhs[3])
            ray = expansion_dot+alpha*alpha+2*beta*beta+pi*pi+z*vel*vel-3*u
            return np.array([max(abs(ham))/(6*branch.u_top),max(abs(mom))/branch.u_top,
                             max(abs(ray))/(3*branch.u_top)])
        vol = e.check_weight*(a*branch.h)**3*np.exp(3*m+w*d)
        weight = vol / branch.u_top
        kp, gp, kd, gd = pi*pi/6, gsq/6, z*vel*vel/6, z*csq/6
        eparent = weight @ (kp+gp+u)
        edaughter = weight @ (kd+gd)
        wp = weight @ (alpha*(kp+gp-u)+2*beta*(kp-gp-u))
        wd = weight @ (alpha*(kd+gd)+2*beta*(kd-gd))
        source = weight @ (z*log_derivative*pi*(csq-vel*vel)/6)
        shift_charge = vol @ (z*vel)/branch.h
        return np.array([eparent, edaughter, wp, wd, source, shift_charge,
                         vol.sum(),vol @ (alpha+2*beta)/(3*vol.sum()),
                         weight @ kd, weight @ gd])


def run_case(model, parent, seed, phase, end, step, rtol=2e-12):
    e = model.e
    initial = model.initial(parent, seed, phase)
    target = model.branch.fitted.tmin-model.branch.top
    sol = solve_ivp(model.equation, (0.,end), initial, method="DOP853", rtol=rtol,
                    atol=2e-15, max_step=step, dense_output=True)
    if not sol.success:
        raise ArithmeticError(sol.message)
    def event(t):
        coeff = sol.sol(t).reshape(8,e.n)[4]
        return max(e.monitor @ coeff) - target
    signed = [event(t) for t in sol.t]
    brackets = np.flatnonzero((np.asarray(signed[:-1]) < 0) & (np.asarray(signed[1:]) >= 0))
    if len(brackets) != 1:
        raise ArithmeticError("Expected one first-passage bracket in the calculated interval")
    index = brackets[0]
    crossing = brentq(event,sol.t[index],sol.t[index+1],xtol=1e-14)
    if crossing+.005 >= end:
        raise ArithmeticError("Event bracket reaches end of calculated interval")
    first = np.argmax(e.monitor @ sol.sol(crossing).reshape(8,e.n)[4])
    if first not in ((0,len(e.monitor_chi)-1) if parent>0 else (len(e.monitor_chi)//2,)):
        raise ArithmeticError("Unexpected first-contact site")
    times = dict(initial=0.,before=crossing-.005,crossing=crossing,after=crossing+.005,fixed=end)
    snapshots = {}
    for name,t in times.items():
        r = model.diagnostics(t,sol.sol(t))
        snapshots[name] = dict(t=t,E_parent=r[0],E_daughter=r[1],transfer_rate_into_daughter=r[4],
                               daughter_shift_charge=r[5],volume=r[6],H_volume=r[7],
                               E_daughter_kinetic=r[8],E_daughter_gradient=r[9],
                               daughter_fraction_of_total_scalar_energy=r[1]/(r[0]+r[1]),
                               daughter_pressure_to_energy=(r[8]-r[9]/3)/r[1] if r[1]>0 else None)
    def balance(start, stop, order):
        points, weights = roots_legendre(order)
        times = start+(stop-start)*(points+1)/2
        integrand = np.stack([model.diagnostics(t,sol.sol(t))[2:5] for t in times])
        wp, wd, transfer = (stop-start)/2 * (weights @ integrand)
        delta = model.diagnostics(stop,sol.sol(stop))[:2] - model.diagnostics(start,sol.sol(start))[:2]
        return dict(parent_energy_change=delta[0],daughter_energy_change=delta[1],
                    parent_pressure_work=wp,daughter_pressure_work=wd,transfer_into_daughter=transfer,
                    parent_balance_residual=delta[0]+wp+transfer,
                    daughter_balance_residual=delta[1]+wd-transfer,
                    total_balance_residual=sum(delta)+wp+wd)
    ranges = dict(to_crossing=(0.,crossing),through_crossing=(crossing-.005,crossing+.005),to_fixed=(0.,end))
    balances = {name:balance(start,stop,56) for name,(start,stop) in ranges.items()}
    fine = balance(0.,end,80)
    quadrature = max(abs(fine[k]-balances['to_fixed'][k]) for k in fine)
    peak = np.max([model.diagnostics(t,sol.sol(t),constraints=True)
                  for t in np.linspace(0.,end,31)],axis=0)
    charges = [model.diagnostics(t,sol.sol(t))[5] for t in np.linspace(0.,end,31)]
    return dict(parent_amplitude=parent,daughter_seed=seed,launch_phase=phase,matter_weight=model.matter,
                degree=e.degree,angular_nodes=len(e.chi),maximum_step=step,snapshots=snapshots,balances=balances,
                constraints=dict(hamiltonian=peak[0],momentum=peak[1],raychaudhuri=peak[2]),
                max_abs_shift_charge_drift=max(abs(np.asarray(charges)-charges[0])),
                time_quadrature_absolute_change=quadrature)


def compare_runs(a,b):
    return dict(crossing_time_delta=a["snapshots"]["crossing"]["t"]-b["snapshots"]["crossing"]["t"],
                fixed_H_volume_delta=a["snapshots"]["fixed"]["H_volume"]-b["snapshots"]["fixed"]["H_volume"],
                fixed_volume_delta=a["snapshots"]["fixed"]["volume"]-b["snapshots"]["fixed"]["volume"],
                fixed_transfer_delta=a["balances"]["to_fixed"]["transfer_into_daughter"]-b["balances"]["to_fixed"]["transfer_into_daughter"],
                fixed_daughter_energy_delta=a["snapshots"]["fixed"]["E_daughter"]-b["snapshots"]["fixed"]["E_daughter"])


def calculate(prior, short=False):
    module = load(prior)
    branch = module.constraint.load_branches(prior/"gravity/previous")["phi_m=-0.5"]
    configs = [("primary",64,192,.004),("space",96,256,.004),("time",64,192,.002)]
    output = dict(branch="phi_m=-0.5",hilltop_H=branch.h,normalization="E and integrated work divided by V_round,0 times u_top",
                  initial_data="matched local daughter stress in the two matter couplings, parent U2, vanishing extrinsic curvature and parent velocity",
                  runs={},refinement={})
    for tag,degree,nodes,step in configs[:1] if short else configs:
        models = {matter:TwoFields(module,branch,degree,nodes,matter) for matter in (True,False)}
        for parent in [.001,-.001]:
            for phase in ("velocity","displacement"):
                for seed in (.01,.005,0.):
                    for matter,model in models.items():
                        if seed == 0. and (not matter or phase != 'velocity'):
                            continue
                        if tag != 'primary' and seed != .01:
                            continue
                        end = .505 if parent>0 else .605
                        key=f"{parent:g}/{phase}/{seed:g}/"+("matter" if matter else "constant")
                        case=run_case(model,parent,seed,phase,end,step)
                        if tag=="primary":
                            output["runs"][key]=case
                        else:
                            output["refinement"].setdefault(key,{})[tag]=dict(
                                differences=compare_runs(case,output["runs"][key]),
                                constraints=case["constraints"],
                                max_balance=max(abs(v[k]) for v in case["balances"].values()
                                                for k in ('parent_balance_residual','daughter_balance_residual')))
    runs=output['runs']
    comparisons={}
    for parent in (.001,-.001):
        for phase in ("velocity","displacement"):
            for seed in (.01,.005):
                key=f"{parent:g}/{phase}/{seed:g}"
                matter=runs[key+"/matter"]; constant=runs[key+"/constant"]
                base=runs[f"{parent:g}/velocity/0/matter"]
                comparisons[key]=dict(matter_minus_constant=compare_runs(matter,constant),
                                      matter_minus_unseeded=compare_runs(matter,base),
                                      initial_energy=matter["snapshots"]["initial"]["E_daughter"],
                                      initial_metric_mismatch=matter["snapshots"]["initial"]["volume"]-constant["snapshots"]["initial"]["volume"],
                                      initial_daughter_energy_mismatch=matter["snapshots"]["initial"]["E_daughter"]-constant["snapshots"]["initial"]["E_daughter"])
    output["paired_controls"]=comparisons
    summary = {}
    for key, control in comparisons.items():
        run=runs[key+"/matter"]
        launch=control["initial_energy"]
        event=run["balances"]["to_crossing"]
        short_interval=run["balances"]["through_crossing"]
        summary[key]=dict(transfer_to_crossing_as_percent_of_daughter_seed=100*event["transfer_into_daughter"]/launch,
                          energy_change_to_crossing_as_percent_of_daughter_seed=100*event["daughter_energy_change"]/launch,
                          transfer_through_crossing_as_percent_of_daughter_seed=100*short_interval["transfer_into_daughter"]/launch,
                          crossing_total_daughter_share=run["snapshots"]["crossing"]["daughter_fraction_of_total_scalar_energy"],
                          crossing_daughter_pressure_to_energy=run["snapshots"]["crossing"]["daughter_pressure_to_energy"])
    output["summary"]=summary
    for case in runs.values():
        if (max(case["constraints"].values())>2e-7 or case["max_abs_shift_charge_drift"]>1e-10
            or any(abs(row[k])>1e-9 for row in case["balances"].values() for k in ("parent_balance_residual","daughter_balance_residual","total_balance_residual"))
            or case["time_quadrature_absolute_change"]>1e-9):
            raise ArithmeticError(f"Physical controls failed for {case['parent_amplitude']}, {case['launch_phase']}, {case['daughter_seed']}: {case['constraints']}")
        if not case["matter_weight"] and any(r["transfer_into_daughter"] != 0 for r in case["balances"].values()):
            raise ArithmeticError("Geometric control shows direct transfer")
    for parent in (.001,-.001):
        for phase in ("velocity","displacement"):
            key=f"{parent:g}/{phase}"
            full=summary[key+"/0.01"]["transfer_to_crossing_as_percent_of_daughter_seed"]
            half=summary[key+"/0.005"]["transfer_to_crossing_as_percent_of_daughter_seed"]
            expected_sign=(1 if parent>0 else -1)*(1 if phase=="velocity" else -1)
            if expected_sign*half<=0 or abs(full/half-1)>.001:
                raise ArithmeticError("Transfer sign or matched weak-seed scaling failed")
        zero=runs[f"{parent:g}/velocity/0/matter"]
        if any(row["E_daughter"]!=0 for row in zero["snapshots"].values()):
            raise ArithmeticError("Unseeded classical daughter moved")
    for pair in comparisons.values():
        if max(abs(pair[k]) for k in ("initial_metric_mismatch","initial_daughter_energy_mismatch"))>2e-13:
            raise ArithmeticError("Unequal initial stress in paired runs")
    for tests in output["refinement"].values():
        for row in tests.values():
            if max(abs(v) for v in row["differences"].values())>3e-8 or row["max_balance"]>1e-9:
                raise ArithmeticError("Unconverged refinement")
    output["passed"]=True
    return output


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior",type=Path,default=Path(__file__).resolve().parents[1]/"prior")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--short",action="store_true")
    args=parser.parse_args()
    out=calculate(args.prior,args.short)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
