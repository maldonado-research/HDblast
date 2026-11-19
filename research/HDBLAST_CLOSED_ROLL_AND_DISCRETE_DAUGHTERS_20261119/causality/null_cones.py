"""Whole-slab numerical upper bounds on the radial coordinate null distance."""
import argparse
import json
import math
import os
import sys
from pathlib import Path
os.environ["OPENBLAS_NUM_THREADS"] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import DOP853
from scipy.special import eval_chebyu

root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"rollout"))
import rollout

def run(e, A, end, dt):
    m = rollout.Monitor(e,A,1e-6)
    z = e.initial(A,0.)
    s = DOP853(m.guarded_equation,0.,z,end,max_step=dt,rtol=3e-12,atol=2e-15)
    reach = speed = 0.
    while s.status=="running":
        s.step()
        if s.status=="failed":
            raise ArithmeticError("DOP853 failed")
        dense = s.dense_output()
        # Each coefficient F_j of SciPy's nested DOP853 representation
        # multiplies products from {r,1-r}, all in [0,1] over this step.
        abs_coeff = abs(dense.y_old.reshape(6,-1)[:2])+np.sum(abs(
                    dense.F.reshape(len(dense.F),6,-1)[:,:2]),axis=0)
        # |U_n(cos chi)| <= n+1.  The radial log-metric is m+sin²(chi)d.
        log_bound = abs_coeff.sum(axis=0)@(e.degrees+1)
        # a_bg(t) increases on the forward branch.
        max_speed = math.exp(log_bound)/e.background(s.t_old)[0]
        reach += (s.t-s.t_old)*max_speed
        speed = max(speed,max_speed)
    mcoef,dcoef,_,_,qcoef,_ = s.y.reshape(6,-1)
    endpoint_log_bound = (abs(mcoef)+abs(dcoef))@(e.degrees+1)
    length_lower = np.pi/2*e.background(end)[0]*math.exp(-endpoint_log_bound)
    target = 0. if A<0 else np.pi/2
    chi = np.linspace(target, target+(.3 if A<0 else -.3), 501)
    vals = eval_chebyu(e.degrees, np.cos(chi)[:,None])@qcoef
    error = abs(qcoef)@m.curvature_weights*(chi[1]-chi[0])**2/8
    return dict(roll_target_chi=target,guard_x=end,
                coordinate_null_reach_upper=reach, max_coordinate_speed_upper=speed,
                initial_domain_gap_lower=np.pi/2-2*reach,
                other_flank_initial_distance_lower=np.pi/2-reach,
                guard_proper_flank_separation_lower=length_lower,
                target_q=float(vals[0]), target_q_minus_minimum=float(vals[0]-m.target),
                target_radius_0p3_field_enclosure=[float(min(vals)-error),float(max(vals)+error)])

def calculate(gravity, rollout_data):
    arch = json.loads(rollout_data.read_text())
    runs=arch['branches']['phi_m=0.5']['runs']
    module=rollout.load(gravity)
    branch=module.constraint.load_branches(gravity/"previous")['phi_m=0.5']
    out={}
    for deg,step in [(112,.003),(192,.0015)]:
        e=module.Evolution(branch,deg,2*deg+32,6*deg+1)
        for a,r in runs.items():
            end=r["near_edge"]["end"]["x"]
            rec=run(e,float(a),end,step)
            out.setdefault(a,{})[str(deg)]=rec
            if rec['initial_domain_gap_lower']<.5 or rec['target_q_minus_minimum']<=.04:
                raise ArithmeticError("causal or field separation failed")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gravity', type=Path, default=root/'gravity')
    parser.add_argument('--rollout', type=Path, default=root/'rollout/ROLLOUT.json')
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('NULL.json'))
    args = parser.parse_args()
    out = calculate(args.gravity.resolve(), args.rollout)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
