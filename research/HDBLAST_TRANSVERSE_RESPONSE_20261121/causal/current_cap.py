"""Maximum sideways force from a seedless Fourier pair driven by the wall particles."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import kve

GAUGE = 2.54


def vector_kernel_limit():
    mass2 = 8300.**2
    width = .30
    half = mass2/2.
    below = 24.*math.pi*width**2*(half+2.)**2/(half*(mass2-half)**2)
    delta = width**2*(1.-8./mass2)
    above = 12.*math.pi/mass2*width**2*(1.+4./mass2)**2*(math.sqrt(1.+delta)+1.)/delta
    return max(below,above)


def thermal_cdf(p, temp):
    w = np.hypot(1.,p)
    log_survival = 2.*np.log(w)+(1.-w)/temp+np.log(kve(2,w/temp)/kve(2,1./temp))
    return np.maximum(0.,-np.expm1(log_survival))


def max_gap(floors, threshold, temp):
    """Difference between guaranteed CDF jumps and the desired continuous CDF."""
    deficit = floors-thermal_cdf(threshold,temp)
    index = int(np.argmax(deficit))
    return max(0.,float(deficit[index])), index


def smallest_count(floors, tube, kicks, temp, tolerance):
    def insufficient(m):
        return max_gap(floors,tube+m*kicks,temp)[0]>tolerance
    nonpass, upper = 0, 1
    if not insufficient(0):
        return 0
    while insufficient(upper):
        nonpass,upper=upper,2*upper
        if upper>1<<26:raise RuntimeError("Unresolved Fourier multiplicity")
    while upper-nonpass>1:
        mid=(upper+nonpass)//2
        if insufficient(mid):nonpass=mid
        else:upper=mid
    return upper


def analyze(trace, marker, refs, branch, charge, mass, action, enclose_marker=False):
    end=marker["time"]
    earlier=trace["t"]<end
    right_edge = np.searchsorted(trace["t"],end)
    def series(name):
        if enclose_marker:
            return trace[name][:right_edge+1]
        return np.append(trace[name][earlier],np.interp(end,trace["t"],trace[name]))
    t,a,d=map(series,("t","a","D"))
    dt=np.diff(t)
    if np.any(dt<=0.) or np.any(np.diff(a)<0.) or np.any(np.diff(d)<0.):
        raise ValueError("Monotone radius and birth count required")
    # D/a <= D_right/a_left on each archived cell. Partial inner integrals
    # have a linear upper envelope, which the outer integral integrates exactly.
    inner_slope=d[1:]/a[:-1]
    inner_left=np.r_[0.,np.cumsum(inner_slope*dt)[:-1]]
    endpoint_radius_lower=a[-2] if enclose_marker else a[-1]
    impulse_cell=2.*charge**2/(GAUGE*endpoint_radius_lower)*(inner_left+.5*inner_slope*dt)*dt/a[:-1]
    kicks=np.cumsum(impulse_cell[::-1])[::-1]
    # For a marker known only to fall in the final cell, use opposite endpoint
    # enclosures for upper RMS and lower cumulative participation fractions.
    upper_weights=np.diff(d)/(d[-2] if enclose_marker else d[-1])
    ever_born=np.diff(d)/d[-1]
    if enclose_marker:ever_born[-1]=0.
    tube=np.sqrt(max(0.,float(np.max(series("R")))**2-1.))
    temp=min(series("theta")[-2:]) if enclose_marker else marker["theta"]
    # R and the earlier electric impulse sets horizontal birth support. The
    # strongest global-field upper is compatible with all later source kicks
    # being in phase at every observation location.
    upper_rms=float(np.sqrt(upper_weights@(kicks*kicks)))
    smooth_inner=cumulative_trapezoid(d/a,t,initial=0.)
    instantaneous_field=2.*charge/GAUGE*smooth_inner/a**2
    smooth_prefix=charge/a[-1]*cumulative_trapezoid(a*instantaneous_field,t,initial=0.)
    smooth=smooth_prefix[-1]-smooth_prefix
    smooth_rms=float(np.sqrt(ever_born@(.5*(smooth[:-1]**2+smooth[1:]**2))))
    hazard_cell=branch*action*vector_kernel_limit()/(6.*mass**3)*dt*d[1:]/a[:-1]**3
    earliest_survival=np.exp(-np.cumsum(hazard_cell[::-1])[::-1])
    total_floor=np.cumsum(ever_born[::-1])[::-1]
    unhit_floor=np.cumsum((ever_born*earliest_survival)[::-1])[::-1]
    hit_floor,worst_hit=max_gap(total_floor,tube+kicks,temp)
    conditional_gap,idx=max_gap(unhit_floor,tube+kicks,temp)
    counts={str(tol):smallest_count(unhit_floor,tube,kicks,temp,tol) for tol in (0.,.001,.01,.05,.1)}
    if enclose_marker:
        return {"rms_upper":upper_rms,"required_hit_fraction_lower":hit_floor,
                "cdf_gap_lower":conditional_gap,"harmonics_within_one_point":counts["0.01"]}
    out={"time":end,"temperature":temp,"birth_radius":float(tube),"birth_cells":len(kicks),
         "field_amplitude_upper_at_marker":float(2.*charge/GAUGE/a[-1]**2*sum(inner_slope*dt)),
         "maximum_impulse_upper":float(kicks[0]),"population_rms_impulse_upper":upper_rms,
         "linear_history_trapezoid_rms":smooth_rms,
         "minimum_distinct_hit_fraction_required_for_thermal_cdf":hit_floor,
         "cutoff_for_required_hits":float(tube+kicks[worst_hit]),
         "unhit_probability_floor_under_homogeneous_hazard":float(unhit_floor[0]),
         "mandatory_cdf_excess_under_homogeneous_hazard":conditional_gap,"momentum_at_excess":float(tube+kicks[idx]),
         "minimum_harmonics_for_cdf_tolerance_under_homogeneous_hazard":counts}
    if refs is not None:
        target=refs["rms_extra_impulse_lower_over_m_all_present"]
        out.update(previous_conditional_rms_impulse_lower=target,target_over_cap=target/upper_rms,
                   minimum_harmonics_from_rms=math.ceil(target/upper_rms))
    return out


def real_mode_control():
    angle=np.linspace(0.,2.*math.pi,49)
    largest=0.
    rng=np.random.default_rng(205621)
    local_tight=0.
    for _ in range(100):
        khat=rng.normal(size=3); khat/=np.linalg.norm(khat)
        # The guide and the particle's transverse radial direction can be z and x.
        ex=np.array([1.,0.,0.]); ez=np.array([0.,0.,1.])
        na=ex@khat; za=ez@khat; radial_l=na*khat
        for delta in angle:
            A=radial_l+math.cos(delta)*(ex-radial_l)
            B=math.sin(delta)*(na*ez-za*ex)
            gram=np.array([[A@A,A@B],[A@B,B@B]])
            largest=max(largest,float(np.linalg.eigvalsh(gram)[1]))
            check=math.sin(delta)**2*(1.-na**2-za**2)*(1.-math.sin(delta)**2*na**2)
            local_tight=max(local_tight,abs(np.linalg.det(np.eye(2)-gram)-check))
    if largest>1.+1.e-12 or local_tight>1.e-12:raise ValueError("Single-source radial bound failed")
    return {"largest_squared_operator_norm":largest,"largest_absolute_determinant_residual":local_tight}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args(); meta=json.loads((args.history/"markers.json").read_bytes())
    # Data from the published November 21 parent study, ballistic/KICK.json.
    reference={"dedicated_half":541.645547427229,"dedicated_native":540.9609444060337,
               "shared_half":674.0982875141613,"shared_native":673.7102609380573}
    runs={}
    for name,entry in meta["runs"].items():
        with np.load(args.history/(name+".npz")) as data:
            trace={k:data[k] for k in ("t","a","D","R","theta")}
        runs[name]={}
        for mark,marker in entry["markers"].items():
            refs={"rms_extra_impulse_lower_over_m_all_present":reference[name]} if mark=="0.01" else None
            params=(trace,marker,refs,entry["channel_fraction"],meta["charge_to_mass"],meta["daughter_mass_L"],meta["action"])
            row=analyze(*params)
            row["enclosing_entire_marker_cell"]=analyze(*params,enclose_marker=True)
            runs[name][mark]=row
    args.output.write_text(json.dumps({"controls":real_mode_control(),"runs":runs},indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":main()
