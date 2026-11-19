"""Compute switched pair yields on the returned walls and retained-wall controls."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("pair_spectrum",HERE/"pair_spectrum.py")
pairs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pairs)
ARMS = ("quiet", "plus", "minus", "far")


def csv_data(stream):
    return np.genfromtxt(stream, delimiter=",", names=True)


def fork_states(fork):
    q = np.asarray(fork["q"])
    weights = np.asarray(fork["weight"])
    x,w = np.polynomial.legendre.leggauss(len(q))
    if not np.allclose(q,4*(1+x),rtol=1e-13,atol=1e-13) or not np.allclose(
            weights,4*w*q**2/(2*np.pi**2),rtol=2e-12,atol=1e-15):
        raise ValueError("Native phase-space normalization")
    states = {}
    for label, row in fork["states"].items():
        X,P,C = (np.asarray(row[k]) for k in ("X","P","C"))
        n,c = pairs.source_covariance(q,fork["a"],100*(fork["phi_b"]-.5),X,P,C)
        if not np.allclose(row["n"],n,atol=2e-13) or not np.allclose(row["Re_c"],c.real,atol=2e-13):
            raise ValueError("Fork state normalization")
        states[label] = (n,c)
    return q,weights,states


def check_source_equivalence(prior_walls, grid, post, ref, q, weights, result):
    # Without a fork operation the original preparation was a vacuum. Evolving its
    # mode through the continuous y=0 wall must reproduce the count reconstructed
    # from the retained fork density matrix. Test both windows on each pulse arm.
    with prior_walls.open(f"{grid}/y0_b1/pre.csv") as source:
        pre = csv_data(source)
    common = ("H0tau","A_b","phi_b","H_over_H0","v_over_H0")
    history = {k: np.r_[pre[k][pre["H0tau"] < post["H0tau"][0]],post[k]] for k in common}
    start,bg = pairs.wall_interpolants(history,ref)
    resp = pairs.fourier_response(start,bg,q,1024,32)
    vacuum = pairs.evaluate(resp,np.zeros_like(q),np.zeros_like(q),weights,np.exp(bg[1](pairs.CUT)))
    return max(abs(vacuum[w][b][k]["total"]/outputs["total"]-1)
               for w,windows in result.items() for b,bands in windows.items() for k,outputs in bands.items())


def max_solver_difference(one,two):
    # Express both total and excess uncertainties against the positive total. A
    # matched-vacuum excess can approach zero, making its relative error singular.
    return max(abs(values[key]-two[window][band][observable][key])/values["total"]
               for window,windows in one.items() for band,bands in windows.items()
               for observable,values in bands.items() for key in ("total","excess"))


def contrast_summary(walls):
    # Normalize every difference by N or E of the same preparation on its
    # retained quiet wall. This remains positive when the excess crosses zero.
    output = {}
    for state in walls:
        impacts = {}
        for arm in ARMS[1:]:
            windows = {}
            for win in pairs.WINDOWS:
                windows[win] = {}
                for band in ("full","screen"):
                    table = {}
                    for obs in ("number","energy"):
                        fixed = walls["retain"][arm]["states"][state][win][band][obs]
                        quiet = walls["retain"]["quiet"]["states"][state][win][band][obs]
                        active = walls[state][arm]["states"][state][win][band][obs]
                        active_quiet = walls[state]["quiet"]["states"][state][win][band][obs]
                        for key in ("total","excess","total_density"):
                            scale = quiet["total_density" if key=="total_density" else "total"]
                            fixed_signal = 100*(fixed[key]-quiet[key])/scale
                            active_signal = 100*(active[key]-active_quiet[key])/scale
                            table[obs+"_"+key] = dict(fixed=fixed_signal, responding=active_signal,
                                geometry_effect=active_signal-fixed_signal)
                    windows[win][band] = table
            impacts[arm] = windows
        output[state] = impacts
    return output


def midpoint_test(walls):
    # All three preparations share n and Re(c). At fixed wall, N is affine
    # in Im(c), so the second difference in C must vanish for every pulse.
    output = {}
    for arm in ARMS:
        windows = {}
        for window in pairs.WINDOWS:
            bands = {}
            for band in ("full","screen"):
                result = {}
                for observable in ("number","energy"):
                    for key in ("total","excess","total_density"):
                        label = observable+"_"+key
                        values = []
                        for geom_fixed in (True,False):
                            x = [walls["retain" if geom_fixed else state][arm]["states"][state]
                                 [window][band][observable][key] for state in ("retain","erase","reflect")]
                            reference = walls["retain"][arm]["states"]["retain"][window][band][observable]
                            unit = reference["total_density" if key=="total_density" else "total"]
                            values.append(float((x[0]-2*x[1]+x[2])/unit))
                        result[label] = {"fixed":values[0],"responding":values[1]}
                bands[band] = result
            windows[window] = bands
        output[arm] = windows
    return output


def midpoint_mismatch(first, second):
    return max(abs(v["responding"]-second[arm][win][band][obs]["responding"])
               for arm,windows in first.items() for win,bands in windows.items()
               for band,row in bands.items() for obs,v in row.items())


def contrast_errors(first, second):
    absolute = {k:0. for k in ("fixed","responding","geometry_effect")}
    for state,arms in first.items():
        for arm,windows in arms.items():
            for win,bands in windows.items():
                for band,observables in bands.items():
                    for k,components in observables.items():
                        other = second[state][arm][win][band][k]
                        for key,x in components.items():
                            absolute[key]=max(absolute[key],abs(x-other[key]))
    return absolute


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("native",type=Path,help="Directory containing FORK.json and WALLS.zip")
    args=ap.parse_args()
    repo=args.repository.resolve()
    method=pairs.published_method(repo)
    forks=json.loads((args.native/"FORK.json").read_bytes())
    results={"switches":pairs.WINDOWS, "proper_cut":pairs.CUT, "grids":{}}
    prior=method.base.pinned(repo/method.base.WALL,method.base.WALL_SHA)
    with zipfile.ZipFile(args.native/"WALLS.zip") as walls, zipfile.ZipFile(io.BytesIO(prior)) as old:
        for grid,fork in forks.items():
            q,weights,states=fork_states(fork)
            walls_grid={}
            independent_grid={}
            retained_checks={}
            worst=0.
            for geom in states:
                arms={}
                independent_arms={}
                for arm in ARMS:
                    with walls.open(f"{grid}/{geom}/{arm}/wall.csv") as f:
                        d=csv_data(f)
                    start,bg=pairs.wall_interpolants(d,fork["reference_A"])
                    if abs(start-fork["H0tau"])>2e-10 or abs(np.exp(bg[1](start))/fork["a"]-1)>1e-8:
                        raise ValueError("Wall starts on a different fork")
                    af=float(np.exp(bg[1](pairs.CUT)))
                    first=pairs.fourier_response(start,bg,q,1024,32)
                    other=pairs.direct_response(method,bg,q,16384,16)
                    cases={}
                    checks={}
                    selected = states if geom=="retain" else {k:states[k] for k in ("retain",geom)}
                    for k,(n,c) in selected.items():
                        cases[k]=pairs.evaluate(first,n,c,weights,af)
                        independent=pairs.evaluate(other,n,c,weights,af)
                        worst=max(worst,max_solver_difference(cases[k],independent))
                        checks[k]=independent
                    arms[arm]={"a_final":af,"states":cases}
                    independent_arms[arm]={"states":checks}
                    if geom=="retain":
                        retained_checks[arm]=check_source_equivalence(old,grid,d,fork["reference_A"],q,weights,cases["retain"])
                walls_grid[geom]=arms
                independent_grid[geom]=independent_arms
            midpoint=midpoint_test(walls_grid)
            midpoint_error=midpoint_mismatch(midpoint,midpoint_test(independent_grid))
            summary=contrast_summary(walls_grid)
            errors=contrast_errors(summary,contrast_summary(independent_grid))
            if max(errors.values())>2e-4:
                raise ValueError((grid,"Pulse contrast disagreement",errors))
            if worst>2e-6 or max(retained_checks.values())>1e-5:
                raise ValueError((grid,"solver comparison",worst,retained_checks))
            results["grids"][grid]={"walls":walls_grid,"contrasts":summary,
                "midpoint_fraction":midpoint,"independent_midpoint_error_fraction":midpoint_error,
                "independent_contrast_error_pp":errors,"max_independent_difference_over_total":worst,
                "original_vacuum_to_retained_fork_relative_discrepancy":retained_checks,
                "comoving_initial_number":{k:float(weights@n) for k,(n,c) in states.items()}}
    print(json.dumps(results,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
