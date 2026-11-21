"""Resolve continuous injection immediately before the reference 1%-work marker."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

BUDGET = runpy.run_path(str(Path(__file__).with_name("budget.py")))


def solver(prior, scratch):
    fb = runpy.run_path(str(prior / "feedback/feedback.py"))
    old = prior / "prior"
    source = (old / "coupled/joined_forward.cpp").read_text()
    source = runpy.run_path(str(old / "envelope/build.py"))["transform"](source, centered=True)
    source = runpy.run_path(str(prior / "boost/evolve.py"))["make_source"](source, fb, (prior / "feedback/angular.hpp").read_text())
    change = fb['change']
    source = '#include <fstream>\n' + source
    source = change(source, "if (argc != 34)", "if (argc != 39)")
    source = change(source, "int collision_mode = 0;", "double branching=1.;\nint collision_mode = 0;")
    source = change(source, "collision_mode=std::stoi(argv[33]);", """collision_mode=std::stoi(argv[33]);
    branching=std::stod(argv[34]);
    const double mark=std::stod(argv[35]), min_focus_step=std::stod(argv[36]), width_focus=std::stod(argv[37]);
    std::ofstream trace(argv[38],std::ios::binary);
    trace.exceptions(std::ios::badbit|std::ios::failbit);""")
    source = change(source, "return proper_density*boosted_rate(", "return branching*proper_density*boosted_rate(")
    source = change(source, "        bool due = want_samples", """        // Resolve the terminal window, preserving the inherited pre-window evolution.
        if (t >= mark-width_focus && t<mark)
            dt=std::min({dt,std::max(min_focus_step,(mark-t)/64.),mark-t});
        bool due = want_samples""")
    anchor = "        const double margin=gas.intrinsic_random-.9*wnew.e;"
    source = change(source, anchor, """        if (t+dt>=mark-width_focus) {
            const double row[]={t+dt,a,gas.comoving_rest,gas.theta,oscillator(y).mass_ratio};
            trace.write(reinterpret_cast<const char*>(row),sizeof(row));
        }
""" + anchor)
    cpp = scratch / 'focus.cpp'
    cpp.write_text(source)
    out = scratch / 'focus'
    subprocess.run(['g++','-std=c++17','-O3','-I',str(prior/'boost'),str(cpp),'-o',str(out)],check=True)
    return out


def evaluate(job, exe, folder, prior, meta, config, driver, history):
    name, delta, step = job
    spec = meta['runs'][name]
    ref_time = spec['markers']['0.01']['time']
    mark = ref_time + delta
    tag = f'{name}_{delta:g}_{step:g}'
    path = folder / (tag + '.bin')
    wrapper = folder / tag
    width = 1.5e-3
    args = ['2', repr(meta['action']), str(prior / 'boost/BOOST_TABLE.dat'), '1',
            repr(spec['channel_fraction']), repr(mark), repr(step), repr(width), str(path)]
    wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(args)+'\n')
    wrapper.chmod(0o755)
    try:
        driver(wrapper, dict(config, cfl=config['cfl']*(.5 if name.endswith('half') else 1.), end=mark+5.e-8, sample=0))
    except subprocess.CalledProcessError as err:
        raise RuntimeError(err.stderr) from err
    raw = np.fromfile(path, dtype='=f8').reshape(-1, 5)
    data = {key: raw[:, i] for i, key in enumerate(('t','a','D','theta','R'))}
    # These backgrounds start after many particles have formed. Account only for new ones.
    t = data['t']
    j = np.searchsorted(t, mark)
    if t[j] != mark:
        raise ValueError("Focused replay did not land on requested time")
    rest = data['D']
    factor = spec['channel_fraction'] * meta['action'] * BUDGET['pair_kernel_upper']() / (6. * meta['daughter_mass_L']**3)
    start, age, mass, theta, right = BUDGET['enclosed_suffix'](data, mark, factor)
    want = BUDGET['inverse_cdf'](mass, theta)
    force = np.maximum(0., want-BUDGET['TUBE'])/age
    best = int(force.argmax())
    field = force[best] / meta['charge_to_mass']
    slope_left = np.searchsorted(t, mark - 2.e-8)
    slope = (rest[j]-rest[slope_left])/(t[j]-t[slope_left])/rest[j]
    with np.load(history / (name + '.npz'), allow_pickle=False) as old:
        end_original = int(np.searchsorted(old['t'], mark))
        peak_x = float(np.max(np.abs(old['X'][:end_original+1])))
    quantile_bounds = {}
    for tol in (.0001, .001, .01):
        lower = np.maximum(0., BUDGET['inverse_cdf'](np.maximum(0.,mass-tol),theta)-BUDGET['TUBE'])/age
        where = int(np.argmax(lower))
        quantile_bounds[str(tol)] = dict(equivalent_electric_field_lower=float(lower[where]/meta['charge_to_mass']),
            witness_age_L=mark-start[where])
    gap_at_bounded_force = {}
    for relative in (.1,.25,1.,2.,3.):
        effective = relative * peak_x
        limit = BUDGET['TUBE'] + (effective * meta['charge_to_mass']) * age
        mismatch = mass - BUDGET['thermal_cdf'](limit,theta)
        where = int(np.argmax(mismatch))
        min_birth_cap = float(np.max(want - (effective * meta['charge_to_mass']) * age))
        gap_at_bounded_force[str(relative)] = dict(direction_sum_cap=effective, tv_floor=max(0.,mismatch[where]),
            min_newborn_transverse_cap_over_mc=max(0.,min_birth_cap),
            min_mass_ratio_for_cold_two_body_parent=math.hypot(1.,max(0.,min_birth_cap)),
            witness_age_L=mark-start[where], minimum_recent_unhit_fraction=mass[where], momentum_cutoff=limit[where])
    theta0 = data['theta'][j]
    # Hot-gas, steady newborn-rate asymptote with a uniform upper bound on birth p_perp.
    estimate = theta0**2 * slope / (BUDGET['TUBE']*meta['charge_to_mass'])
    return tag, dict(run=name, offset=delta, marker_time=mark, minimum_step=step, num_native_cells=len(t)-1,
        terminal_fractional_birth_rate=slope, theta=theta0, rest_at_mark=rest[j],
        equivalent_electric_field_lower=field, peak_force_lower_L_over_mc=force[best], peak_parallel_field=peak_x,
        force_at_cdf_tolerance=quantile_bounds, tv_for_force_in_units_of_peak_parallel_electric=gap_at_bounded_force,
        lookback_L_at_witness=mark-start[best], impulse_time_upper=age[best],
        thermal_mass_at_cutoff=BUDGET['thermal_cdf'](BUDGET['TUBE'] + age[best]*force[best],theta),
        guaranteed_recent_uncollided_share=mass[best], transverse_cutoff_over_mc=BUDGET['TUBE']+age[best]*force[best],
        ultrarelativistic_constant_influx_estimate=estimate)


def main():
    pa = argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--prior', type=Path, required=True)
    pa.add_argument('--history', type=Path, required=True)
    pa.add_argument('--output', type=Path, required=True)
    pa.add_argument('--jobs', type=int, default=2)
    args = pa.parse_args()
    prior = args.prior.resolve()
    meta = json.loads((args.history / 'markers.json').read_text())
    config = json.loads((prior / 'prior/envelope/ENVELOPE.json').read_text())['parameters']['mass_4e12']
    driver = runpy.run_path(str(prior / 'prior/coupled/joined.py'))['run']
    jobs = [(name,0.,step) for name in meta['runs'] for step in (1.e-8,1.e-9,1.e-10,2.e-11)]
    jobs.extend((name,delta,2.e-11) for name in meta['runs'] if name.endswith('half') for delta in (-2.e-5,2.e-5))
    with tempfile.TemporaryDirectory(prefix='hdblast-focused-injection-') as parent:
        scratch = Path(parent)
        exe = solver(prior, scratch)
        with ThreadPoolExecutor(args.jobs) as pool:
            answers = dict(pool.map(lambda j: evaluate(j,exe,scratch,prior,meta,config,driver,args.history), jobs))
    args.output.write_text(json.dumps({'near_one_percent':answers},sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
