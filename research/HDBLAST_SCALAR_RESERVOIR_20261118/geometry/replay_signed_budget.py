#!/usr/bin/env python3
"""Replay two signed local budgets using only the accompanying derived CSVs.

No raw histories, old producers, or third-party packages are required.
The boundary flux is inferred using the stated junction, not measured independently.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
HEADER=('x','tension_energy','tension_expansion','tension_work',
        'matter_energy','matter_expansion','matter_work')


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def trapezoid(x,y):
    return math.fsum((x[i+1]-x[i])*(y[i]+y[i+1])/2 for i in range(len(x)-1))


def evaluate(path):
    with path.open(newline='') as f:
        reader=csv.DictReader(f); require(tuple(reader.fieldnames)==HEADER,'unexpected columns')
        rows=[{k:float(v) for k,v in row.items()} for row in reader]
    require(len(rows)>=3,'insufficient rows')
    require(all(math.isfinite(v) for r in rows for v in r.values()),'nonfinite value')
    x=[r['x'] for r in rows]
    require(x[0]==1 and x[-1]==4 and all(a<b for a,b in zip(x,x[1:])),'wrong window/clock')
    values={k:[r[k] for r in rows] for k in HEADER[1:]}
    de_t=values['tension_energy'][-1]-values['tension_energy'][0]
    de_m=values['matter_energy'][-1]-values['matter_energy'][0]
    expansion_t=trapezoid(x,values['tension_expansion'])
    expansion_m=trapezoid(x,values['matter_expansion'])
    force_t=trapezoid(x,values['tension_work'])
    force_m=trapezoid(x,values['matter_work'])
    tension_by_parts=de_t-expansion_t
    matter_by_parts=de_m+expansion_m
    net_force=force_t+force_m
    net_by_parts=tension_by_parts+matter_by_parts
    require(force_m>0 and force_t<0 and net_force<0 and net_by_parts<0,'signed result changed')
    require(abs(tension_by_parts-force_t)/abs(force_t)<1e-3,'tension forms fail finite tolerance')
    require(abs(matter_by_parts-force_m)/abs(force_m)<1e-3,'matter forms fail finite tolerance')
    every_other=list(range(0,len(x),2))
    if every_other[-1]!=len(x)-1:every_other.append(len(x)-1)
    xt=[x[i] for i in every_other]
    thin_force=trapezoid(xt,[values['tension_work'][i]+values['matter_work'][i] for i in every_other])
    return dict(csv=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=len(rows),
        tension_energy_start=values['tension_energy'][0],tension_energy_end=values['tension_energy'][-1],
        matter_energy_start=values['matter_energy'][0],matter_energy_end=values['matter_energy'][-1],
        delta_tension_energy=de_t,delta_matter_energy=de_m,
        positive_tension_expansion_term=expansion_t,matter_pressure_work=expansion_m,
        tension_work_from_scalar_force=force_t,tension_work_by_parts=tension_by_parts,
        matter_work_from_scalar_force=force_m,matter_work_by_parts=matter_by_parts,
        tension_forms_relative_difference=(tension_by_parts-force_t)/abs(force_t),
        matter_forms_relative_difference=(matter_by_parts-force_m)/abs(force_m),
        inferred_two_copy_bulk_to_shell_flux=net_force,
        total_shell_energy_plus_pressure=net_by_parts,
        net_forms_relative_difference=(net_by_parts-net_force)/abs(net_force),
        magnitude_tension_work_over_matter_work=abs(force_t)/force_m,
        inferred_flux_after_thinning=thin_force,
        omitted_tension_wrong_inferred_flux=force_m,
        missing_one_copy_wrong_inferred_flux=net_force/2,
        reversed_orientation_wrong_inferred_flux=-net_force)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()
    report=dict(scope='Finite trapezoid comparison of local shell budgets. Inferred flux is not an independent bulk-reservoir measurement.',
                units='Code energies multiplied by kappa5^2, per common a(0)=1 comoving volume; x=H0*tau.',
                endpoint_rule='Separate linear interpolation of each weighted saved column at x=1,4; all interior stored nodes retained.',
                runs=[evaluate(HERE/'data'/('b1_'+g+'.csv')) for g in ('coarse','fine')])
    if a.output: a.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
