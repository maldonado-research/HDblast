"""Independent sign, resolution and extended-cutoff tests for the energy account."""
import argparse
import sys
import json
from pathlib import Path

import numpy as np


def methods():
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    import budget
    return budget


def collect(engine,repo,tables,work,jobs,rebuild=True):
    selected=(("fine","quiet"),("fine","minus"))
    if rebuild:
        engine.generate(repo,tables,jobs,cases=selected,models=("opposite",))
    polarization={}
    for grid,arm in selected:
        for scheme in ("bare","static"):
            first,second,both,opposite=(engine.load(tables,grid,arm,scheme,g)
                                       for g in ("first","second","both","opposite"))
            np.testing.assert_array_equal(both["s"],opposite["s"])
            relative={}
            for f in ("rho","p","j","rad","drho","dp","dj","drad"):
                err=opposite[f]-(2*(first[f]+second[f])-both[f])
                relative[f]=float(np.max(abs(err))/np.max(abs(opposite[f])))
            if max(relative.values())>5e-12:raise ValueError("Sign-flipped quadratic response")
            polarization[arm+"_"+scheme]=relative
    quiet=(("fine","quiet"),)
    extension={}
    refined={}
    for name,cut,steps in (("primary",120,4096),("refined120",120,8192),
                           ("uv240",240,4096),("uv480",480,4096),("refined480",480,8192)):
        directory=tables if name=="primary" else work/name
        if name!="primary" and rebuild:
            engine.generate(repo,directory,jobs,steps,quiet,cutoff=cut)
        node=engine.summarize(directory,quiet)["cases"]["fine_quiet"]
        by_scheme={}
        for scheme in ("bare","static"):
            added=node[scheme]["added"]["late"]
            by_scheme[scheme]={k:added[k] for k in ("change_energy","total_proper_change",
                "scalar_mass_work","gate_interaction_work","gate_completion_work","gate_work",
                "gate_work_path","expansion_work","energy_residual","max_pointwise_ward")}
            by_scheme[scheme]["initial_pulse_proper_daughter"]=node[scheme]["first"]["full"]["proper_energy"]["daughter"]
            by_scheme[scheme]["both_pulses_proper_daughter"]=node[scheme]["both"]["full"]["proper_energy"]["daughter"]
            if abs(added["energy_residual"])>1e-5:raise ValueError("Integrated late-work imbalance")
            if added["max_pointwise_ward"]>1e-8:raise ValueError("Pointwise late-work imbalance")
        if name.startswith("refined"):
            ref=extension[str(cut)]
            deltas={}
            for scheme,r in by_scheme.items():
                deltas[scheme]=dict(daughter_increment=r["change_energy"]["daughter"]-
                                    ref[scheme]["change_energy"]["daughter"],
                                    parent_increment=r["change_energy"]["parent"]-
                                    ref[scheme]["change_energy"]["parent"],
                                    gate_work=r["gate_work"]-ref[scheme]["gate_work"],
                                    total_proper_change=r["total_proper_change"]-
                                    ref[scheme]["total_proper_change"])
                if max(map(abs,deltas[scheme].values()))>2e-3:
                    raise ValueError("Step refinement of UV energy account")
            refined[str(cut)]=deltas
        else:
            extension[str(cut)]=by_scheme
    return dict(model="Quiet fixed wall with 64-node frequency panels, full proper energy coefficients per y^2",
                frequency_panels=[[0,60],[60,120],[120,240],[240,480]],
                first_panel_lower_per_mode="q",static_completion="C=sum(g^2/omega^2) for actual included oscillators",
                direct_opposite_sign_maximum_normalized_residuals=polarization,
                quiet_late_increments_per_cutoff=extension,
                quiet_difference_8192_minus_4096_subdivisions=refined)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo",type=Path,required=True)
    cli.add_argument("--tables",type=Path,required=True)
    cli.add_argument("--work",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    cli.add_argument("--jobs",type=int,default=2)
    cli.add_argument("--read-only",action="store_true")
    a=cli.parse_args()
    result=collect(methods(),a.repo,a.tables,a.work,a.jobs,not a.read_only)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
