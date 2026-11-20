"""Check the coupled energy, particle-rest budgets and sustained interval from the recorded state."""
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent


def check(actual, expected, tolerance, label):
    if not math.isfinite(actual) or abs(actual-expected)>tolerance:
        raise ValueError(f"{label}: {actual:.12g} versus {expected:.12g} (absolute tolerance {tolerance:g})")


def crossing(left,right,field,a,b):
    weight=-a/(b-a)
    return left[field]+weight*(right[field]-left[field])


def verify(data):
    if hashlib.sha256((HERE/"joined_forward.cpp").read_bytes()).hexdigest()!=data["solver_sha256"]:
        raise ValueError("Solver SHA mismatch")
    result={}
    for key in ["primary","historical","oldK_drive13"]:
        p=data["parameters"][key]
        s=data["runs"][key]
        rows=[dict(zip(data["trace_columns"],r)) for r in s["trace"]]
        first=rows[0]
        if any(first[k]!=0 for k in ("chi","chi_dot","D","S","U","gas","decay_energy","pairs","e","X","B","pi","psi","H")):
            raise ValueError("Nonempty initial wall excitations")
        check(first["mass"],4.,0.,"initial bulk mass")
        for r in rows:
            a3=r["a"]**3
            # Gamma cannot be accurately inferred from a speed rounded close to unity.
            # The conservative norm and EOS below avoid cancellation in this limit.
            energy=math.hypot(r["D"],r["S"])+r["U"]
            check(energy/a3,r["gas"],2e-10,"daughter energy")
            thermal=1.5*r["theta"]+2.25*r["theta"]**2/(1+math.hypot(1,1.5*r["theta"]))
            check(r["D"]*thermal/a3,r["intrinsic"],2e-10,"EOS intrinsic energy")
            check(r["D"]*p["charge_to_mass"]*r["velocity"]/a3,r["current"],
                  1e-7*max(1.,abs(r["current"])),"charge current")
            chi_kinetic=.5*r["chi_dot"]**2
            chi_potential=.5*(r["M"]*(r["chi"]-r["v"]))**2
            check(chi_kinetic+chi_potential,r["rho_chi"],2e-10,"parent energy")
            expected=.5*r["pi"]**2+r["V"]+r["rho_chi"]+r["sheet"]+r["gas"]+r["rad"]
            check(r["e"],expected,2e-10,"total energy")
            check(r["p_total"],.5*r["pi"]**2-r["V"]+chi_kinetic-chi_potential+(r["sheet"]+r["rad"])/3+r["pgas"],
                  2e-10,"total pressure")
            check(r["decay_power"],r["Gamma"]*r["chi_dot"]**2,2e-10,"local decay power")
            check(r["D"]+r["birth_kinetic"],r["decay_energy"],2e-10,"birth energy decomposition")
            check(r["D"],2*p["daughter_mass"]*r["pairs"],2e-10,"birth count")
        last=rows[-1];a3=last["a"]**3
        check(a3*last["gas"]+last["gas_pressure_work"],last["electric_work"]+last["decay_energy"],4e-6,"gas energy ledger")
        check(a3*last["rho_chi"]+last["chi_pressure_work"]+last["decay_energy"],last["drive_work"]+last["mass_work"],2e-6,"parent ledger")
        check(last["a"]**4*last["sheet"]+last["electric_work_a4"],last["comoving_wave_to_EM"]+last["comoving_scalar_to_EM"],3e-5,"Maxwell ledger")
        if abs(last["electric_work"])>p["charge_to_mass"]*last["field_exposure"]+1e-7:
            raise ValueError("Electric-current speed limit violated")
        result[key]={"peak_intrinsic_share":s["peak_hot_fraction"],"hot_efolds":s["hot_efolds"],
                     "electric_to_birth_work":s["peak_electric_work"]/s["peak_birth_energy"]}
        if key=="primary":
            def margin(r):return r["intrinsic"]-.9*r["e"]
            i=next(i for i,x in enumerate(rows) if x["H"]>0 and x["psi"]>=p["exit_field"]+p["exit_width"] and margin(x)>=0)
            j=next(j for j in range(i+1,len(rows)) if rows[j]["H"]<=0)
            if any(margin(r)<0 or r["V"]!=0 or r["rad"]!=0 for r in rows[i:j]):
                raise ValueError("Interrupted sampled hot period")
            ai=crossing(rows[i-1],rows[i],"a",margin(rows[i-1]),margin(rows[i]))
            af=crossing(rows[j-1],rows[j],"a",rows[j-1]["H"],rows[j]["H"])
            ef=math.log(af/ai)
            check(ef,s["hot_efolds"],.003,"independently sampled hot growth")
            if ef<=math.log(2):raise ValueError("No sampled hot doubling")
            target=2*s["first_hot_radius"]
            k=next(k for k in range(i+1,j) if rows[k]["a"]>=target)
            v0,v1=rows[k-1],rows[k]
            er=[crossing(v0,v1,f,v0["a"]-target,v1["a"]-target) for f in ("t","intrinsic","e")]
            check(er[0],s["doubling_time"],.002,"doubling time")
            check(er[1]/er[2],s["doubling_hot_fraction"],3e-5,"doubling share")
            result[key]["sampled_hot_efolds"]=ef
            result[key]["sampled_double_share"]=er[1]/er[2]
    for key in ("no_pulse","no_decay","no_drive"):
        check(data["runs"][key]["end_pair_count"],0.,0.,key)
    check(data["runs"]["no_charge"]["end_electric_work"],0.,0.,"uncharged")
    check(data["runs"]["reference18"]["hot_efolds"],data["runs"]["flipped_drive"]["hot_efolds"],1e-11,"drive symmetry")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data",nargs="?",type=Path,default=HERE/"JOINED.json")
    print(json.dumps(verify(json.loads(parser.parse_args().data.read_text())),indent=2))
