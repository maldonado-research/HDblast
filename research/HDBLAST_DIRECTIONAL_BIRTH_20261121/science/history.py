"""Regenerate read-only native-step histories from the archived #192 equations."""
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
from concurrent.futures import ThreadPoolExecutor

CASES = ("mass_4e10", "mass_4e12")
STEP_COLUMNS = ("t", "dt", "rest0", "rest1", "k0", "k1", "electric0", "electric1", "center0", "center1")
WINDOW = (.00001, .99999)


def prepare(prior, scratch, name):
    env = json.loads((prior/"envelope/ENVELOPE.json").read_bytes())
    params = env["parameters"][name]
    centered = params["daughter_mass"] >= 1.e11
    transform = runpy.run_path(str(prior/"envelope/build.py"))["transform"]
    src = transform((prior/"coupled/joined_forward.cpp").read_text(), centered=centered)
    expected = env["centered_solver_sha256" if centered else "transformed_solver_sha256"]
    if hashlib.sha256(src.encode()).hexdigest() != expected:
        raise ValueError("Archived solver fingerprint differs")
    anchor = "        const double margin=gas.intrinsic_random-.9*wnew.e;"
    opening = "int main(int argc, char **argv) {"
    closing = '    out["end_tau"] = t;'
    if any(src.count(mark)!=1 for mark in (anchor, opening, closing)):
        raise ValueError("Unexpected source layout")
    rest = env["runs"][name]["end_rest"]
    # Separate output paths for the two native time resolutions.
    folder = scratch/name
    folder.mkdir(exist_ok=True, parents=True)
    prologue = r'''
    std::ofstream stream(std::string(argv[0]) + "." + argv[7] + ".bin",std::ios::binary);
    if (!stream) throw std::runtime_error("Unable to open native readout");
'''
    code = r'''
        if(gas.comoving_rest > RESTLOW && previous_gas.comoving_rest < RESTHIGH &&
           gas.comoving_rest>previous_gas.comoving_rest) {
            const auto shell = [&](double radius, double psi) {
                const double ratio=mass_low+(mass_high-mass_low)*step((psi-threshold_phi)/field_width+.5)[0];
                return daughter_mass*radius*std::sqrt(std::max(0.,(ratio-1.)*(ratio+1.)));
            };
            const double scale=daughter_mass*charge_mass;
            const double row[]={t, dt, previous_gas.comoving_rest, gas.comoving_rest,
                shell(olda,oldphi), shell(a,y[PHI]), scale*olda*previous_wall.x,scale*a*wnew.x,
                .5*scale*olda*olda*previous_wall.b, .5*scale*a*a*wnew.b};
            stream.write(reinterpret_cast<const char*>(row),sizeof(row));
        }
'''
    code = code.replace("RESTLOW", repr(WINDOW[0]*rest)).replace("RESTHIGH",repr(WINDOW[1]*rest))
    src = "#include <fstream>\n"+src.replace(opening,opening+prologue).replace(anchor,code+anchor)
    src = src.replace(closing,'    stream.close();\n    if (!stream) throw std::runtime_error("Unable to write native readout");\n'+closing)
    cpp, binary = folder/"forward.cpp", folder/"forward"
    cpp.write_text(src)
    subprocess.run([os.environ.get("CXX","g++"),"-std=c++17","-O3",str(cpp),"-o",str(binary)],check=True)
    return params,env["runs"][name],env["screen"][name],binary


def build_histories(prior, scratch, jobs):
    runner = runpy.run_path(str(prior/"coupled/joined.py"))["run"]
    builds = {}
    for name in CASES:
        builds[name] = prepare(prior, scratch, name)
    def evolve(spec):
        name, label, multiplier = spec
        params, ref, _, exe = builds[name]
        cfl = params["cfl"]*multiplier
        p = dict(params, sample=0, cfl=cfl)
        vals = runner(exe,p)
        delta = abs(vals["end_rest"]/ref["end_rest"]-1)
        if delta>(3.e-8 if label=="native" else 3.e-4):
            raise ValueError("Unreproduced reference rest production")
        return name,label,dict(path=str(exe)+"."+str(cfl)+".bin",end_rest=vals['end_rest'],
                                hot_factor=math.exp(vals['hot_efolds']))
    specs = [(name,l,m) for name in CASES for l,m in (("native",1.),("half_step",.5))]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows=list(pool.map(evolve,specs))
    result = {}
    for name,label,fields in rows:
        p,r,s,_=builds[name]
        if name not in result:
            result[name]=dict(mass=p['daughter_mass'],q=p['charge_to_mass'],reference_rest=r['end_rest'],
                              action=s['witness_action'] if name==CASES[0] else 5*s['elastic_drag_action_floor'],
                              floor=s['elastic_drag_action_floor'],runs={})
        result[name]['runs'][label]=fields
    return result
