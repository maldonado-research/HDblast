"""Integrate finite-time transverse mode and birth-noise response on supplied wall histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import struct
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parent
SCALES = (16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1280,
          1536, 1792, 2048, 2304, 2560, 3072, 4096, 6144, 8192, 12288, 16384, 24576, 32768)
PROBES = (32, 256, 2048, 8192)
TARGETS = {"dedicated_native":540.9609444060337,"dedicated_half":541.645547427229,
           "shared_native":673.7102609380573,"shared_half":674.0982875141613}


def sample(path, state, marker, cohorts, steps, groups, mode, action, mass, p_birth=0.):
    """Use midpoint birth times and exact comoving mass increments on each time bin."""
    time = state["t"]
    t0 = time[0]
    tf = marker["time"]
    if steps % (2*cohorts) or groups and cohorts % groups:
        raise ValueError("Birth and source intervals must fit integration grid")
    grid = np.linspace(t0, tf, 2*steps+1)
    background = np.column_stack([np.interp(grid,time,state[c]) for c in ("a","B","I")])
    bin_edges = np.linspace(t0,tf,cohorts+1)
    weights = np.diff(np.interp(bin_edges,time,state["D"]))
    if np.any(weights < 0) or not np.sum(weights)>0: raise ValueError("Invalid source increments")
    midpoints = ((np.arange(cohorts)*2+1)*steps)//(2*cohorts)
    group_ids = (np.arange(cohorts)*groups)//cohorts
    volume = (2*math.pi/mode)**3
    counts = np.bincount(group_ids,weights=weights)*action*volume/(6*mass)
    with path.open("wb") as out:
        out.write(struct.pack("=iiid",cohorts,steps,groups,(tf-t0)/steps))
        background.astype("=f8").tofile(out)
        for i,g in enumerate(group_ids):
            b = midpoints[i]
            # Optional axial birth velocity stress test, in units of the particle mass.
            ib = background[2*b,2]-background[2*b,0]*p_birth
            scale = math.sqrt(2/counts[g]) if groups and counts[g]>0 else 0.
            out.write(struct.pack("=dddii",weights[i],ib,scale,int(b),int(g)))
    return float(sum(counts))


def integrate(exe, path, mode, guide=1.):
    line=subprocess.check_output([str(exe),str(path),str(mode),str(guide)],text=True)
    out=np.array([float(v) for v in line.split()])
    if len(out)!=41 or not np.isfinite(out).all(): raise ValueError("Incomplete mode response")
    gram=out[9:25].reshape(4,4)
    fields=out[25:41].reshape(4,4)
    return {
        "unit_field_impulse_upper":out[0],"current_induced_impulse_upper":out[1],
        "magnetic_seed_impulse_upper":out[2],"unit_birth_modulation_impulse_upper":out[3],
        "poisson_rms_impulse_upper":out[4],"unit_birth_modulation_net_rms":out[5],
        "poisson_net_rms":out[6],"peak_perpendicular_field_operator_norm":out[7],
        "vacuum_phase":out[8],"net_rms_impulse_operator_norm":math.sqrt(max(np.linalg.eigvalsh(gram))),
        "terminal_field_operator_norm":np.linalg.svd(fields,compute_uv=False)[0],
        "terminal_field_transfer":fields.tolist()}


def configurations():
    for ch in ("dedicated","shared"):
        for mode in SCALES:
            yield ch,"half","scan",mode,512,4096,0,1.,0.
        for mode in PROBES:
            for label,n,steps,res in [("coarse",256,2048,"half"),("more_cohorts",1024,4096,"half"),
                    ("fine",1024,8192,"half"),("original_wall",512,4096,"native")]:
                yield ch,res,label,mode,n,steps,0,1.,0.
            yield ch,"half","no_guide",mode,512,4096,0,0.,0.
        for mode in (32,256,2048,8192,32768):
            for ngrp in (16,64,256):
                yield ch,"half","shot_"+str(ngrp),mode,256,2048,ngrp,1.,0.
        for mode in (2048,32768):
            yield ch,"half","shot_512",mode,512,4096,512,1.,0.
        for mode in (256,2048):
            for p in (-.32015621187164256,.32015621187164256):
                yield ch,"half",("negative" if p<0 else "positive")+"_birth",mode,512,4096,0,1.,p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history",type=Path,default=ROOT.parent/"history")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--pilot",action="store_true")
    args=parser.parse_args()
    meta=json.loads((args.history/"markers.json").read_text())
    states={}
    for name,run in meta["runs"].items():
        with np.load(args.history/(name+".npz")) as table: states[name]={k:table[k] for k in table.files}
    tasks=list(configurations())
    if args.pilot: tasks=[x for x in tasks if x[2]=="scan" and x[3] in PROBES]
    with tempfile.TemporaryDirectory(prefix="hd-transverse-modes-") as folder:
        temp=Path(folder);exe=temp/"modes"
        subprocess.run(["g++","-O3","-std=c++17",str(ROOT/"evolve.cpp"),"-o",str(exe)],check=True)
        def one(cfg):
            ch,mesh,name,mode,n,steps,groups,guide,p=cfg
            key=f"{ch}_{name}_{mode}";path=temp/key
            marker=meta["runs"][f"{ch}_{mesh}"]["markers"]["0.01"]
            pairs=sample(path,states[ch+"_"+mesh],marker,n,steps,groups,mode,meta["action"],meta["daughter_mass_L"],p)
            result=integrate(exe,path,mode,guide)
            target=TARGETS[ch+"_"+mesh]
            result.update(channel=ch,wall_mesh=mesh,mode=mode,cohorts=n,time_steps=steps,
                          source_bins=groups,guide_scale=guide,axial_birth_p_over_m=p,
                          pairs_per_wavelength_cube=pairs,target_rms_impulse=target,
                          minimum_seed_rms=target/result["unit_field_impulse_upper"],
                          minimum_seed_rms_induced_only=target/result["current_induced_impulse_upper"],
                          minimum_fractional_birth_modulation=target/result["unit_birth_modulation_impulse_upper"])
            path.unlink()
            return key,result
        with ThreadPoolExecutor(args.jobs) as pool: results=dict(pool.map(one,tasks))
    args.output.write_text(json.dumps({"seed_order":["B_y","E_z","B_z","E_y"],"runs":results},
                                      indent=2,sort_keys=True,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
