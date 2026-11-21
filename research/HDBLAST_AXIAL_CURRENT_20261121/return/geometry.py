"""Record the charged background collar geometry and a manufactured wall-driven k=2 wave."""
import argparse
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np


def instrument(study, temp):
    source = study/"parent"
    prior = source/"parent/parent/prior"
    build = runpy.run_path(str(source/"geometry/probe.py"))["build"]
    build(prior, temp, (2,))
    cpp = temp/"probe.cpp"
    text = cpp.read_text()
    changes = [
        ("struct Solver {", """static std::ofstream node_stream;
static double wave_time=0.;
static bool capture_geometry=false;
static double test_drive(double t){
    if(t<=6.252) return 0.;
    double x=(t-6.252)/.0015;
    return -std::expm1(-x*x)*std::cos(35.*(t-6.252));
}
struct Solver {"""),
        ("(m==unforced ? 0. : emitted)", "0."),
        ("y[PROBE_MOM+m]/gauge-inc", "(m==unforced ? 0. : test_drive(wave_time))-inc"),
        ("            probe_cell(y,dy,cell,geom,emitted,w.k);",
         """            if(capture_geometry){
                for(const auto &z:geom){
                    const double record[]={z.r,z.rx,z.u,z.n,z.nx};
                    node_stream.write(reinterpret_cast<const char*>(record),sizeof(record));
                }
            }
            probe_cell(y,dy,cell,geom,emitted,w.k);"""),
        ("double last_probe_print=-1.;", """double last_probe_print=-1., previous_geometry=-1.;
    node_stream.open(std::string(argv[35])+".nodes",std::ios::binary);
    node_stream.exceptions(std::ios::badbit|std::ios::failbit);"""),
        ("        s.rhs(y, dy, check ? &z : nullptr);", """        wave_time=t;
        capture_geometry=t>=6.248 && (t-previous_geometry>=.00012 || t+2.*dt0>=end);
        if(capture_geometry){
            const double head[]={t,y[s.index(0,3)],y[s.index(0,4)],test_drive(t)};
            node_stream.write(reinterpret_cast<const char*>(head),sizeof(head));
            previous_geometry=t;
        }
        s.rhs(y, dy, check ? &z : nullptr);
        capture_geometry=false;"""),
    ]
    for old,new in changes:
        if text.count(old)!=1:
            raise ValueError("Expected inherited source at "+old)
        text=text.replace(old,new)
    stages="        s.rhs(stage, dy);"
    if text.count(stages)!=3:
        raise ValueError("Unknown time integrator")
    parts=text.split(stages)
    text=parts[0]
    for delta,part in zip(("dt/2","dt/2","dt"),parts[1:]):
        text+="        wave_time=t+"+delta+";\n"+stages+part
    cpp.write_text(text)
    exe=temp/"radial"
    subprocess.run(["g++","-std=c++17","-O3","-ffp-contract=off","-I",str(prior/"boost"),str(cpp),"-o",str(exe)],check=True)
    return exe,prior


def export(study, out, end, shared=False):
    out.mkdir(parents=True,exist_ok=True)
    metadata=json.loads((study/"parent/parent/history/markers.json").read_text())
    with tempfile.TemporaryDirectory(prefix="collar-geometry-") as folder:
        tmp=Path(folder)
        exe,prior=instrument(study,tmp)
        config=json.loads((prior/"prior/envelope/ENVELOPE.json").read_bytes())["parameters"]["mass_4e12"]
        for name in (["dedicated","shared"] if shared else ["dedicated"]):
            trace=tmp/"wall"
            run=tmp/"run"
            args=[2,metadata["action"],prior/"boost/BOOST_TABLE.dat",1,metadata["runs"][name+"_half"]["channel_fraction"],trace]
            run.write_text("#!/bin/sh\nexec "+shlex.quote(str(exe))+' "$@" '+shlex.join(list(map(str,args)))+"\n")
            run.chmod(0o755)
            config_run=dict(config,cfl=.5*config["cfl"],end=end+.0005,sample=0)
            report=runpy.run_path(str(prior/"prior/coupled/joined.py"))["run"](run,config_run)
            if abs(report["end_tau"]-config_run["end"])>1e-9:
                raise ValueError("Background failed to reach required endpoint")
            cells=config["late_cells"]
            rows=np.fromfile(str(trace)+".nodes").reshape(-1,4+cells*4*5)
            t,reflection,zero,drive=rows[:,:4].T
            if np.any(np.diff(t)<=0.) or t[0]>=6.25 or t[-1]<end or np.any(zero!=0.):
                raise ValueError("Bad geometry coverage or zero-wave control")
            geo=rows[:,4:].reshape(-1,cells,4,5)
            if np.any(geo[:,:,:,0]<=0) or np.any(geo[:,:,:,3]<=0):
                raise ValueError("Lapse or radius not positive")
            np.savez_compressed(out/(name+".npz"),t=t,geo=geo,reference=reflection,drive=drive,depth=config["depth"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polarization",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--end",type=float,default=6.67)
    parser.add_argument("--shared",action="store_true")
    options=parser.parse_args()
    export(options.polarization.resolve(),options.output_dir,options.end,options.shared)
