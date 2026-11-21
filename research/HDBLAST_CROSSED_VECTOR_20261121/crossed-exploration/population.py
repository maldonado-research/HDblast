"""All-pair equal-time CM angular transport for the local massive-vector completion of the isolated pole."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
PHYS=runpy.run_path(str(HERE/"check.py"))
M=PHYS["M"];G2=PHYS["G2"]
NAMES=("pole_total","pole_quadrupole","unlike_total","unlike_quadrupole","unlike_transfer","like_total","like_quadrupole")
CUT=140379347.50895345  # Exact finite-mass |Re a0|=1/2 threshold from the helicity check.


def areas(s,nodes=48):
    """Exact angular kernels projected by quadrature, with a coordinate resolving both forward poles."""
    s=np.asarray(s).reshape(-1,1)
    x,w=np.polynomial.legendre.leggauss(nodes)
    k=(s-4.)/(2*M**2);width=np.log1p(k)
    # The integration variable is log(1+K*(1-|z|)). Use the continuous threshold limit.
    b=(1.+x)/2.
    h=np.divide(np.expm1(width*b),k,out=np.broadcast_to(b,(len(s),len(b))).copy(),where=k>0.)
    dh=np.divide(width,k,out=np.ones_like(k),where=k>0.)*(1.+k*h)*w[None,:]/2.
    z=1.-h
    values=np.zeros((len(s),len(NAMES)))
    for sign in (1.,-1.):
        a,b,c,d,e,f=PHYS["terms"](s,sign*z)
        unlike=a+b+c;like=(d+e+f)/2.
        quadrupole=1.5*h*(2.-h);trans=1.-sign*z
        for j,col in enumerate((a,a*quadrupole,unlike,unlike*quadrupole,unlike*trans,like,like*quadrupole)):
            values[:,j]+=(col*dh).sum(axis=1)
    return values*(G2**2/(32.*np.pi*s))


def table(exact_limit,power,directory):
    n=int(np.ceil(np.log(exact_limit/4.)*power))+3
    grid=4.*np.exp(np.arange(n)/power)
    out=np.concatenate([areas(grid[i:i+256]) for i in range(0,n,256)])
    path=directory/("sigma-"+str(power))
    out.tofile(path)
    # Compare floats at the table midpoints (every 13th bin) to fresh high-order integrals.
    ix=np.arange(0,len(grid)-1,13)
    center=np.sqrt(grid[ix]*grid[ix+1])
    check=np.concatenate([areas(center[i:i+128],96) for i in range(0,len(ix),128)])
    error=np.max(abs((out[ix]+out[ix+1])/(2*check)-1.),axis=0)
    if np.max(error)>6/(power**2):raise ValueError(("Table accuracy",error,power))
    return path,error.tolist()


def measurements(total,factor,fraction):
    r,h,area,fast,hfast,slow=total
    if r==0:return None
    return dict(rate_per_particle_wall_time=float(factor*r),hot_participation_share=float(h/(2*r)),
        hot_rate_over_all_particle_rate=float(h/(2*r*fraction)),cross_section_weighted_flux=float(r/area),
        from_pair_cm_boost_over_10=float(fast/r),hot_participation_from_boost_over_10=float(hfast/h) if h else 0.,
        from_relative_flux_below_point_one=float(slow/r))


def one(params):
    file,power,folder,binary,tab,cut=params
    with np.load(file) as d:
        p=d["p_final"];eg=np.hypot(1.,np.linalg.norm(p,axis=1));N=int(d["pairs"]);time=float(d["time"])
        if len(p)!=2*N or float(d["channel_fraction"])!=1.:raise ValueError("Wrong population")
        hot=np.zeros(len(p));n_hot=(len(p)+9)//10;hot[np.argsort(-eg,kind="stable")[:n_hot]]=1.
        factor=float(d["action"])*float(d["D"])/(6.*float(d["daughter_mass_L"])**3*float(d["endpoint_a"])**3*N*N)
        inp=folder/(file.stem+"-"+str(power)+".bin");np.column_stack((p,hot)).tofile(inp)
    out=inp.with_suffix(".txt")
    subprocess.run([str(binary),str(inp),str(tab),str(power),str(cut),str(out)],check=True)
    blocks=[np.fromstring(row,sep=" ") for row in out.read_text().splitlines()]
    if len(blocks)!=16 or not np.array_equal(blocks[0]+blocks[8],[N*N,N*(N-1)]):raise ValueError("Pair count")
    low=np.array(blocks[1:8]);hi=np.array(blocks[9:]);frac=n_hot/(2*N)
    def section(x):
        values={name:measurements(row,factor,frac) for name,row in zip(NAMES,x)}
        for i,name in ((2,"complete_total"),(3,"complete_quadrupole")):
            values[name]=measurements(x[i]+x[i+3],factor,frac)
        return values
    full=low+hi
    def masspart(x):return float(x[3,0]+x[6,0]-x[1,0])
    return dict(tape=file.name,time=time,pairs_per_sign=N,top_energy_share=float(hot@eg/eg.sum()),
        top_particles=n_hot,bins_per_log_s=power,cut_s_over_M2=cut/(M*M),unordered_below=blocks[0].astype(int).tolist(),
        unordered_above=blocks[8].astype(int).tolist(),all_pairs=section(full),below_limit=section(low),above_limit=section(hi),
        enhancement_fraction_above_limit=masspart(hi)/masspart(full),combined_over_pole_quadrupole=float((full[3,0]+full[6,0])/full[1,0]),
        below_combined_over_pole_quadrupole=float((low[3,0]+low[6,0])/low[1,0]),
        combined_fraction_above_limit=float((hi[3,0]+hi[6,0])/(full[3,0]+full[6,0])),
        completion_below_plus_isolated_above_over_all_isolated=float(1.+masspart(low)/full[1,0]))


def run(args):
    files=sorted(args.artifacts.glob("*.npz"))
    if len(files)!=8:raise ValueError("Expected eight checked endpoint tapes")
    maximum=0.
    for path in files:
        with np.load(path) as d:
            p=d["p_final"]
            maximum=max(maximum,4*(1.+np.max((p*p).sum(axis=1))))
    with tempfile.TemporaryDirectory() as temp:
        work=Path(temp);binary=work/"census"
        subprocess.run(["g++","-O3","-ffp-contract=off","-std=c++17",str(HERE/"population.cpp"),"-o",str(binary)],check=True)
        tables={};error={}
        for n in (256,512):
            tables[n],error[str(n)]=table(maximum,n,work)
        tasks=[(f,512,work,binary,tables[512],CUT) for f in files]
        tasks.extend((f,256,work,binary,tables[256],CUT) for f in files if "p14" in f.stem)
        with ThreadPoolExecutor(args.jobs) as tp:data=list(tp.map(one,tasks))
    if args.isolated_reference:
        refs={r["snapshot"]:r for r in json.loads(args.isolated_reference.read_bytes())}
        for row in data:
            r=refs[Path(row["tape"]).stem]
            if not np.isclose(row["all_pairs"]["pole_total"]["hot_participation_share"],r["pole_top_participation_share"],rtol=2e-5):
                raise ValueError("Previously measured all-pair isolated pole cohort differs")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(summary=data,kernel_midpoint_checks=error),sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts",required=True,type=Path);parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--isolated-reference",type=Path,help="Optional previous pairs/RESULTS.json for an independent baseline check")
    parser.add_argument("--jobs",type=int,default=4)
    run(parser.parse_args())
