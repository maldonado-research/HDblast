"""Verify or recompute the HDblast initial-source-matching comparisons."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
PRIOR = Path("research/HDBLAST_MOVING_PARENT_20261119")
# Exact manifest of the moving-parent archive at HDblast revision 930e1b76.
PRIOR_MANIFEST = "a2281d51980e0e1fb9ca60155a516991ee3baaef6bd36114bfcd54a567e819ad"
DEPENDENCIES = ("SOURCES.json", "phase_feasibility/feasibility.py", "coherence/coherence.py",
                "independent/transfer.py", "independent/squeeze.py")
RUNS = {"baseline":("source_match/BASELINE.json",[]),
        "refined":("source_match/REFINED.json",["--nodes","48","--samples","2048","--padding","32"]),
        "high":("source_match/HIGH_Q.json",["--nodes","96","--samples","2048","--padding","32"])}
ARMS = ("plus","minus","far")
ROOT2 = math.sqrt(2.)


def require(condition, label):
    if not condition:
        raise ValueError(label)


def close(a,b,label,rel=1e-7,atol=1e-11):
    require(math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=rel,abs_tol=atol),
            f"{label}: {a!r} vs {b!r}")


def checked(path,digest):
    data=path.read_bytes()
    require(hashlib.sha256(data).hexdigest()==digest, f"SHA-256 mismatch for {path}")
    return data


def authenticate(repo):
    archive=repo/PRIOR
    manifest=checked(archive/"SHA256SUMS",PRIOR_MANIFEST).decode()
    pinned={name:digest for digest,name in (line.split(maxsplit=1) for line in manifest.splitlines())}
    for name in DEPENDENCIES:
        checked(archive/name,pinned[name])
    sources=json.loads((archive/"SOURCES.json").read_bytes())["sources"]
    for item in sources:
        path=(repo/item["path"]).resolve()
        require(path.is_relative_to(repo),"Historical path outside repository")
        data=checked(path,item["sha256"])
        require(len(data)==item["bytes"],f"Historical byte count for {path}")


def check_independent(report):
    require(report.keys()=={"32","64","128"},"Independent refinements")
    for gridsize, control in report.items():
        require(control.keys()==set(ARMS),"Independent arms")
        for arm,v in control.items():
            label=f"Independent {gridsize}/{arm}"
            s=v["separable_bound"]
            q=v["quantum_bound"]
            close(ROOT2*s,q,label+" sqrt(2)")
            close(s,v["separable_witness"],label+" primal/dual",atol=1e-7)
            close(q,v["quantum_witness"],label+" quantum primal/dual",atol=1e-7)
            ratio=abs(v["diagonal_excess_contrast"]/v["diagonal_absolute_contrast"])
            close(s*ratio,v["separable_absolute_bound"],"Absolute support normalization")
            close(q*ratio,v["quantum_absolute_bound"],"Absolute quantum normalization")
            close(ROOT2*v["imag_separable"],v["imag_quantum"],"Imaginary support")
            require(v["imag_quantum"]<1 and v["imag_separable"]<=s<=v["unconstrained_separable"],
                    label+" nested supports")
            require(s<1 and v["quantum_absolute_bound"]<1 and (q>1)==(arm!="far"), label+" threshold signs")
            require(v["max_moment_residual"] < 1e-13 and v["max_factor_change_12_to_24"] < 1e-8,
                    label+" source and integration checks")
    for arm in ARMS:
        for key in ("separable_bound","quantum_bound","quantum_absolute_bound","imag_quantum"):
            close(report["64"][arm][key],report["128"][arm][key],f"Independent refinement {arm}/{key}",
                  rel=0.,atol=1e-7)


def run(script,repo,args):
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    print("Rebuilding",script,*args,file=sys.stderr,flush=True)
    result=subprocess.run([sys.executable,"-I","-B",str(HERE/script),str(repo),*args],
                          env=env,stdout=subprocess.PIPE,check=True)
    return json.loads(result.stdout)


def compare_direct(a,b):
    for grid,controls in b.items():
        for arm,y in controls.items():
            x=a[grid][arm]
            for k in ("diagonal_excess_contrast","diagonal_absolute_contrast","separable_bound",
                      "quantum_bound","imag_separable","imag_quantum","unconstrained_separable",
                      "separable_absolute_bound","quantum_absolute_bound"):
                close(x[k],y[k],f"Transfer replay {grid}/{arm}/{k}",rel=2e-6,atol=1e-10)


def quadrature_weights(q):
    """Recover the Legendre weights independently from the saved momentum nodes."""
    n=len(q)
    weights=[]
    for k in q:
        x=k/4.-1
        previous,current=1.,x
        for j in range(2,n+1):
            previous,current=current,((2*j-1)*x*current-(j-1)*previous)/j
        require(abs(current)<5e-12, "Off-grid Legendre node")
        derivative=n*(previous-x*current)/(1-x*x)
        weights.append(4*2*k*k/((1-x*x)*derivative**2))
    close(sum(weights),8**3/3,"Integrated parent measure")
    return weights


def certificates(group):
    for name in ("fine_moments","imaginary_only","both_mesh_moments"):
        if name in group:
            yield name,group[name]
    for arm,report in group.get("individual",{}).items():
        yield arm,report


def check_primary(data):
    labels=data["labels"]
    expected={f"{mesh}/{win}/{band}/{arm}" for mesh in ("fine","coarse") for win in ("wide","narrow")
              for band in ("number","screened_number") for arm in ARMS}
    require(len(labels)==len(expected) and set(labels)==expected,"Twenty-four pulse controls")
    q=data["momentum"]
    require(len(q)==data["nodes"] and q==sorted(set(q)) and 0<q[0]<q[-1]<8,"Momentum support")
    weights=quadrature_weights(q)
    physical={}
    for mesh in ("fine","coarse"):
        row=data["initial"][mesh]
        mu=abs(row["a"]*row["mass"])
        first=[w/math.hypot(k,mu) for w,k in zip(weights,q)]
        physical[mesh]=[first,[w*k*k for w,k in zip(first,q)]]
    for label,group in data["groups"].items():
        rows=group["rows"]
        target=list(range(24)) if label=="all_24" else [i for i,l in enumerate(labels) if l.startswith("fine/wide/number/")]
        require(rows==target,"Constraint indices")
        for kind,cert in certificates(group):
            lower,upper=cert["lower"],cert["upper"]
            close(ROOT2*lower,cert["quantum_lower"],"Quantum witness normalization")
            close(ROOT2*upper,cert["quantum_upper"],"Quantum dual normalization")
            close(lower,min(cert["margins"]),"Primal margin")
            require(-2e-10 <= upper-lower < 1e-4,"Primal-dual bracket")
            count=1 if kind in ARMS else len(rows)
            require(len(cert["margins"])==len(cert["dual"])==count,"Certificate dimensions")
            require(min(cert["dual"])>=0,"Positive dual weights")
            close(sum(cert["dual"]),1.,"Dual simplex")
            real,imag=cert["z_real"],cert["z_imag"]
            require(len(q)==len(real)==len(imag) and max(map(math.hypot,real,imag))<=1+2e-12,"Unit disk")
            meshes=("fine","coarse") if kind=="both_mesh_moments" else ("fine",)
            residual=max(abs(math.fsum(w*z for w,z in zip(f,real))/sum(f))
                         for mesh in meshes for f in physical[mesh])
            require(max(residual,cert["moment_residual"],cert["orthonormal_moment_residual"])<1e-9,
                    "Source-matched witness")
            if kind=="imaginary_only":
                require(max(map(abs,real))<1e-12,"Imaginary covariance")
            if kind in ("plus","minus"):
                require(cert["quantum_lower"]>1,"Individual near-arm quantum reversal")
            else:
                require(cert["quantum_upper"]<1,"Excluded joint or far count reversal")
            require(upper<1,"Separable threshold")


def crosscheck(grid,direct):
    require((grid["baseline"]["nodes"],grid["refined"]["nodes"],grid["high"]["nodes"])==(24,48,96),
            "Source-matching momentum refinements")
    for data in grid.values():
        check_primary(data)
    check_independent(direct)
    fine=grid["high"]["groups"]["fine_wide_full"]
    for arm in ARMS:
        bounds=fine["individual"][arm]
        measured=direct["128"][arm]["quantum_bound"]
        require(bounds["quantum_lower"]-2e-5 <= measured <= bounds["quantum_upper"]+2e-5,
                "Independent quantum support "+arm)
    close(fine["imaginary_only"]["quantum_upper"],direct["128"]["far"]["imag_quantum"],
          "Imaginary far control",rel=0.,atol=.001)
    for resolution in ("baseline","refined","high"):
        fine=grid[resolution]["groups"]["fine_wide_full"]
        close(fine["fine_moments"]["quantum_upper"],fine["individual"]["far"]["quantum_upper"],
              "Limiting far arm",atol=1e-5)
    high,ref=(grid[k]["groups"] for k in ("high","refined"))
    for group in high:
        for kind,cert in certificates(high[group]):
            other=ref[group]["individual"][kind] if kind in ARMS else ref[group][kind]
            close(cert["quantum_upper"],other["quantum_upper"],"Momentum-grid stability",
                  rel=0.,atol=.0015)


def compare_primary(actual,expected):
    for name in ("labels","nodes","samples","padding","vertices"):
        require(actual[name]==expected[name],"Replay settings "+name)
    for g in ("fine_wide_full","all_24"):
        for key,case in certificates(expected["groups"][g]):
            other=actual["groups"][g]["individual"][key] if key in ARMS else actual["groups"][g][key]
            # HiGHS can choose different points on an optimum face.
            for edge in ("lower","upper","quantum_lower","quantum_upper"):
                close(other[edge],case[edge],f"Replay {g}/{key}/{edge}",rel=0.,atol=2e-5)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository",type=Path,default=HERE.parent.parent,
                        help="HDblast root containing the moving-parent release and its pinned sources")
    parser.add_argument("--full",action="store_true",help="repeat all three primary grids and the independent transfer")
    args=parser.parse_args()
    repo=args.repository.resolve()
    authenticate(repo)
    saved={name:json.loads((HERE/path).read_bytes()) for name,(path,_) in RUNS.items()}
    direct=json.loads((HERE/"source_check/SOURCE_TRANSFER.json").read_bytes())
    crosscheck(saved,direct)
    if args.full:
        computed={name:run("source_match/source_match.py",repo,options)
                  for name,(_,options) in RUNS.items()}
        independent=run("source_check/source_transfer.py",repo,[])
        crosscheck(computed,independent)
        for name in saved:
            compare_primary(computed[name],saved[name])
        compare_direct(independent,direct)
    print(f"Verified initial stress and force moments, 24 pulse controls and independent "
          f"far-arm quantum support {direct['128']['far']['quantum_bound']:.9f}. Full numerical replay: {args.full}.")


if __name__=="__main__":
    main()
