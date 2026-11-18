"""Reproduce the selected B1 runs with read-only collar output at the original record calls.

Pass the public HDblast checkout and a fresh private output directory. No pickle
cache is imported. All original numerical modules are authenticated before use.
"""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys

PREFIX = "research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/"
PINS = {
    PREFIX+"evolve_b1.py":"a0ae79e12ff90f549a0d1ca674fa3285b8ad456b780880df405b9597d3434f0f",
    PREFIX+"static_w.py":"ac4f5df5d8a1fcedfbd9503caf4717c1473fde6699388e095618f269b6577b1f",
    PREFIX+"shell_matter.py":"58c4c4cc18132f86aaa45ac9f0930490ef3debe160fe8f6ed59274451aadc7d1",
    "research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/rolloff5d_v1.py":"9ea4d8d6ef36ce3d242fbb0e09f27a944b571fdfb8646e3d71195b60b089080b",
}
INJECT = '''
            # Read-only diagnostic at the archived solver's ordinary record calls.
            mask = self.z >= -0.081
            dest = Path(out_dir) / "collar"
            dest.mkdir(parents=True,exist_ok=True)
            fields = {name:F[name][mask].copy() for name in ("A","B","phi","A_T","B_T","phi_T")}
            np.savez(dest / ("sample_%06d.npz" % len(rows)), z=self.z[mask],
                     T=T, s=tau, sigma1=row["sigma1"], J=row["J"], R=row["R"],
                     X=row["X"], H_over_H0=row["H_over_H0"], **fields)
'''
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo",type=Path)
    parser.add_argument("directory",type=Path)
    parser.add_argument("grid",choices=["coarse","fine"])
    p=parser.parse_args()
    dest=p.directory.resolve()
    dest.mkdir(parents=True,exist_ok=False)
    for rel,digest in PINS.items():
        raw=(p.repo/rel).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=digest:
            raise ValueError("published input hash: "+rel)
        if rel.endswith("/evolve_b1.py"):
            key=b"            return row, F\n"
            if raw.count(key)!=1:
                raise ValueError("original read-only diagnostic location")
            raw=raw.replace(key,INJECT.encode("ascii")+key)
        (dest/Path(rel).name).write_bytes(raw)
    fine=p.grid=="fine"
    args=[sys.executable,"-B","-W","ignore","evolve_b1.py","--d","-3.106933495673783",
          "--source","chi","--dc","0.01","--phistar","0.5","--G","100","--y","1","--b","6.49e-6",
          "--L","16","--Tf","2.43","--kappa","10","--project","15.9","--xc","1.9",
          "--dzf","0.0005" if fine else "0.001","--dzc","0.002" if fine else "0.004",
          "--out","run","--tag","b1_"+p.grid,"--wall","3600"]
    env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (dest/"log.txt").open("w") as log:
        subprocess.run(args,cwd=dest,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)

if __name__=="__main__":
    main()
