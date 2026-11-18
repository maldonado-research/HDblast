"""Independent endpoint and integral certificates for two source accelerations."""
import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import argparse

PINS = {
    "inputs/SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "METRIC.json": "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42",
    "inputs/SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "inputs/FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
}


def mul(x, y):
    products = [u*v for u in x for v in y]
    return min(products), max(products)


def source_bounds(parent):
    for name, sha in PINS.items():
        if hashlib.sha256((parent/name).read_bytes()).hexdigest() != sha:
            raise ValueError(name)
    metric = json.loads((parent/"METRIC.json").read_text())
    b, a, ell = F(99,64), F(49,32), F(1,10**6)
    cmax, transfer = 1+F(1,10**11), F(3,10**6)
    m_error = F(metric["absolute_mean_error"])
    qm = [F(0), F(0)]
    im = [F(0), F(0)]
    count = 0
    edge = F(1)
    with (parent/"inputs/FUTURE_CELLS.csv").open() as fp:
        for row in csv.DictReader(fp):
            left=max(F(int(row["left_grid"]),4096),1)
            right=min(F(int(row["right_grid"]),4096),b)
            if right<=left:
                continue
            if left!=edge:
                raise ValueError("coverage")
            lo=F(int(row["H_lo"]),2**112)-transfer
            hi=F(int(row["H_hi"]),2**112)+transfer
            # Integral of u times (b-t) or (b-t)^2/2 across a whole certified cell.
            end_kernel=((b-left)**2-(b-right)**2)/2
            int_kernel=((b-left)**3-(b-right)**3)/6
            for i, acc in enumerate([lo,hi]):
                qm[i]+=end_kernel*acc
                im[i]+=int_kernel*acc
            edge=right
            count+=1
    if edge!=b:
        raise ValueError("endpoint")
    m=(qm[0]-m_error,qm[1]+m_error)
    integ=(im[0]-(b-1)*m_error,im[1]+(b-1)*m_error)
    tq=(1/(b*cmax**5), (1+cmax**2*ell**2*a*b)/b)
    span=(b-a,cmax*(b-a))
    raw=mul(span,m)
    static=mul(tq,(-raw[1],-raw[0]))
    velocity=mul((F(2),2*cmax),integ)
    free=mul(tq,(velocity[0]-raw[1],velocity[1]-raw[0]))
    # X=ell/(mu*rho_b^2)=ell B_b/rho_b converts back to delta D_A / L.
    box=(parent/"inputs/SHELL_BOX.json")
    bins=json.loads(box.read_text())["boxes"]
    blo,bhi=map(F,bins["B_b"]["exact"])
    rlo,rhi=map(F,bins["rho_b"]["exact"])
    gain=(ell*blo/rhi,ell*bhi/rlo)
    cap=(F(18,25),)*2
    # Baseline z = -(cosh(q1)-cosh(q0))/cosh(q1).
    kernel=ell**2*(b*b-a*a)/2
    bgz=(-cmax*kernel,-kernel/cmax)
    K0=1/(ell**2*b*(b-a))
    K=(K0/cmax**2,K0*cmax**2)
    return {"cells":count,"fixed_redshift_amplification":K,"endpoint_M":m,"integral_M":integ,
            "normalized_supported":static,"normalized_free":free,
            "fractional_supported_cap":mul(cap,mul(gain,static)),
            "fractional_free_cap":mul(cap,mul(gain,free)),
            "baseline_z":bgz}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("parent", type=Path, nargs="?", default=Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118")
    result=source_bounds(parser.parse_args().parent)
    print(json.dumps(result,indent=2,default=str))
