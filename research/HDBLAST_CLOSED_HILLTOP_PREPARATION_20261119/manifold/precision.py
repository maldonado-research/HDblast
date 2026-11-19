"""Direct absolute-coordinate decimal Taylor shoots on the two spline cells adjoining the top."""
import argparse
import importlib.util
import json
import math
from decimal import Decimal as D, getcontext
from pathlib import Path

ZERO, ONE = D(0), D(1)


def product(a,b):
    values=[ZERO]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            values[i+j]+=x*y
    return values


def quotient(num,den,count):
    numbers=[]
    for i in range(count):
        val=num[i] if i<len(num) else ZERO
        for j in range(1,min(i,len(den)-1)+1):
            val-=den[j]*numbers[i-j]
        numbers.append(val/den[0])
    return numbers


def make_chart(row,k,degree):
    # Convert each binary spline coefficient exactly to decimal before shifting.
    off=D(row["top"])-D(row["spline_knots"][k])
    p,f=(list(reversed([D(c[k][i]) for c in row["phi_f_coefficients"]])) for i in (0,1))
    def shift(poly):
        return [sum(poly[j]*math.comb(j,i)*(off**(j-i) if j>i else ONE)
                    for j in range(i,4)) for i in range(4)]
    p,f=shift(p),shift(f)
    pm,b,v0=map(D,(row["minimum"],row["b"],row["v0"]))
    z=[-v/pm for v in p]; z[0]+=ONE
    w=[v*b for v in p]; w[0]+=ONE
    num=[v/v0 for v in product(product(z,z),w)]
    u=quotient(num,product(f,f),degree+2)
    # Same centered force convention as the two binary64 adjacent spline charts.
    force=[-3*(i+1)*u[i+1] for i in range(degree+1)]
    corr=quotient([force[0]*f[0]**3],product(product(f,f),f),degree+1)
    return u[:degree+1],[v-corr[i] for i,v in enumerate(force)]


def evaluate(coeffs,x):
    s=coeffs[-1]
    for i in range(len(coeffs)-2,-1,-1):
        s=s*x+coeffs[i]
    return s


def launch(charts,start,velocity,end,rate,order):
    # Taylor coefficients are in powers of the fractional step x/dx. This
    # avoids large coefficient ranges when constructing nonlinear products.
    top,left,right=charts
    chart=left if start < -top else right
    degree=len(chart[0])-1
    rho=evaluate(chart[0],start)+velocity**2/6
    y=[start,velocity,ZERO,-rho.ln()/2]
    nsteps=int(math.ceil(float(end)*rate))
    step=end/nsteps
    weights=[step/(n+1) for n in range(order)]
    for it in range(nsteps):
        q,p,h,l=([val] for val in y)
        if abs(q[0])>D(".0005"):
            raise ValueError("Shoot left the adjacent spline cells")
        u,force=left if q[0] < -top else right
        powers=[[ONE],q]+[[q[0]**k] for k in range(2,degree+1)]
        for n,w in enumerate(weights):
            if n:
                for k in range(2,degree+1):
                    powers[k].append(sum(a*b for a,b in zip(q,reversed(powers[k-1]))))
            ucoef=sum(u[k]*powers[k][n] for k in range(1,degree+1))
            fcoef=sum(force[k]*powers[k][n] for k in range(1,degree+1))
            if n==0:
                ucoef+=u[0]
                fcoef+=force[0]
            pp=sum(a*b for a,b in zip(p,reversed(p)))
            hh=sum(a*b for a,b in zip(h,reversed(h)))
            hp=sum(a*b for a,b in zip(h,reversed(p)))
            q.append(w*p[n])
            p.append(w*(fcoef-3*hp))
            h.append(w*(ucoef-pp/3-hh))
            l.append(w*h[n])
        y=[sum(s) for s in (q,p,h,l)]
    return y


def check(row,reference,rate,order,degree):
    top=D(row["top"])
    center_ix=next(i-1 for i,k in enumerate(row["spline_knots"]) if D(k)>top)
    charts=(top,make_chart(row,center_ix-1,degree),make_chart(row,center_ix,degree))
    v=D("-.001")
    point=next(s for s in reference["samples"] if s["v"]==-.001)
    center=D(point["center"])
    decay=(-3-(D(9)+4*D(row["source_mu_squared"])).sqrt())/2
    # At successive endpoints suppress the growing solution. The endpoint
    # shifts are exponentially small once the decaying field is close to zero.
    ends=[D("1.3"),D("2.25"),D("2.85")]
    differences=[]
    for end in ends:
        old=center
        def residual(c):
            q,p,_,_=launch(charts,c,v,end,rate,order)
            return p-decay*q
        inc=D("1e-28")
        r=residual(center)
        derivative=(residual(center+inc)-r)/inc
        for _ in range(4):
            delta=r/derivative
            center-=delta
            if abs(delta)<D("1e-34"):break
            r=residual(center)
        differences.append(str(center-old))
    sign=D(-1 if row["minimum"]>0 else 1)
    width=D(point["conditional_width"])
    displaced=center+sign*width
    assert displaced!=center
    q,p,h,l=launch(charts,displaced,v,D(point["match_time"]),rate,order)
    return dict(center_decimal=str(center),signed_separation_decimal=str(displaced-center),
                finite_T_changes=differences,match_q_relative_error=float(q/D(reference["actual_zero_v_match_q"])-1),
                match_loga_absolute_error=float(l-D(reference["match_loga"])),
                match_q=float(q),match_p=float(p))


def main():
    getcontext().prec=50
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior",type=Path,default=here.parent/"prior")
    parser.add_argument("--manifold",type=Path,default=here/"MANIFOLD.json")
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    ref=json.loads(args.manifold.read_text())["branches"]
    spec=importlib.util.spec_from_file_location("extract_local_spline",here/"extract.py")
    ex=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    result={}
    for k,row in ex.extract(args.prior)["branches"].items():
        coarse=check(row,ref[k],20,28,10)
        fine=check(row,ref[k],24,32,12)
        gap=float(D(fine["center_decimal"])-D(coarse["center_decimal"]))
        if (abs(gap)>ref[k]["critical_zero_v_distance"]*1e-5 or
            abs(fine["match_q_relative_error"])>2e-7 or abs(fine["match_loga_absolute_error"])>2e-9):
            raise ArithmeticError((k,gap,fine))
        result[k]=dict(v=-.001,coarse=coarse,fine=fine,center_refinement_difference=gap)
        print(k,fine,flush=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
