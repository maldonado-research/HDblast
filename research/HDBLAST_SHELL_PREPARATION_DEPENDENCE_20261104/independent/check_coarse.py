#!/usr/bin/env python3
"""Independent standard-library audit of the elementary P1 certificate.

This reconstructs exact bounds from pinned retained JSON and compares the
separately produced coarse report. It does not import the coarse producer.
"""
import argparse
import copy
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path


class Parser(argparse.ArgumentParser):
    def error(self,message):
        raise ValueError(message)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def endpoints(value):
    return tuple(map(F,value["exact"]))


def enc(value):
    if isinstance(value,tuple):return {"exact":list(map(str,value))}
    if isinstance(value,F):return {"exact":str(value)}
    if isinstance(value,dict):return {k:enc(v) for k,v in value.items()}
    return value


def compare(calc, produced, path=""):
    count=0
    for key,value in calc.items():
        label=path+"/"+key
        reported=produced[key]
        if isinstance(value,F):
            need(value==F(reported["exact"]),"coarse rational disagreement: "+label)
            count+=1
        elif isinstance(value,tuple):
            need(value==endpoints(reported),"coarse interval disagreement: "+label)
            count+=1
        elif isinstance(value,dict):
            count+=compare(value,reported,label)
        else:
            need(value==reported,"coarse exact disagreement: "+label)
            count+=1
    return count


def exp_partial(x,n):
    term=F(1)
    total=term
    for k in range(1,n+1):
        term*=x/k
        total+=term
    return total,term*x/(n+1)/(1-x/(n+2))


def bump_proof(third_coefficients=(1,-10,24,-12)):
    # Equality at degree+1 distinct rational nodes proves polynomial equality.
    for d in [F(0),F(1,16),F(1,8),F(3,16)]:
        p2=1-6*d+6*d*d
        second_direct=(1-4*d)*(1-2*d)-2*d*d
        need(p2==second_direct,"Q second derivative identity")
        differentiated=p2+d*d*(-6+12*d)-4*d*p2
        claimed=sum(F(c)*d**j for j,c in enumerate(third_coefficients))
        need(differentiated==claimed,"Q third derivative identity")
    dmax=F(3,16)
    need(-6+12*dmax<0 and -10+48*dmax<0,"continuous polynomial monotonicity")
    need(1-6*dmax+6*dmax**2==F(11,128),"Q second numerator lower bound")
    need(1-10*dmax+24*dmax**2-12*dmax**3==-F(113,1024),"Q third numerator lower bound")
    low,remainder=exp_partial(F(4,3),8)
    need(low>F(15,4) and low+remainder<4,"continuous exp(-4/3) bounds")
    low2,_=exp_partial(F(2),9)
    need(low2>F(6**6,6400),"continuous third derivative maximum")
    # w>=16/3; w^j exp(-w) decreases for j=2,4, and peaks at w=6 for j=6.
    need(F(4,15)<F(1,3),"Q bound")
    need(F(4,15)*F(16,3)**2<8,"Q prime bound")
    need(F(4,15)*F(16,3)**4<216,"Q second bound")
    need(F(5632,81)/4==F(1408,81),"target scalar factor")
    need(F(128,9)/4==F(32,9),"unweighted primitive factor")


def recompute(static,base):
    bump_proof()
    ell,T,delta=F(1,1000000),F(5,4),F(1,4)
    q,eta=F(1201,3200000),F(1,10**11)
    rho=(F(129),F(131)); pabs=(F(3,25000),F(1,7000))
    amax,aymax,pymax=F(14,125),F(1,10000),F(1,4000)
    def within(v,lo,hi,label):
        a,b=endpoints(v)
        need(lo<=a<=b<=hi,"coarse background premise "+label)
    for field,order,lo,hi in [("rho",0,*rho),("p",0,-pabs[1],-pabs[0]),
                             ("p",1,-pymax,pymax),("a",0,F(111,1000),amax),
                             ("a",1,-aymax,aymax)]:
        within(static["y_jets"][field][order],lo,hi,field+str(order))
    within(base["first"]["rho"],*rho,"rho_b")
    within(base["first"]["Wprime"],-F(1,1000),-F(1,1200),"Fprime")
    need(F(static["domain"]["t"])==0 and F(static["domain"]["ell"])==ell,"static domain normalization")
    need(F(static["domain"]["z_min"])<=ell and F(static["domain"]["z_max"])>=ell*T,"static domain covers cone")
    need(F(base["q"]["exact"])==q,"full coupled q premise")
    rmax=(rho[1]/rho[0])**3
    bmax=ell*rho[1]*amax
    dmax=ell*rho[1]*pymax/pabs[0]
    cmin=F(3,2)*rho[0]**3/(ell*rho[1]**4*pabs[1])
    cmax=F(3,2)*rho[1]**3/(ell*rho[0]**4*pabs[0])
    q0,q1,q2,q3=F(1,3),F(8),F(216),F(6400)
    hmax=rmax*q1
    hrmax=rmax*(q2+3*bmax*q1)
    vlmax=rmax*(q2+4*bmax*q1)
    gmax=q2+4*bmax*q1+3*ell**2*q0
    smax=cmax*gmax
    # Expand the product-rule derivative into its seven nonnegative allowances.
    srmax=cmax*(q3+4*bmax*q2+3*ell**2*q1+4*bmax*gmax+dmax*gmax
                 +4*ell**2*rho[1]**2*amax**2*q1+4*ell**2*rho[1]**2*aymax*q1)
    D0,Dr,Dv=2*hmax+smax,2*hrmax+srmax,2*vlmax
    N=max((D0+T*Dv)/3,Dr+Dv)
    E=q*N/(1-q)
    hv=((1-eta)*cmin/4-delta**3*E/6,delta**2*smax/2+delta**3*E/6)
    ht=((1-eta)*cmin*F(32,9)-delta**2*E/2,delta*smax+delta**2*E/2)
    m=(F(1408,81)*cmin-delta*E,smax+delta*E)
    need(min(m[0],hv[0],ht[0])>0,"strict coarse histories")
    lv=(-(rho[1]*ell)**2*hv[1]/2000,-(rho[0]*ell)**2*hv[0]/2400)
    lt=(-rho[1]*ell*ht[1]/2000,-rho[0]*ell*ht[0]/2400)
    bv=(F('-5.4520e16'),F('-2.6521e16'))
    bt=(F('-9.0208e21'),F('-4.4458e21'))
    within(base["proper_v2"],*bv,"base volume")
    within(base["proper_Theta2"],*bt,"base expansion")
    threshold=max(bv[0]/lv[1],bt[0]/lt[1])
    n,L=0,F(1)
    while L<=threshold:n,L=n+1,L*10
    nv=(bv[0]-L*lv[1],bv[1]-L*lv[0])
    nt=(bt[0]-L*lt[1],bt[1]-L*lt[0])
    need(nv[0]>0 and nt[0]>0,"strict coarse simultaneous reversal")
    return {"constants":{"q":q,"eta":eta,"c_lo":cmin,"c_hi":cmax},
            "analytic_bump_caps":{"Q":"1/3","Q_r":"8","Q_rr":"216","Q_rrr":"6400"},
            "initial_norm_caps":{"D0":D0,"Dr":Dr,"Dv":Dv,"N":N},
            "full_resolvent_error":E,"positive_history_lower_bounds":{"v":hv[0],"Theta":ht[0]},
            "slopes":{"M2":m,"v2":lv,"Theta2":lt},
            "deterministic_choice":{"threshold":threshold,"n":n,"lambda":str(-L)},
            "new_coefficients":{"v2":nv,"Theta2":nt}}


def main():
    ap=Parser()
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--coarse-report",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    destination=args.output.resolve()
    need(not destination.exists(),"fresh output required")
    need(not destination.is_relative_to(args.repo_root.resolve()),"output must be outside source repository")
    need(destination!=args.coarse_report.resolve(),"output must differ from coarse input")
    source_specs={
      "static":("research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json","3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22"),
      "base":("research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json","1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3")}
    sources={}
    for key,(path,sha) in source_specs.items():
        payload=(args.repo_root/path).read_bytes()
        need(hashlib.sha256(payload).hexdigest()==sha,"coarse source bytes "+key)
        sources[key]=json.loads(payload)
    produced=json.loads(args.coarse_report.read_text())
    calc=recompute(sources["static"],sources["base"])
    matched=compare(calc,produced)
    mutants={}
    def reject(name,fn):
        try:fn()
        except ValueError as error:mutants[name]={"status":"REJECTED","reason":str(error)}
        else:raise ValueError("coarse mutant survived "+name)
    reject("wrong_Q_third_polynomial",lambda:bump_proof((1,-12,36,-24)))
    for name,path,value in [
      ("omit_scalar_derivative_cap",("analytic_bump_caps","Q_rrr"),"0"),
      ("omit_mixed_velocity_norm",("initial_norm_caps","Dv"),{"exact":"0"}),
      ("omit_coupled_error",("full_resolvent_error",),{"exact":"0"}),
      ("omit_kernel_lower_loss",("constants","eta"),{"exact":"0"}),
      ("wrong_finite_lambda",("deterministic_choice","lambda"),"1000000000000000000000"),
      ("extra_volume_factorial",("slopes","v2"),enc(tuple(2*x for x in calc["slopes"]["v2"]))),
      ("missing_proper_clock",("slopes","Theta2"),enc(tuple(x*10**6 for x in calc["slopes"]["Theta2"])) )]:
        changed=copy.deepcopy(produced)
        target=changed
        for key in path[:-1]:target=target[key]
        target[path[-1]]=value
        reject(name,lambda:compare(calc,changed))
    result={"status":"PASS_INDEPENDENT_COARSE_REVIEW","matched_scientific_fields":matched,
            "source_sha256":{key:sha for key,(_,sha) in source_specs.items()},
            "coarse_report_sha256":hashlib.sha256(args.coarse_report.read_bytes()).hexdigest(),
            "producer_imported":False,"oracle_calls":0,
            "continuous_bound_review":"Exact polynomial identities, continuous extrema, Taylor exponential enclosures, positive primitive moments, full Gram caps and coupled error checked.",
            "negative_controls":mutants,"science":enc(calc)}
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open("x",encoding="utf-8") as stream:
        stream.write(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"matched_scientific_fields":matched,"controls":len(mutants)}))


if __name__=="__main__":
    try:
        main()
    except (ValueError,KeyError,TypeError,OSError) as error:
        print(json.dumps({"status":"FAIL","error":str(error)},sort_keys=True))
        sys.exit(2)
