#!/usr/bin/env python3
"""Independent exact-rational witness review. Standard library only.

No producer import, oracle call, float endpoint or numerical field solver.
The retained background and freshly enclosed bump jets remain shared premises.
"""
import argparse
import copy
import hashlib
import json
import sys
from fractions import Fraction as F
from math import isqrt
from pathlib import Path


class Failure(RuntimeError):
    pass


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Failure(message)


def need(ok, label):
    if not ok:
        raise Failure(label)


class I:
    def __init__(self, lo=0, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        need(self.lo <= self.hi, "reversed interval")

    @classmethod
    def read(cls, value):
        need(isinstance(value, dict) and "exact" in value, "missing exact endpoint")
        e = value["exact"]
        return cls(*e) if isinstance(e, list) else cls(e)

    def __add__(self, other):
        b = other if isinstance(other, I) else I(other)
        return I(self.lo + b.lo, self.hi + b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -(other if isinstance(other, I) else I(other))

    def __rsub__(self, other):
        return I(other) + -self

    def __mul__(self, other):
        b = other if isinstance(other, I) else I(other)
        p = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return I(min(p), max(p))
    __rmul__ = __mul__

    def __truediv__(self, other):
        b = other if isinstance(other, I) else I(other)
        need(not b.lo <= 0 <= b.hi, "division through zero")
        return self * I(1/b.hi, 1/b.lo)

    def __rtruediv__(self, other):
        return I(other) / self

    def __pow__(self, n):
        need(isinstance(n, int) and n >= 0, "invalid power")
        if n == 0:
            return I(1)
        if n % 2:
            return I(self.lo**n, self.hi**n)
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**n,self.hi**n),
                 max(self.lo**n,self.hi**n))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def contains(self, other):
        return self.lo <= other.lo and self.hi >= other.hi

    def exact(self):
        return [str(self.lo), str(self.hi)]


def scalar(raw):
    return F(raw["exact"])


def sqrt_cap(x):
    """Smallest dyadic upper square root on the frozen 2^-112 lattice."""
    scale = 1 << 112
    n = isqrt((x.numerator * scale**2) // x.denominator)
    if F(n*n, scale**2) < x:
        n += 1
    out = F(n, scale)
    need(out*out >= x and (out-F(1,scale))**2 < x, "sqrt rounding")
    return out


def outward_grid(interval):
    scale=1<<112
    lo=(interval.lo*scale).__floor__()
    hi=(interval.hi*scale).__ceil__()
    return I(F(lo,scale),F(hi,scale))


def read(path):
    return json.loads(path.read_text())


def digest(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def check_profile(row, shared, mutant=None):
    rho, a, p = (I.read(row["background"][k]) for k in ("rho","a","p"))
    rb, py = I.read(shared["rho_b"]), I.read(shared["p_y_global"])
    need(rho.lo > 0 and rb.lo > 0 and p.hi < 0, "background divisions/sign")
    ell = F(1, 1000000)
    q0,q1,q2,q3 = map(I.read, row["Q_r_jets"])
    R, b = (rb/rho)**3, ell*rho*a
    ay = -1/rho**2-p**2/3
    c, d = -3*rb**3/(2*ell*rho**4*p), ell*rho*py/p
    if mutant == "drop_p_y":
        d = I(0)
    if mutant == "drop_a_y":
        ay = I(0)
    g = q2+4*b*q1+3*ell**2*q0
    calculated = {
        "h": -R*q1,
        "h_r": -R*(q2+3*b*q1),
        "V_L": -R*(q2+4*b*q1),
        "s": c*g,
        "s_r": c*(q3+4*b*q2+3*ell**2*q1
                   -4*ell**2*rho**2*(ay+a**2)*q1+(4*b+d)*g),
    }
    for key, value in calculated.items():
        need(I.read(row["profiles"][key]).contains(value), "profile inclusion "+key)
    for key, value in {"R":R,"b":b,"c":c,"d":d,"a_y":ay}.items():
        need(I.read(row["auxiliary"][key]).contains(value), "auxiliary inclusion "+key)
    profiles = {k: I.read(v) for k,v in row["profiles"].items()}
    z = I(0)
    expected = {"U0":[profiles["h"],z,z,profiles["h"],profiles["s"]],
                "U0_r":[profiles["h_r"],z,z,profiles["h_r"],profiles["s_r"]],
                "V0":[z,profiles["V_L"],z,z,z]}
    for key, vector in expected.items():
        need(len(row[key]) == 5, "five initial entries")
        need([I.read(v).exact() for v in row[key]] == [v.exact() for v in vector],
             "five-component adapter "+key)


def aggregate(raw, base, old, profile_checks=True):
    ell, T, delta, eta, q = F(1,1000000), F(5,4), F(1,4), F(1,10**11), F(1201,3200000)
    rows = raw["cells"]
    need(len(rows)==241, "registered cell count")
    expected = [(F(1),F(65,64))]+[(F(k,1024),F(k+1,1024)) for k in range(1040,1280)]
    actual = [tuple(map(F,row["r"])) for row in rows]
    need(actual==expected, "registered contiguous cover including flat tail")
    need(rows[0]["tail"] and not any(row["tail"] for row in rows[1:]), "flat-tail placement")
    need(raw["zero_extension"]=={"r_ge_2":True,"r_le_1":True,"reason":"Analytic compactness, exact moment derivative identity"}, "zero wedge")
    target = raw["target"]
    # The frozen raw schema identifies this row by the fixed campaign target.
    # Its background must equal the retained exact T row below.
    old_by_r = {tuple(row["r"]):row for row in old["profile_cells"]}
    for row in rows:
        source = old_by_r[tuple(row["r"])]["profiles"]
        for k in ("rho","p","a"):
            need(I.read(row["background"][k]).exact()==I.read(source[k]).exact(), "retained cell background "+k)
    old_target = next(row for row in old["targets"] if F(row["r"])==T)
    for k in ("rho","p","a"):
        need(I.read(target["background"][k]).exact()==I.read(old_target["profiles"][k]).exact(), "retained target background "+k)
    if profile_checks:
        for row in rows+[target]:
            check_profile(row, raw["shared_inputs"])
    caps = {key:[max(I.read(row[key][j]).mag() for row in rows) for j in range(5)]
            for key in ("U0","U0_r","V0")}
    gram = [1,2,1,3,1]
    norms = {k:sqrt_cap(sum(w*v*v for w,v in zip(gram,values))) for k,values in caps.items()}
    N = max((norms["U0"]+T*norms["V0"])/3,norms["U0_r"]+norms["V0"])
    E = q*N/(1-q)
    jv,jt,av,at = I(0),I(0),F(0),F(0)
    for row,(l,u) in zip(rows,actual):
        s = I.read(row["profiles"]["s"])
        wt = u-l
        wv = wt*(T-(u+l)/2)
        need(wt>0 and wv>0,"positive history weights")
        jv,jt = jv+wv*s,jt+wt*s
        av,at = av+wv*s.mag(),at+wt*s.mag()
    ev,et = delta**3*E/6+eta*av,delta**2*E/2+eta*at
    hv,ht = jv+I(-ev,ev),jt+I(-et,et)
    ms = I.read(target["profiles"]["s"])+I(-delta*E,delta*E)
    shell = raw["shared_inputs"]["shell"]
    sig,sp,spp,up = (I.read(shell[k]) for k in ("sigma_b","sigma_prime_b","sigma_second_b","U_prime_b"))
    fp = 2*up-sp*spp/2+2*sig*sp/3
    rb = I.read(raw["shared_inputs"]["rho_b"])
    lv,lt = (rb*ell)**2*fp*hv/2,rb*ell*fp*ht/2
    need(ms.lo>0 and hv.lo>0 and ht.lo>0 and fp.hi<0 and lv.hi<0 and lt.hi<0,"strict registered signs")
    bv,bt = I.read(base["proper_v2"]),I.read(base["proper_Theta2"])
    threshold = max((-bv.lo)/(-lv.hi),(-bt.lo)/(-lt.hi))
    power,L = 0,F(1)
    while L<=threshold:
        power,L = power+1,L*10
    nv,nt = bv-L*lv,bt-L*lt
    need(nv.lo>0 and nt.lo>0,"finite-lambda positivity")
    return {"initial_component_caps":caps,"initial_Gram_norm_caps":{"D0":norms["U0"],"Dr":norms["U0_r"],"Dv":norms["V0"]},
            "N_h":N,"E_h":E,"J_volume":jv,"J_expansion":jt,"A_volume":av,"A_expansion":at,
            "volume_history_error":ev,"expansion_history_error":et,"kernel_volume_history":hv,"kernel_expansion_history":ht,
            "free_scalar_target":I.read(target["profiles"]["s"]),"M2_slope":ms,"Fprime":fp,"volume_slope":lv,"expansion_slope":lt,
            "finite_lambda_witness":{"lambda":-L,"L":L,"threshold":threshold,"smallest_power_of_ten_n":power,
                "base_volume":bv,"base_expansion":bt,"new_volume_second_coefficient":nv,"new_expansion_second_coefficient":nt}}


def comparison(calc, producer):
    count = 0
    def visit(a,b,path):
        nonlocal count
        if isinstance(a,I):
            need(a.exact()==I.read(b).exact(),"producer endpoint disagreement: "+path)
            count += 1
        elif isinstance(a,F):
            need(a==scalar(b),"producer rational disagreement: "+path)
            count += 1
        elif isinstance(a,dict):
            for key,value in a.items():visit(value,b[key],path+"/"+key)
        elif isinstance(a,list):
            need(len(a)==len(b),"producer list length: "+path)
            for j,value in enumerate(a):visit(value,b[j],path+"/"+str(j))
        else:
            need(a==b,"producer integer disagreement: "+path)
            count += 1
    visit(calc,producer,"")
    return count


def encode(value):
    if isinstance(value,I):return {"exact":value.exact()}
    if isinstance(value,F):return {"exact":str(value)}
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,list):return [encode(v) for v in value]
    return value


def controls(raw, base, old, calc, producer):
    results = {}
    def rejected(name, operation):
        try:
            operation()
        except Failure as exc:
            results[name] = {"status":"REJECTED", "reason":str(exc)}
        else:
            raise Failure("mutant survived: "+name)
    m=copy.deepcopy(raw)
    m["cells"]=m["cells"][1:]
    rejected("missing_flat_tail",lambda:aggregate(m,base,old))
    m=copy.deepcopy(raw)
    m["cells"][1]["r"][0]="1"
    rejected("overlapping_history_cover",lambda:aggregate(m,base,old))
    for name, key, slot in [("omit_mixed_velocity","V0",1),
                            ("omit_scalar_pullback","U0",4),
                            ("spurious_lapse_difference","U0",2)]:
        m=copy.deepcopy(raw)
        m["target"][key][slot]={"exact":["1","1"]} if slot==2 else {"exact":["0","0"]}
        rejected(name,lambda:aggregate(m,base,old))
    rejected("drop_p_y_in_scalar_derivative",lambda:check_profile(raw["target"],raw["shared_inputs"],"drop_p_y"))
    rejected("drop_a_y_in_scalar_derivative",lambda:check_profile(raw["target"],raw["shared_inputs"],"drop_a_y"))
    for name, key, value in [
        ("omit_coupled_resolvent_error","E_h",{"exact":"0"}),
        ("wrong_Gram_L_weight","initial_Gram_norm_caps",{**producer["initial_Gram_norm_caps"],"Dv":{"exact":"0"}}),
        ("wrong_positive_volume_weight","J_volume",producer["J_expansion"]),
        ("extra_second_factorial","volume_slope",encode(2*calc["volume_slope"])),
        ("missing_expansion_clock_factor","expansion_slope",encode(calc["expansion_slope"]/F(1,1000000))),
        ("wrong_shell_derivative_sign","Fprime",encode(-calc["Fprime"]))]:
        p=copy.deepcopy(producer)
        p[key]=value
        rejected(name,lambda:comparison(calc,p))
    p=copy.deepcopy(producer)
    p["finite_lambda_witness"]["lambda"]={"exact":"1000000000000000000000"}
    rejected("wrong_finite_lambda_sign",lambda:comparison(calc,p))
    return results


def main():
    ap=Parser()
    ap.add_argument("--numeric-root",type=Path,required=True)
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--raw",type=Path)
    ap.add_argument("--producer-report",type=Path)
    ap.add_argument("--comparison-raw",type=Path)
    ap.add_argument("--comparison-report",type=Path)
    args=ap.parse_args()
    nr,repo=args.numeric_root,args.repo_root
    destination=args.output.resolve()
    need(not destination.exists(),"fresh output required")
    need(not destination.is_relative_to(repo.resolve()),"output must be outside source repository")
    need(not destination.is_relative_to(nr.resolve()),"output must be outside numeric source component")
    capsule=repo/"research/HDBLAST_PREPARED_SHELL_VOLUME_20261101"
    # Bind exact independent review inputs; do not execute retained authentication modules.
    paths={"campaign":nr/"CAMPAIGN_FREEZE.json","postprocess":nr/"POSTPROCESS_FREEZE.json",
           "raw":args.raw or nr/"normal/RAW_WITNESS.json","producer_report":args.producer_report or nr/"normal/REPORT.json",
           "retained_initial_boxes":capsule/"initial/references/RAW_BOXES.json",
           "retained_base_response":capsule/"response/FROZEN_SCIENTIFIC_REFERENCE.json",
           "retained_static_jets":repo/"research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json"}
    pins={"campaign":"83907b72a6b92c41e40568266077f6db2fbaecd43d7b5db3d72babf4c3d821d1",
          "postprocess":"2325969b88bcbfb1464424b60f228433343905e3fb147b072a2f20a231be0d36",
          "retained_initial_boxes":"20c2e33ab7721860171ee7931a3c4d07cb18c597d98848f0096d5b117a60faf3",
          "retained_base_response":"2c9bc1d376d6fc4e77a853fdc0548c814a5253e0b5d78c0a1fd4d5ea4695a22e",
          "retained_static_jets":"3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22"}
    for key,sha in pins.items():need(digest(paths[key])["sha256"]==sha,"source bytes: "+key)
    comparison_raw=args.comparison_raw or (nr/"optimized/RAW_WITNESS.json" if not args.raw else None)
    comparison_report=args.comparison_report or (nr/"optimized/REPORT.json" if not args.producer_report else None)
    if comparison_raw:need(paths["raw"].read_bytes()==comparison_raw.read_bytes(),"producer normal/O raw equality")
    if comparison_report:need(paths["producer_report"].read_bytes()==comparison_report.read_bytes(),"producer normal/O report equality")
    raw,report,old,base=(read(paths[k]) for k in ("raw","producer_report","retained_initial_boxes","retained_base_response"))
    need(raw["campaign_freeze_sha256"]==digest(paths["campaign"])["sha256"],"campaign binding")
    need(F(read(paths["campaign"])["target"])==F(5,4),"frozen target")
    st=read(paths["retained_static_jets"])
    need(I.read(raw["shared_inputs"]["p_y_global"]).exact()==outward_grid(I.read(st["y_jets"]["p"][1])).exact(),"retained p_y")
    need(I.read(raw["shared_inputs"]["rho_b"]).exact()==I.read(base["first"]["rho"]).exact(),"retained rho_b")
    calc=aggregate(raw,base,old)
    need(calc["Fprime"].exact()==I.read(base["first"]["Wprime"]).exact(),"retained Fprime")
    matched=comparison(calc,report)
    out={"schema":"HDBLAST_PREPARATION_INDEPENDENT_RATIONAL_REAGGREGATION_V1","status":"PASS",
         "source_inputs":{k:digest(p) for k,p in paths.items()},"producer_scientific_fields_matched":matched,
         "whole_cell_and_target_profile_checks":242,"cells":241,
         "negative_controls":controls(raw,base,old,calc,report),
         "independence":{"producer_modules_imported":False,"oracle_calls":0,"arithmetic":"stdlib Fraction and integer square root",
          "shared_premises":["Retained background boxes and fresh whole-cell bump derivative enclosures","Analytic constraint, gauge, five-mean PDE, finite-propagation and development theorems"],
          "review_limit":"Independent AI arithmetic and derivation, not external specialist review or a new PDE observation"},
         "science":encode(calc)}
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open("x",encoding="utf-8") as stream:
        stream.write(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":out["status"],"matched_fields":matched,"negative_controls":len(out["negative_controls"])}))


if __name__=="__main__":
    try:
        main()
    except (Failure,ValueError,KeyError,TypeError,OSError) as error:
        print(json.dumps({"status":"FAIL","error":str(error)},sort_keys=True))
        sys.exit(2)
