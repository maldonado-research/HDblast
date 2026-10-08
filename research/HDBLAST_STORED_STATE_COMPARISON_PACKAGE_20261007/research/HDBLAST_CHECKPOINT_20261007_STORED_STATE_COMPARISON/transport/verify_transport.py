#!/usr/bin/env python3
"""Stdlib exact Laurent identities and fabricated transport controls; no real input."""
from __future__ import annotations
import argparse
from dataclasses import replace
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

from discrete_transport import (CAPSULES, EPSILON, PI, LA, LB, Interval, NodeError,
                                PrefixAccumulator, aggregate_capsule, aggregate_registered,
                                encode_exact, from_saved_target, node_terms)

NAMES = ("k", "L", "x", "v", "z", "c", "p", "r", "d", "h")
N = len(NAMES)

class Poly:
    """Canonical finite Laurent polynomial over Q; algebra is exact."""
    def __init__(self, val=0):
        if isinstance(val, Poly):
            self.t = dict(val.t)
        elif isinstance(val, dict):
            self.t = {m: F(c) for m, c in val.items() if c}
        else:
            self.t = {(0,)*N: F(val)} if val else {}
    @classmethod
    def symbol(cls, name):
        e = [0]*N; e[NAMES.index(name)] = 1
        return cls({tuple(e): F(1)})
    def __add__(self, b):
        b, terms = Poly(b), dict(self.t)
        for m, c in b.t.items(): terms[m] = terms.get(m, F(0))+c
        return Poly(terms)
    __radd__ = __add__
    def __neg__(self): return Poly({m: -c for m, c in self.t.items()})
    def __sub__(self, b): return self+-Poly(b)
    def __rsub__(self, b): return Poly(b)+-self
    def __mul__(self, b):
        terms = {}
        for m, c in self.t.items():
            for n, d in Poly(b).t.items():
                key = tuple(a+b for a, b in zip(m, n))
                terms[key] = terms.get(key, F(0))+c*d
        return Poly(terms)
    __rmul__ = __mul__
    def __pow__(self, e):
        if type(e) is not int: raise TypeError("integer exponent required")
        if e < 0:
            if len(self.t) != 1: raise ValueError("only monomials invertible")
            m, c = next(iter(self.t.items()))
            return Poly({tuple(e*x for x in m): c**e})
        out = Poly(1)
        for unused in range(e): out = out*self
        return out
    def __truediv__(self, b): return self*Poly(b)**-1
    def diff(self, name):
        i, terms = NAMES.index(name), {}
        for m, c in self.t.items():
            if m[i]:
                new = list(m); new[i] -= 1
                terms[tuple(new)] = c*m[i]
        return Poly(terms)
    def subs(self, replacements):
        out = Poly()
        for m, c in self.t.items():
            term = Poly(c)
            for name, exponent in zip(NAMES, m):
                if exponent:
                    term *= Poly(replacements.get(name, Poly.symbol(name)))**exponent
            out += term
        return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--full-fabricated", action="store_true")
    args = parser.parse_args()
    if args.output.exists(): raise RuntimeError("fresh output required")
    checks, mutations = [], []
    def check(name, condition):
        if not condition: raise RuntimeError(f"check failed: {name}")
        checks.append(name)
    def zero(name, poly): check(name, not poly.t)
    def rejects(name, function):
        try: function()
        except (ValueError, TypeError, KeyError): mutations.append(name)
        else: raise RuntimeError(f"mutation escaped: {name}")
    def different(name, a, b):
        if a == b: raise RuntimeError(f"mutation escaped: {name}")
        mutations.append(name)
    k,L,x,v,z,c,p,r,d,h = [Poly.symbol(name) for name in NAMES]
    R = ((2*k*k+3*L*L)*x-k*z-L*v)/(2*k)
    P = ((2*k*k/3-L*L)*x-k*z-L*v)/(2*k)
    def deriv(f):
        return f.diff("x")*v+f.diff("v")*(-2*k*z)+f.diff("z")*(2*k*v)+f.diff("L")*L*L
    zero("direct homogeneous continuity", deriv(R)-L*(R-3*P))
    zero("direct pressure primitive", deriv(-R/(3*L))-P)
    zero("canonical c invariant", deriv(x-z/(2*k)))
    zero("canonical direct density", R.subs({"x":c+p,"z":2*k*p,"v":-2*k*r})-((k+3*L*L/(2*k))*c+3*L*L*p/(2*k)+L*r))
    zero("canonical direct pressure", P.subs({"x":c+p,"z":2*k*p,"v":-2*k*r})-((k/3-L*L/(2*k))*c-(2*k/3+L*L/(2*k))*p+L*r))
    zero("pressure majorant square positive gap", (2*k/3+5*L*L/(4*k))**2-((2*k/3+L*L/(2*k))**2+L*L)-21*L**4/(16*k*k))
    zero("density triangle majorant square positive gap", (L+3*L*L/(2*k))**2-(L*L+9*L**4/(4*k*k))-3*L**3/k)
    # Here h is the exact retained weight/(4 represented-Pi^2), mu=2*h*k^2.
    fusedR = h*((k*k+3*L*L/2)*d+3*L*L*z/2-L*k*v)
    fusedP = h*((k*k/3-L*L/2)*d-(2*k*k/3+L*L/2)*z-L*k*v)
    zero("fused direct density anchor", 2*h*k*k*R.subs({"x":(d+z)/(2*k)})-fusedR)
    zero("fused direct pressure anchor", 2*h*k*k*P.subs({"x":(d+z)/(2*k)})-fusedP)
    zero("fused density canonical majorant", 2*h*k*k*(k+3*L*L/(2*k))*d/(2*k)-h*(k*k+3*L*L/2)*d)
    zero("fused density amplitude majorant", 2*h*k*k*(L+3*L*L/(2*k))*z/(2*k)-h*(L*k+3*L*L/2)*z)
    zero("fused pressure canonical majorant", 2*h*k*k*(k/3+L*L/(2*k))*d/(2*k)-h*(k*k/3+L*L/2)*d)
    zero("fused pressure amplitude majorant", 2*h*k*k*(2*k/3+5*L*L/(4*k))*z/(2*k)-h*(2*k*k/3+5*L*L/4)*z)
    Rc = h*(k*k+3*L*L/2)*d
    zero("work cancels leading canonical k squared term", Rc.subs({"L":LB})-Rc.subs({"L":LA})-h*3*(LB*LB-LA*LA)*d/2)
    zero("time pressure canonical signed coefficient", (Rc.subs({"L":LA})/LA-Rc.subs({"L":LB})/LB)/3-h*(k*k/3-(LB-LA)/2)*d)
    zero("time pressure amplitude triangle coefficient", (h*(LB*k+3*LB*LB/2)/LB+h*(LA*k+3*LA*LA/2)/LA)*z/3-h*(2*k/3+(LA+LB)/2)*z)
    wrongR, wrongP = R-3*L*L*c/(2*k), P+L*L*c/(2*k)
    zero("paired geometry omission Ward blind spot", deriv(wrongR)-L*(wrongR-3*wrongP))
    check("paired geometry omission rejected by direct density", bool((wrongR-R).t))
    check("paired geometry omission rejected by direct pressure", bool((wrongP-P).t))
    check("density majorant nonnegative gap for positive k L", F(3)>0)
    check("pressure majorant nonnegative gap for positive k L", F(21,16)>0)
    check("represented Pi is retained, not mathematical pi or three", PI>3 and PI<F(22,7))

    def fixture(kval=F(3,8), weight=F(17,64), dval=F(1,16)):
        wr, wi = Interval(F(-7,64), F(1,32)), Interval(F(3,64), F(5,64))
        return NodeError(kval, weight, dval, wi.shift(dval).scale(F(1)/(2*kval)),
                         Interval(F(-3,128), F(1,64)), wr, wi)
    sample = fixture()
    for kval in (F(1,1<<60), F(1,8), F(3,8), F(2), F(255)):
        n = fixture(kval)
        ts = node_terms(n)
        mu = n.weight*kval*kval/(2*PI*PI)
        rvals, pvals = [], []
        for wr, wi in itertools.product((n.delta_w_re.lo,n.delta_w_re.hi),(n.delta_w_im.lo,n.delta_w_im.hi)):
            u = (n.d+wi)/(2*kval)
            directR = mu*((2*kval*kval+3*LA*LA)*u-kval*wi-LA*wr)/(2*kval)
            directP = mu*((2*kval*kval/3-LA*LA)*u-kval*wi-LA*wr)/(2*kval)
            rvals.append(directR); pvals.append(directP)
        check(f"exact density rectangle extrema k={kval}", ts["anchor_R_W"].shift(ts["anchor_R_c"])==Interval(min(rvals),max(rvals)))
        check(f"exact pressure rectangle extrema k={kval}", ts["anchor_P_W"].shift(ts["anchor_P_c"])==Interval(min(pvals),max(pvals)))
        # Rational unit-circle parametrization tests honest coupled phases.
        for f in (F(-10),F(-1),F(0),F(2),F(19)):
            cs,sn = (1-f*f)/(1+f*f),2*f/(1+f*f)
            for Lval in (LA,F(1,4),LB):
                for wr, wi in itertools.product((n.delta_w_re.lo,n.delta_w_re.hi),(n.delta_w_im.lo,n.delta_w_im.hi)):
                    av, bv = wi/(2*kval), -wr/(2*kval)
                    pp, qq = cs*av-sn*bv, sn*av+cs*bv
                    rr = mu*((kval+3*Lval*Lval/(2*kval))*n.c+3*Lval*Lval*pp/(2*kval)+Lval*qq)
                    pv = mu*((kval/3-Lval*Lval/(2*kval))*n.c-(2*kval/3+Lval*Lval/(2*kval))*pp+Lval*qq)
                    if abs(rr)>ts["R_c_abs_sum"]+ts["R_A_abs_sum"] or abs(pv)>ts["P_c_abs_sum"]+ts["P_A_abs_sum"]:
                        raise RuntimeError("rational phase fixture exceeds uniform bound")
        check(f"rational phase fixtures bounded k={kval}", True)
    for name, call in (
        ("float momentum rejected",lambda:replace(sample,k=0.5)),
        ("nonpositive momentum rejected",lambda:replace(sample,k=F(0))),
        ("negative inherited weight rejected",lambda:replace(sample,weight=F(-1))),
        ("zero inherited weight rejected",lambda:replace(sample,weight=F(0))),
        ("reversed rectangle rejected",lambda:Interval(F(2),F(1))),
        ("canonical projection inconsistent with saved state rejected",lambda:replace(sample,d=F(0))),
        ("wrong target exponent rejected",lambda:from_saved_target(F(1),F(1),(F(0),F(0)),(F(0),F(0)),{},1)),
        ("missing registered capsule rejected",lambda:aggregate_registered({})),
        ("incomplete capsule rejected",lambda:aggregate_capsule("fake",[sample],2,{1:2})),
        ("excess capsule node rejected",lambda:aggregate_capsule("fake",[sample,sample],1,{1:1})),
        ("wrong positional cutoff rejected",lambda:aggregate_capsule("fake",[fixture(F(2))],1,{1:1})),
        ("cutoff equality rejected",lambda:aggregate_capsule("fake",[fixture(F(1))],1,{1:1})),
        ("prefix count ordering rejected",lambda:aggregate_capsule("fake",[],2,{1:2,2:1})),
        ("prefix cutoff ordering rejected",lambda:aggregate_capsule("fake",[],2,{2:1,1:2})),
    ): rejects(name,call)
    def duplicate():
        accumulator=PrefixAccumulator(); accumulator.add(sample); accumulator.add(sample)
    rejects("duplicate momentum rejected",duplicate)
    first = fixture(F(1,8), F(17,64), F(1,16))
    second = fixture(F(3,8), F(19,64), F(-1,16))
    res = aggregate_capsule("fabricated",[first,second],2,{1:2})[0]
    ss=res["sums"]
    check("actual weight sum kept unequal to cutoff",ss["weight_sum"]==F(9,16) and ss["weight_sum"]!=res["K"])
    check("signed anchor interval finite positive sum",ss["anchor_R_interval"]==node_terms(first)["anchor_R_W"].shift(node_terms(first)["anchor_R_c"])+node_terms(second)["anchor_R_W"].shift(node_terms(second)["anchor_R_c"]))
    check("canonical cancellation improves or equals triangle work",ss["work_abs_upper"]<=ss["work_c_abs_sum"]+ss["work_A_abs_sum"])
    check("time-pressure bound distinguished from absolute pressure integral",ss["int_pressure_abs_upper"]!=ss["integral_absolute_pressure_abs_upper"])
    wrong_sign = node_terms(sample)["anchor_R_W"].shift(-node_terms(sample)["anchor_R_c"])
    different("canonical sign mutation rejected",wrong_sign,node_terms(sample)["anchor_R_W"].shift(node_terms(sample)["anchor_R_c"]))
    different("pressure/density substitution mutation rejected",ss["anchor_R_interval"],ss["anchor_P_interval"])
    altered = replace(first,delta_u_im=Interval(F(-1000),F(1000)))
    check("imaginary U retained in norm but does not enter stress",node_terms(altered)["anchor_R_W"]==node_terms(first)["anchor_R_W"] and node_terms(altered)["delta_U_weighted_L1_upper"]>node_terms(first)["delta_U_weighted_L1_upper"])
    targetzero={"U":{"real":[0,0],"imag":[0,0]},"W":{"real":[0,0],"imag":[0,0]}}
    bridge=from_saved_target(F(3,8),F(1,16),(F(7,16),F(-1,8)),(F(5,32),F(3,64)),targetzero)
    check("integer target bridge retains exact saved residual",bridge.d==2*F(3,8)*F(7,16)-F(3,64))
    check("Wronskian defect uses represented epsilon",res["suprema"]["first_order_Wronskian_defect_abs"]==2*EPSILON*res["suprema"]["c_abs"])

    benchmark = {"status":"NOT_REQUESTED"}
    if args.full_fabricated:
        def nodes(grid):
            count = 8192 if grid=="coarse" else 16384
            spacing = 32 if grid=="coarse" else 64
            for index in range(count):
                kval=F(2*index+1,2*spacing)
                # Binary80-scale dyadics and a shared represented-epsilon denominator.
                savedu=(F((index%17)-8,1<<66)/EPSILON,F((index%13)-6,1<<68)/EPSILON)
                savedw=(F((index%11)-5,1<<65)/EPSILON,F((index%19)-9,1<<67)/EPSILON)
                yield from_saved_target(kval,F(1,spacing)+F((index%3)-1,1<<80),savedu,savedw,targetzero)
        start=time.monotonic()
        full=aggregate_registered({key:nodes(key.rsplit("/",1)[1]) for key in CAPSULES})
        duration=time.monotonic()-start
        encoded=json.dumps(encode_exact(full),sort_keys=True,separators=(",",":"))
        check("full fabricated registration has twelve rows",full["row_count"]==12)
        check("full fabricated registration all 49152 nodes",sum(row["node_count"] for row in full["rows"] if row["K"]==256)==49152)
        check("full fabricated registration prefix appearances",sum(row["node_count"] for row in full["rows"])==86016)
        benchmark={"status":"PASS_FABRICATED_ONLY","distinct_nodes":49152,"rows":12,
                   "seconds":round(duration,3),"exact_summary_bytes":len(encoded),
                   "exact_summary_sha256":hashlib.sha256(encoded.encode()).hexdigest()}
    scientific={"result":"PASS_EXACT_DISCRETE_TRANSPORT_FABRICATED_ONLY",
                "checks":checks,"mutation_controls":mutations,
                "retained_array_reads":0,"retained_array_decodes":0,"real_source_callbacks":0,
                "full_integral_pressure_contact_gate":"UNRESOLVED",
                "higher_dimensional_cause":"NOT_ESTABLISHED","novelty":"NOT_ASSESSED"}
    digest=hashlib.sha256(json.dumps(scientific,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    receipt={"scientific":scientific,"scientific_sha256":digest,
             "python_optimization":sys.flags.optimize,"benchmark":benchmark,
             "sources_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(__file__).with_name("discrete_transport.py"))}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"checks":len(checks),"mutations":len(mutations),"scientific_sha256":digest,"benchmark":benchmark},sort_keys=True))

if __name__=="__main__": main()
