#!/usr/bin/env python3
"""Independent exact-rational all-cell audit of coordinate shell sources.

No producer, prior scientific, or other reviewer program is imported.  The
upstream whole-cell H and ell*S certificates remain conditional inputs.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile


class AuditFailure(Exception):
    pass


def require(ok, message):
    if not ok:
        raise AuditFailure(message)


ELL = Q(1, 10**6)
L = Q(3)
M = Q(5, 2) * ELL**2
CONTRACTION = M * L**2 / 2
VALUE_DEN = 2**144
FIELDS = ("y", "yprime", "u", "utau", "uy", "v", "vy", "Sp", "SK", "EN", "dSp", "dSK", "dEN")
SOURCES = ("Sp","SK","EN","dSp","dSK","dEN","Sp_material","SK_material","EN_material")
INPUT_PINS_SHA256 = "2507bbb4f279ec146fcea7342bbd17322f3ded6eb86eba346aeb3bd985765481"
PINS = {
    "INERTIAL_WAVEFORM_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "INERTIAL_WAVEFORM_ENCLOSURE.json": "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8",
    "PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
}
REFERENCE_PINS = {
    "inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md":"8b4d08bf3c4b1704be99e291049a61c902feb39439f8edbda54a5d7d3f9735c8",
    "inputs/reference/PHYSICAL_SHELL_RESPONSE.md":"a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b",
}


def interval(lo, hi=None):
    lo, hi = Q(lo), Q(lo if hi is None else hi)
    require(lo <= hi, "reversed interval")
    return (lo, hi)


def plus(x, y):
    return (x[0] + y[0], x[1] + y[1])


def minus(x, y):
    return (x[0] - y[1], x[1] - y[0])


def times(x, y):
    p = [a*b for a in x for b in y]
    return (min(p), max(p))


def scale(c, x):
    return times((Q(c), Q(c)), x)


def squared(x):
    return (Q(0) if x[0] <= 0 <= x[1] else min(x[0]**2, x[1]**2), max(x[0]**2, x[1]**2))


def inv(x):
    require(x[0] > 0, "inverse expects a strictly positive interval")
    return (1/x[1], 1/x[0])


def hull(xs):
    return (min(x[0] for x in xs), max(x[1] for x in xs))


def qrange(c0, c1, c2, a, b):
    """Closed-cell range of c0+c1*x+c2*x*x, including interior vertex."""
    candidates = [a, b]
    if c2:
        z = -c1/(2*c2)
        if a < z < b:
            candidates.append(z)
    values = [c0 + z*(c1 + z*c2) for z in candidates]
    return (min(values), max(values))


def primitive(p0, p1, h, a, left, right, omit_half=False):
    k = Q(1) if omit_half else Q(1, 2)
    # Completed-cell kernel is x*p0-p1. Current-cell term is h*(x-a)^2/2.
    return qrange(k*h*a*a-p1, p0-2*k*h*a, k*h, left, right)


def box_json(x):
    return {"exact": [str(x[0]), str(x[1])], "display_only": [float(x[0]), float(x[1])]}


def exact_json(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, tuple):
        return [exact_json(x) for x in value]
    if isinstance(value, list):
        return [exact_json(x) for x in value]
    if isinstance(value, dict):
        return {k: exact_json(v) for k, v in value.items()}
    return value


def authenticate(root):
    pinfile=root/"inputs/INPUT_PINS.json"
    require(hashlib.sha256(pinfile.read_bytes()).hexdigest()==INPUT_PINS_SHA256,"immutable input inventory SHA256")
    entries=json.loads(pinfile.read_text())["files"]
    require(len(entries)==18 and len({x["path"] for x in entries})==18,"complete unique input inventory")
    for entry in entries:
        relative=Path(entry["path"])
        require(relative.parts[0]=="inputs" and ".." not in relative.parts,"input path stays within capsule")
        blob=(root/relative).read_bytes()
        require(len(blob)==entry["bytes"] and hashlib.sha256(blob).hexdigest()==entry["sha256"],"full input inventory identity: "+str(relative))
    found = {}
    for name, expected in PINS.items():
        matches = [root / "inputs" / "prior" / name, root / "inputs" / name]
        matches = [p for p in matches if p.is_file()]
        require(len(matches) == 1, "need exactly one pinned input " + name)
        p = matches[0]
        data = p.read_bytes()
        require(hashlib.sha256(data).hexdigest() == expected, "input SHA256 mismatch: " + name)
        found[name] = p
    for relative, expected in REFERENCE_PINS.items():
        p = root/relative
        require(p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == expected, "derivation source SHA256 mismatch: "+relative)
    return found


def load_inputs(root):
    paths = authenticate(root)
    report = json.loads(paths["INERTIAL_WAVEFORM_ENCLOSURE.json"].read_text())
    shell = json.loads(paths["PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json"].read_text())
    require(Q(report["state"]["ell"]) == ELL, "wrong ell")
    den = report["arithmetic"]["value_denominator"]
    grid = report["cell_data"]["time_grid_denominator"]
    cells = []
    with paths["INERTIAL_WAVEFORM_CELLS.csv"].open(newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            a, b = Q(row["left_grid"])/grid, Q(row["right_grid"])/grid
            require(a < b, "non-positive cell width")
            require(a == (cells[-1]["b"] if cells else 0), "cell gap or overlap")
            cells.append({"a": a, "b": b,
                          "H": interval(Q(row["H_lo"])/den, Q(row["H_hi"])/den),
                          "v": interval(Q(row["ell_S_lo"])/den/ELL, Q(row["ell_S_hi"])/den/ELL)})
    require(len(cells) == 4353 and cells[-1]["b"] == L, "wrong continuous 4353-cell cover")
    params = {k: interval(*v["exact"]) for k, v in shell["boxes"].items()}
    # The actual waveform box, slightly outward from the analytic shell box.
    params["mu"] = interval(*report["shell_parameter_intervals"]["mu"])
    true_mu = interval(*shell["boxes"]["mu"]["exact"])
    require(params["mu"][0] <= true_mu[0] <= true_mu[1] <= params["mu"][1], "actual mu does not enclose shell mu")
    require(params["p_b"][1] < 0 and params["a_b"][0] > 0 and params["sigma_second_b"][0] > 0, "shell signs")
    require(CONTRACTION == Q(9, 800000000000) < 1, "global contraction")
    return cells, params, report, paths


def reconstruct(cells, params, control=None):
    c = scale(ELL, inv(params["mu"]))
    if control == "wrong_scaled_force":
        c = inv(params["mu"])
    p0, p1 = [Q(0), Q(0)], [Q(0), Q(0)]
    raw = []
    B_H = Q(0)
    B_H1 = Q(0)
    for cell in cells:
        a, b = cell["a"], cell["b"]
        h = cell["H"]
        bounds = [primitive(p0[k], p1[k], h[k], a, a, b, control == "omit_triangular_half") for k in range(2)]
        I = (bounds[0][0], bounds[1][1])
        F = times(c, I)
        if control == "wrong_negative_mu" and I[1] < 0:
            F = (c[0]*I[0], c[1]*I[1])
        derivs = [(p0[k]+h[k]*(a-a), p0[k]+h[k]*(b-a)) for k in range(2)]
        Fp = times(c, (min(derivs[0]), max(derivs[1])))
        raw.append({**cell, "I": I, "F": F, "Fp": Fp})
        absH = max(abs(h[0]), abs(h[1]))
        B_H += absH*((L-a)**2-(L-b)**2)/2
        B_H1 += absH*(b-a)
        for k in range(2):
            p0[k] += h[k]*(b-a)
            p1[k] += h[k]*(b*b-a*a)/2
    signed_Y0 = max(max(abs(x) for x in z["F"]) for z in raw)
    coarse_Y0 = c[1]*B_H
    By = signed_Y0/(1-CONTRACTION)
    E = CONTRACTION*By
    Ep = M*L*By
    if control == "omit_resolvent":
        By, E, Ep = signed_Y0, CONTRACTION*signed_Y0, M*L*signed_Y0
    H2 = params["H_squared"]
    rho = params["rho_b"]
    p, a, sig2 = params["p_b"], params["a_b"], params["sigma_second_b"]
    rows = []
    for z in raw:
        left, right = z["a"], z["b"]
        y = plus(z["F"], (-E, E))
        yp = plus(z["Fp"], (-Ep, Ep))
        clock = (1-(ELL*right)**2, Q(1))
        # Elementary monotonicity and tanh(t) >= t-t^3/3 for t>=0.
        tanh = (Q(0), ELL*right)
        u = times(times(H2, clock), y)
        velocity_core = minus(yp, scale(Q(3, 2)*ELL, times(tanh, y)))
        du = times(times(times(H2, clock), scale(1/ELL, inv(rho))), velocity_core)
        v = z["v"]
        if control == "wrong_metric_radius":
            u = times(times(inv(rho), clock), y)
        coefficients=source_intervals(u,v,params)
        uy = minus(scale(-2,times(a,u)), scale(Q(1,3),times(p,v)))
        vy = minus(scale(-2,times(p,u)), scale(Q(1,2),times(sig2,v)))
        row = {**z, "y": y, "yprime": yp, "clock": clock, "tanh": tanh, "u": u, "utau": du,
               "uy":uy, "vy":vy, **coefficients}
        for name in ("Sp","SK"):
            source=row[name]
            for degree, coefficient in ((0,Q(1,6)),(2,Q(1,2)),(4,Q(1))):
                row[f"{name}_l{degree}"] = scale(coefficient, source)
        rows.append(row)
    return rows, {"c":c, "B_H": B_H, "B_H1": B_H1, "coarse_Y0": coarse_Y0,
                  "coarse_By": coarse_Y0/(1-CONTRACTION), "coarse_Bu": H2[1]*coarse_Y0/(1-CONTRACTION),
                  "signed_Y0": signed_Y0, "By": By, "E": E, "Ep": Ep,
                  "q": CONTRACTION, "M": M}


def source_intervals(u,v,params):
    a,p,sig2=params["a_b"],params["p_b"],params["sigma_second_b"]
    u2,v2,uv=squared(u),squared(v),times(u,v)
    result={"dSp":plus(scale(-4,times(p,u2)),scale(2,times(sig2,uv))),
            "dSK":scale(-12,times(a,u2)),
            "dEN":plus(scale(-16,times(a,u2)),scale(Q(-16,3),times(p,uv))),
            "Sp_material":scale(-2,v2),
            "SK_material":scale(Q(1,6),times(sig2,v2)),
            "EN_material":scale(Q(2,3),times(sig2,v2))}
    for name in ("Sp","SK","EN"):
        result[name]=plus(result["d"+name],result[name+"_material"])
    return result


@lru_cache(maxsize=100000)
def rounded(x):
    lo=x[0].numerator*VALUE_DEN//x[0].denominator
    hi=-((-x[1].numerator*VALUE_DEN)//x[1].denominator)
    return (Q(lo,VALUE_DEN),Q(hi,VALUE_DEN))


def rational_pi():
    """Machin identity, alternating atan series with explicit rational tail."""
    def atan_recip(k, n):
        s = sum((Q((-1)**j, (2*j+1)*k**(2*j+1)) for j in range(n)), Q(0))
        nxt = Q((-1)**n, (2*n+1)*k**(2*n+1))
        return (min(s,s+nxt), max(s,s+nxt))
    return rounded(minus(scale(16, atan_recip(5,64)), scale(4, atan_recip(239,16))))


def sphere_moment(a,b):
    # Squared coordinates on S^3 have Dirichlet parameters (1/2,1/2,1/2,1/2).
    numerator=Q(1)
    for count in (a,b):
        for j in range(count):
            numerator*=Q(1,2)+j
    for j in range(a+b):
        numerator/=2+j
    return numerator


def angular_checks():
    e_y4=sum(((-1)**j*comb(4,j)*sphere_moment(4-j,j) for j in range(5)),Q(0))
    e_z2sq=sphere_moment(2,0)+2*sphere_moment(1,1)+sphere_moment(0,2)-2*sphere_moment(1,0)+Q(1,4)
    e_h4sq=e_y4-Q(1,36)-e_z2sq/4
    require(e_y4==Q(3,40) and e_z2sq==Q(1,12) and e_h4sq==Q(19,720),"angular harmonic normalization")
    return {"Y4_mean":e_y4,"Z2_squared_mean":e_z2sq,"H4_squared_mean":e_h4sq,
            "Y_squared_projection":(Q(1,6),Q(1,2),Q(1)),
            "constant_angular_factor_over_pi_squared":Q(1,18),
            "degree2_angular_factor_over_pi_squared":Q(1,24),
            "degree4_angular_factor_over_pi_squared":Q(19,360)}


def bracket_result(cells, params, constants, bracket):
    left, right = bracket
    p0,p1 = [Q(0),Q(0)],[Q(0),Q(0)]
    pieces = []
    for cell in cells:
        a,b,h = cell["a"],cell["b"],cell["H"]
        if max(left,a) <= min(right,b):
            ranges = [primitive(p0[k],p1[k],h[k],a,max(left,a),min(right,b)) for k in range(2)]
            pieces.append((ranges[0][0],ranges[1][1]))
        for k in range(2):
            p0[k] += h[k]*(b-a)
            p1[k] += h[k]*(b*b-a*a)/2
    I = hull(pieces)
    F = times(constants["c"],I)
    E = constants["E"]
    y = plus(F,(-E,E))
    u = times(times(params["H_squared"],(1-(ELL*right)**2,Q(1))),y)
    return {"bracket":bracket,"I":I,"F":F,"y":y,"u":u,
            "Sp_at_v_zero":scale(-4,times(params["p_b"],squared(u))),
            "SK_at_v_zero":scale(-12,times(params["a_b"],squared(u)))}


def norms(rows, params):
    pi = rational_pi()
    # rho^4 is positive. For z=t^2<=9e-12, cosh^3(t)<=1+2*z.
    prefactor = times(scale(Q(3,20)*ELL, squared(pi)), squared(squared(params["rho_b"])))
    prefactor1=scale(Q(20,9),prefactor)
    prefactor4=scale(Q(19,54),prefactor)
    result = {}
    for name in SOURCES:
        total,absolute,signed = (Q(0),Q(0)),(Q(0),Q(0)),(Q(0),Q(0))
        for z in rows:
            weight = (Q(1), 1+2*(ELL*z["b"])**2)
            value=rounded(z[name])
            absv=(Q(0) if value[0]<=0<=value[1] else min(map(abs,value)),max(map(abs,value)))
            width=z["b"]-z["a"]
            total = plus(total, scale(width, times(weight, squared(value))))
            absolute=plus(absolute,scale(width,times(weight,absv)))
            signed=plus(signed,scale(width,times(weight,value)))
        fourth=name in ("EN","dEN","EN_material")
        result[name]={"angular_field":"H4" if fourth else "Y^2","weighted_x_squared_integral":total,
                      "squared_L2":times(prefactor4 if fourth else prefactor,total)}
        if not fourth:
            result[name].update(L1=times(prefactor1,absolute),signed_integral=times(prefactor1,signed))
    return {"pi":pi,"L2_prefactor":prefactor,"L1_prefactor":prefactor1,"EN_L2_prefactor":prefactor4,"norms":result}


def self_checks(rows, constants, params, control=None):
    # These are independently stated identities/endpoint witnesses, not asserts.
    require(ELL**2*(1/ELL) == ELL, "scaled force coefficient identity")
    require(constants["c"][0]*params["mu"][1] == ELL, "wrong_scaled_force: ell/mu factor missing")
    expected_coarse = Q("203118884062410186041920423796832536543834836817363247669963296595318905770391/2349828472120094682435900240404606958989727491460798474881531730796345294848000000000")
    require(constants["coarse_Bu"] == expected_coarse, "coarse metric bound or ell scaling mismatch")
    test = primitive(Q(0),Q(0),Q(2),Q(0),Q(0),Q(1),control == "omit_triangular_half")
    require(test == (Q(0),Q(1)), "omit_triangular_half: integral of 2*(x-s) is x^2")
    # This fixture changes sign; assigning mu endpoints by position loses enclosure.
    c = scale(ELL, inv(params["mu"]))
    negative = times(c, (Q(-2), Q(-1)))
    attempted = (c[0]*-2, c[1]*-1) if control == "wrong_negative_mu" else negative
    require(attempted[0] <= c[1]*-2 and c[0]*-1 <= attempted[1], "wrong_negative_mu: signed endpoint witness excluded")
    require(constants["By"] >= constants["signed_Y0"] + CONTRACTION*constants["By"], "omit_resolvent: contraction residual not closed")
    zmax=(ELL*L)**2
    require(zmax<Q(1,100) and 4-18*Q(1,100)-2*Q(1,100)**3>0, "cosh-cubed measure upper bound")
    attempted_angular=Q(3,20) if control=="wrong_angular_normalization" else Q(3,40)
    moment=sum(((-1)**j*comb(4,j)*sphere_moment(4-j,j) for j in range(5)),Q(0))
    require(attempted_angular==moment,"wrong_angular_normalization: sphere Y^4 mean is 3/40")
    clock_upper=Q(1) if control=="omit_measure_clock" else 1+2*zmax
    require(clock_upper>=1+Q(3,2)*zmax,"omit_measure_clock: positive Taylor term of cosh^3 excluded")
    for z in rows:
        require(z["u"] == times(times(params["H_squared"],z["clock"]),z["y"]), "wrong_metric_radius: u=rho^-2*C*y")
        for key in ("F","Fp")+FIELDS:
            require(z[key][0] <= z[key][1], "reversed mathematical cell " + key)
    require(max(abs(x) for z in rows for x in z["u"]) < Q("1.12e-8"), "sharper signed global u bound")


def execute_control(name,params):
    """Actually compute a false quantity, then check its mathematical witness.

    These local exact witnesses isolate a particular formula rather than
    repeating unrelated all-cell source arithmetic in each deliberate failure.
    """
    if name=="wrong_scaled_force":
        wrong=inv(params["mu"])
        require(wrong[0]*params["mu"][1]==ELL,"wrong_scaled_force: ell/mu factor missing")
    elif name=="omit_triangular_half":
        wrong=primitive(Q(0),Q(0),Q(2),Q(0),Q(0),Q(1),True)
        require(wrong==(Q(0),Q(1)),"omit_triangular_half: integral of 2*(x-s) is x^2")
    elif name=="wrong_negative_mu":
        c=scale(ELL,inv(params["mu"]))
        wrong=(c[0]*-2,c[1]*-1)
        require(wrong[0]<=c[1]*-2 and c[0]*-1<=wrong[1],"wrong_negative_mu: signed endpoint witness excluded")
    elif name=="omit_resolvent":
        Y0=Q(1)
        wrong=Y0
        require(wrong>=Y0+CONTRACTION*wrong,"omit_resolvent: contraction residual not closed")
    elif name=="wrong_metric_radius":
        positive_y=(Q(1),Q(1))
        wrong=times(inv(params["rho_b"]),positive_y)
        require(wrong==times(params["H_squared"],positive_y),"wrong_metric_radius: u=rho^-2*C*y")
    elif name=="wrong_angular_normalization":
        wrong=Q(3,20)
        actual=sum(((-1)**j*comb(4,j)*sphere_moment(4-j,j) for j in range(5)),Q(0))
        require(wrong==actual,"wrong_angular_normalization: sphere Y^4 mean is 3/40")
    elif name=="omit_measure_clock":
        wrong=Q(1)
        require(wrong>=1+Q(3,2)*(ELL*L)**2,"omit_measure_clock: positive Taylor term of cosh^3 excluded")
    else:
        raise AuditFailure("unknown mathematical control")
    raise AuditFailure("deliberate wrong quantity unexpectedly passed")


def summarize(rows, constants, params, cells, prior):
    result = {"schema":"hdblast-sol-coordinate-shell-audit-v1", "status":"PASS",
              "cell_count":len(rows), "arithmetic":"fractions.Fraction; no assert checks",
              "finite_groups":["input identity and coverage", "signed Volterra primitive and resolvent", "clock/radius derivative reconstruction", "complete coefficient and harmonic rows", "whole-cell source norms"],
              "constants": constants, "norms": norms(rows,params),"angular":angular_checks()}
    result["global_boxes"] = {k:hull([rounded(z[k]) for z in rows]) for k in FIELDS}
    result["sign_counts"] = {k:{"strict_positive":sum(z[k][0]>0 for z in rows), "strict_negative":sum(z[k][1]<0 for z in rows), "contains_zero":sum(z[k][0]<=0<=z[k][1] for z in rows)} for k in ("Sp","SK")}
    tail = [z for z in rows if z["a"]>=2]
    result["tail_boxes"]={k:hull([z[k] for z in tail]) for k in ("u","v","Sp","SK")}
    require(result["tail_boxes"]["Sp"][1] < 0 and result["tail_boxes"]["SK"][0] > 0, "uniform tail signs")
    central=bracket_result(cells,params,constants,(Q("1.49998"),Q("1.49999")))
    require(central["u"][0]>0 and central["Sp_at_v_zero"][0]>0 and central["SK_at_v_zero"][1]<0,"central coordinate source signs")
    witness=prior["sign_witnesses"]
    require(Q(witness["central_left"]["x"])==central["bracket"][0] and Q(witness["central_right"]["x"])==central["bracket"][1],"central witness clocks")
    require(Q(witness["central_left"]["H"][0])>0>Q(witness["central_right"]["H"][1]),"central opposite-sign actual-H witnesses")
    result["central"] = central
    result["scope"] = "Conditional coordinate coefficients in retained longitudinal gauge and background proper-clock slab; no evaluated U2, invariant response sign, upstream PDE error-certificate reproduction, or external peer review."
    return result


def expected_primary(rows,constants,params,cells,prior,normdata):
    C=constants
    central=[]
    for bracket in prior["zero_existence_brackets"]:
        z=bracket_result(cells,params,C,tuple(map(Q,bracket)))
        source=source_intervals(z["u"],(Q(0),Q(0)),params)
        central.append({k:z[k] for k in ("bracket","F","y","u")})
        central[-1].update({name+"_at_scalar_zero":source[name] for name in ("Sp","SK","EN")})
        require(z["u"][0]>0 and source["Sp"][0]>0 and source["SK"][1]<0,"central bracket source identities")
    Usigned=params["H_squared"][1]*C["By"]
    gain=interval(*prior["physical_gain_interval"])
    Htail=interval(*prior["tail_theorem"]["actual_H_interval"])
    require(all(Htail[0]<=z["H"][0]<=z["H"][1]<=Htail[1] for z in cells if z["a"]>=2),"retained tail theorem contains every actual-H tail cell")
    vtail=scale(1/ELL,times(times(gain,(1-(ELL*L)**2,Q(1))),Htail))
    tail=source_intervals((-Usigned,Usigned),vtail,params)
    a,p,sigma=params["a_b"],params["p_b"],params["sigma_second_b"]
    off=4*max(map(abs,p))
    determinant=minus(scale(24,times(a,sigma)),scale(16,squared(p)))
    lower=min(6*a[0]-off,4*sigma[0]-off)
    require(lower>0 and determinant[0]>0,"registered Robin symmetrizer positivity")
    normvalues=normdata["norms"]
    ratios={name:normvalues["d"+name]["squared_L2"][1]/normvalues[name+"_material"]["squared_L2"][0] for name in ("Sp","SK","EN")}
    expected={"ell":ELL,"L":L,"mu":params["mu"],"rho":params["rho_b"],"r2":params["H_squared"],
              "a":a,"p":p,"sigma2":sigma,"M":M,"q":CONTRACTION,
              "Fsup":C["signed_Y0"],"Ybound":C["By"],"E":C["E"],"Eprime":C["Ep"],
              "BH":C["B_H"],"Ycoarse":C["coarse_By"],"Ucoarse":C["coarse_Bu"],"Usigned":Usigned,
              "global_envelopes":{k:hull([rounded(z[k]) for z in rows]) for k in FIELDS},
              "central":central,"tail":{"domain":(Q(2),Q(3)),"v":vtail,"coefficients":tail},
              "metric_correction_squared_norm_ratio_upper":ratios,
              "symmetrizer":{"weights":[3,8],"gershgorin_positive_lower":lower,"determinant_interval":determinant},
              "harmonics":{"basis":["1","Z2","H4"],"source_coefficient_factors":[Q(1,6),Q(1,2),Q(1)],
                           "squared_L2_power_shares":[Q(10,27),Q(5,18),Q(19,54)],
                           "EN_scope":"Retained degree-four scalar Robin row only; not a degree-zero master reduction."}}
    expected.update(normdata)
    return exact_json(expected)


def compare_tree(actual,expected,path):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and set(actual)==set(expected),"JSON exact key coverage: "+path)
        for k,v in expected.items():
            compare_tree(actual[k],v,path+"/"+k)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),"JSON list coverage: "+path)
        for i,(a,b) in enumerate(zip(actual,expected)):
            compare_tree(a,b,path+"/"+str(i))
    else:
        require(type(actual)==type(expected) and actual==expected,"JSON mathematical value mismatch: "+path)


def float_tree(value):
    if isinstance(value,str):
        try:
            return float(Q(value))
        except ValueError:
            return value
    if isinstance(value,dict):
        return {k:float_tree(v) for k,v in value.items()}
    if isinstance(value,list):
        return [float_tree(v) for v in value]
    return value


def verify_primary(root,rows,expected):
    path=root/"computed/COORDINATE_SHELL_SOURCE_CELLS.csv"
    data=path.read_bytes()
    fields=["left_grid","right_grid"]+[f"{name}_{side}" for name in FIELDS for side in ("lo","hi")]
    with path.open(newline="") as f:
        reader=csv.DictReader(f)
        require(reader.fieldnames==fields,"primary CSV field coverage")
        actual=list(reader)
    require(len(actual)==len(rows)==4353,"primary CSV full-cell coverage")
    count=0
    for i,(a,b) in enumerate(zip(actual,rows)):
        require(a["left_grid"]==str(int(b["a"]*4096)) and a["right_grid"]==str(int(b["b"]*4096)),"primary CSV grid: "+str(i))
        for key in FIELDS:
            roundedvalue=rounded(b[key])
            for side,value in zip(("lo","hi"),roundedvalue):
                require(a[key+"_"+side]==str(int(value*VALUE_DEN)),f"primary all-cell endpoint mismatch: row {i} {key}_{side}")
                count+=1
    reportpath=root/"computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json"
    report=json.loads(reportpath.read_text())
    require(set(report)=={"schema","status","runtime","producer_sha256","input_pins_sha256","cell_data","exact","display_only","scope"},"primary report full top-level schema")
    require(report["schema"]=="hdblast-coordinate-shell-sources-enclosure-v1" and report["status"]=="PASS","primary report status")
    require(report["producer_sha256"]==hashlib.sha256((root/"compute_coordinate_sources.py").read_bytes()).hexdigest(),"producer source hash")
    require(report["input_pins_sha256"]==INPUT_PINS_SHA256,"primary immutable input pins hash")
    require(set(report["runtime"])=={"python","optimization"} and isinstance(report["runtime"]["python"],str) and report["runtime"]["optimization"] in (0,1,2),"runtime provenance")
    compare_tree(report["cell_data"],{"path":path.name,"rows":4353,"time_grid_denominator":4096,"value_denominator":VALUE_DEN,"fields":list(FIELDS),"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()},"/cell_data")
    compare_tree(report["exact"],expected,"/exact")
    display={k:expected[k] for k in ("global_envelopes","Usigned","Ucoarse","metric_correction_squared_norm_ratio_upper","central","tail","symmetrizer")}
    display["norms"]={k:v["squared_L2"] for k,v in expected["norms"].items()}
    compare_tree(report["display_only"],float_tree(display),"/display_only")
    scope={"conditional_on_registered_analytic_inputs":True,
           "first_representative":"Retained longitudinal, zero first shell graph and fixed intrinsic identification.",
           "normal_convention":"One-copy outward coordinate physical rows; not the outward-normal sum material scalar row.",
           "taylor_convention":"Derivative coefficients U=U0+epsilon U1+epsilon^2 U2/2.",
           "metric_and_scalar_sources":"Sp and SK multiply Y^2 in the original coordinate scalar and isotropic Israel rows.",
           "EN_scope":"Retained degree-four reduced Robin forcing only.",
           "norms":"Scalar amplitude norms in background proper shell measure, not an Einstein energy.",
           "full_U2_evaluated":False,"nonlinear_remainder_established":False,
           "coordinate_source_sign_is_gauge_invariant":False,"cosmological_evidence":False}
    compare_tree(report["scope"],scope,"/scope")
    return {"status":"PASS","cells_authenticated":4353,"scientific_endpoints_authenticated":count,
            "scientific_report_subtree":"Every exact key and leaf; every display leaf and scope claim",
            "csv_sha256":hashlib.sha256(data).hexdigest(),"report_sha256":hashlib.sha256(reportpath.read_bytes()).hexdigest()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-root",type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument("--output",type=Path)
    ap.add_argument("--reconstruct-only",action="store_true")
    ap.add_argument("--inject",choices=["wrong_scaled_force","omit_triangular_half","wrong_negative_mu","omit_resolvent","wrong_metric_radius","wrong_angular_normalization","omit_measure_clock"])
    args = ap.parse_args()
    try:
        cells, params, prior, paths = load_inputs(args.dataset_root.resolve())
        if args.inject:
            execute_control(args.inject,params)
        rows, constants = reconstruct(cells,params,args.inject)
        self_checks(rows,constants,params,args.inject)
        result = summarize(rows,constants,params,cells,prior)
        if not args.reconstruct_only:
            expected=expected_primary(rows,constants,params,cells,prior,result["norms"])
            result["primary_authentication"]=verify_primary(args.dataset_root.resolve(),rows,expected)
        result["input_pins"] = PINS
        result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result["optimized"] = not __debug__
        out = json.dumps(exact_json(result),indent=2)+"\n"
        if args.output:
            args.output.write_text(out)
        else:
            print(out,end="")
    except AuditFailure as e:
        print("FAIL: "+str(e),file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
