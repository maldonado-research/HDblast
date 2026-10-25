#!/usr/bin/env python3
"""Exact golden controls plus independent profile-primitive checks."""
import importlib.util
import json
import math
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ridge_kernel", HERE/"reflected_ridge.py")
k = importlib.util.module_from_spec(spec)
spec.loader.exec_module(k)


def require(ok, label):
    if not ok:
        raise ValueError(label)


def fixture_profile(data, component):
    polynomial = list(map(Q, data["polynomial"]))
    hinges = [tuple(map(Q, row)) for row in data["hinges"]]
    def evaluate(t, order=0):
        t = Q(t)
        value = sum((c*t**(n+order)*Q(math.factorial(n), math.factorial(n+order))
                     for n, c in enumerate(polynomial)), Q(0))
        for c, e, knot in hinges:
            u = e*t-knot
            if u > 0:
                value += c*(e**order)*u**(order+1)/math.factorial(order+1)
        out = [Q(0)]*5
        out[component] = value
        return out
    return evaluate


def main():
    fixtures = json.loads((HERE/"EXACT_FIXTURES.json").read_text())
    comparisons = 0
    for row in fixtures["cases"]:
        component = 0 if row["reflection"] == 1 else 1
        source = fixtures["sources"][row["source"]]
        result = k.reflected_kernel(fixture_profile(source["F"], component),
                                    fixture_profile(source["G"], component),
                                    Q(row["x"]), Q(row["r"]))
        for name in ("value", "dx", "dr"):
            require(result[name][component] == Q(row[name]), "polygon/endpoint golden "+name)
            require(all(v == 0 for j, v in enumerate(result[name]) if j != component),
                    "inactive means")
            comparisons += 5

    nodes = list(map(Q, ("0", "1", "5/4", "3/2", "7/4")))
    urows = [[0]*5, [0]*5, [1,0,2,3,4], [2,0,-1,5,3], [3,0,1,-2,7]]
    drows = [[0]*5, [0]*5, [-1,0,2,3,1], [5,0,2,-1,4], [3,0,1,2,-3]]
    vrows = [[0]*5, [0]*5, [2,1,-1,3,2], [-1,2,3,1,4], [4,-2,1,3,-1]]
    u = k.PiecewiseLinear(nodes, urows)
    du = k.PiecewiseLinear(nodes, drows, tuple(-s for s in k.PARITY))
    v = k.PiecewiseLinear(nodes, vrows)
    # Direct integrate (t-s)^(order-1) times each affine segment.
    primitive_checks = 0
    for profile in (u, du, v):
        for t in nodes + [Q(9,8), Q(11,8), Q(13,8)]:
            for order in (1,2,3):
                direct = [Q(0)]*5
                for a,b,vl,vr in zip(nodes,nodes[1:],profile.values,profile.values[1:]):
                    end = min(t,b)
                    if end <= a:
                        continue
                    for component in range(5):
                        slope = (vr[component]-vl[component])/(b-a)
                        intercept = vl[component]-slope*a
                        for j in range(order):
                            coeff = Q(math.comb(order-1,j))*(t**(order-1-j))*(-1)**j/math.factorial(order-1)
                            direct[component] += coeff*(intercept*(end**(j+1)-a**(j+1))/(j+1)
                                                       +slope*(end**(j+2)-a**(j+2))/(j+2))
                require(tuple(direct) == profile.primitive(t,order), "independent P1 primitive")
                expected_negative = tuple((-1)**order*s*x for s,x in zip(profile.parity,direct))
                require(profile.primitive(-t,order) == expected_negative, "negative primitive orientation")
                primitive_checks += 10
        for t, row in zip(nodes, profile.values):
            require(profile.primitive(t) == row, "shared node value")

    A = [[Q((-1)**(i+j)*(i+1), 19+j) for j in range(5)] for i in range(5)]
    B = [[Q((i+2)*(j+1), 53) for j in range(5)] for i in range(5)]
    D = [[Q(i-j, 71) for j in range(5)] for i in range(5)]
    M = [[Q((i+1)*(j+1), 101) if i!=1 and j!=1 else Q(0)
          for j in range(5)] for i in range(5)]
    source = k.CoupledRidges(A,B,D,u,du,v)
    coupled_checks = 0
    wrong_reflection_detected = False
    points = [(Q(0),Q(0)), (Q(1),Q(0)), (Q(9,8),Q(0)), (Q(5,4),Q(0)),
              (Q(3,2),Q(1,4)), (Q(3,4),Q(3,4)), (Q(1,2),Q(5,4)),
              (Q(1,4),Q(9,8)), (Q(1,2),Q(1,2))]
    for x,r in points:
        a,b = r+x,r-x
        value = k.scale(Q(1,2),k.va(k.va(u.primitive(a),u.primitive(b)),
                                  k.vs(v.primitive(a,1),v.primitive(b,1))))
        dx = k.scale(Q(1,2),k.va(k.vs(du.primitive(a),du.primitive(b)),
                               k.va(v.primitive(a),v.primitive(b))))
        dr = k.scale(Q(1,2),k.va(k.va(du.primitive(a),du.primitive(b)),
                               k.vs(v.primitive(a),v.primitive(b))))
        direct_source = k.va(k.va(k.mv(A,dx),k.mv(B,dr)),k.mv(D,value))
        require(source.physical_source(x,r) == direct_source, "all five coupled ridge terms")
        out = source.feedback(x,r,M)
        if r == 0:
            projected = tuple(z if s==1 else Q(0) for z,s in
                              zip(k.va(u.primitive(x),v.primitive(x,1)),k.PARITY))
            require(out["value"][1] == out["dx"][1] == 0, "Dirichlet shell")
            require(tuple(out["dr"][i] for i in (0,2,3,4)) ==
                    tuple(k.mv(M,projected)[i] for i in (0,2,3,4)), "Robin lift sign")
        if x+r <= 1:
            require(all(z==0 for row in out.values() for z in row), "exact zero wedge")
        if r>0:
            correct = tuple(s*z for s,z in zip(k.PARITY,source.physical_source(x,r)))
            blind = k.va(source.F(-r+x),source.G(-r-x))
            wrong_reflection_detected |= correct != blind
        coupled_checks += 1
    require(wrong_reflection_detected, "noncommuting coefficient reflection control")
    rejected = 0
    for task in (
        lambda: k.PiecewiseLinear([0,1],[[0]*5,[1]*4]),
        lambda: k.PiecewiseLinear([0,0],[[0]*5,[0]*5]),
        lambda: k.PiecewiseLinear([0,1],[[0,1,0,0,0],[0]*5]),
        lambda: source.feedback(Q(3,2),Q(0),M),
        lambda: k.reflected_kernel(lambda t,n:[1]*5,lambda t,n:[0]*5,1,0),
    ):
        try:
            task()
        except ValueError:
            rejected += 1
        else:
            raise ValueError("invalid representation accepted")
    report = {"status":"PASS_EXACT_PYTHON_RIDGE_AND_PROFILE_CONTROLS",
              "independent_fixture_cases":len(fixtures["cases"]),
              "fixture_scalar_equalities":comparisons,
              "profile_primitive_equalities":primitive_checks,
              "coupled_five_mean_cases":coupled_checks,
              "invalid_representations_rejected":rejected,
              "noncommuting_wrong_reflection_rejected":wrong_reflection_detected,
              "scope":"Exact finite kernel/representation controls; not a retained-field error certificate."}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
