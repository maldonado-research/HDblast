#!/usr/bin/env python3
"""Export independently computed SymPy1.14.0 expectations for stdlib replay."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import sympy as S


def run(repo):
    x, r = S.symbols("x r")
    ell = S.Rational(2, 7)
    rho = 3 + r / 5 + r**2 / 11 + r**3 / 13
    a = -S.diff(rho, r) / (ell * rho**2)
    amp = (S.Rational(3, 7) + 2*x - 3*r + x**2/5 + 7*x*r/11
           - 2*r**2/13 + x**3/17 - 5*x**2*r/19 + 3*x*r**2/23 + 2*r**3/29)

    def taylor(expression):
        return {f"{i},{k}":str(S.diff(expression, x, i, r, k).subs({x:0, r:0}) / (math.factorial(i)*math.factorial(k))) for i in range(4) for k in range(4-i)}

    physical = {"": amp, "y": -S.diff(amp, r)/(ell*rho), "t": S.diff(amp, x)/ell}
    physical["yy"] = -S.diff(physical["y"], r)/(ell*rho)
    physical["yt"] = S.diff(physical["y"], x)/ell
    physical["tt"] = S.diff(physical["t"], x)/ell
    expected = [{"suffix":suffix, "source_r_order":order, "exact":str(S.diff(expression,r,order).subs({x:0,r:0}))} for suffix,expression in physical.items() for order in (0,1)]
    source_path = repo / "research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json"
    raw = source_path.read_bytes()
    terms = json.loads(raw)["maps"]["metric_forcing"]["ab"][5]
    names = sorted({name for term in terms for name in term["powers"]})
    symbols = {name:S.Symbol(name) for name in names}
    expression = S.Add(*(S.Rational(term["coefficient"])*S.prod(symbols[name]**power for name,power in term["powers"].items()) for term in terms))
    gradient = {name:S.diff(expression,symbol) for name,symbol in symbols.items()}
    cases = []
    for case in range(3):
        values = {name:S.Rational(2*i+3+case,i+2) for i,name in enumerate(names)}
        rates = {name:S.Rational((-1)**(i+case)*(i+1),i+case+3) for i,name in enumerate(names)}
        substitution = {symbols[name]:value for name,value in values.items()}
        fv = expression.subs(substitution)
        fr = sum(gradient[name].subs(substitution)*rates[name] for name in names)
        cases.append({"case":case, "values":{name:str(value) for name,value in values.items()}, "rates":{name:str(value) for name,value in rates.items()}, "F":str(fv), "Fr":str(fr)})
    return {
        "schema":"INDEPENDENT_TT_FR_EXACT_FIXTURE_V1",
        "generator_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version":S.__version__,
        "source_map_sha256":hashlib.sha256(raw).hexdigest(),
        "synthetic_ell":str(ell),
        "synthetic_normalized_Taylor_coefficients":{name:taylor(expression) for name,expression in (("amp",amp),("rho",rho),("a",a))},
        "synthetic_physical_expected":expected,
        "independent_full_source_symbolic_gradient_cases":cases,
        "provenance":"Source expression substituted and differentiated by SymPy, independently of production Dual arithmetic. Physical derivatives computed by repeated -1/(ell*rho)*partial_r and partial_x/ell, independently of production physical_expression.",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError("FRESH_OUTPUT_REQUIRED")
    args.output.write_text(json.dumps(run(args.repo_root.resolve()),indent=2,sort_keys=True)+"\n")
    print(hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__=="__main__":
    main()
