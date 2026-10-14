#!/usr/bin/env python3
"""Global exact S3 scalar/TT decomposition and coexact-selection checker.

The x=sin(theta)^2 Hopf chart gives rational metric coefficients. This
program derives Christoffels directly and does not import inherited code.
"""
import argparse
import json
from pathlib import Path
import sys

import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--wrong-formula", choices=["none", "wrong-l4-coefficient", "omit-tensor", "wrong-vector-divergence"], default="none")
    args = parser.parse_args()
    x, phi, psi = s.symbols("x phi psi", real=True)
    coords = [x, phi, psi]
    metric = s.diag(1/(4*x*(1-x)), x, 1-x)
    inverse = metric.inv()
    Gamma = [[[s.simplify(sum(inverse[i, k]*(s.diff(metric[k,j],coords[l])+s.diff(metric[k,l],coords[j])-s.diff(metric[j,l],coords[k])) for k in range(3))/2)
              for l in range(3)] for j in range(3)] for i in range(3)]
    def clean(v):
        return s.factor(s.trigsimp(s.expand_trig(v)))
    def hessian(F):
        return s.Matrix(3, 3, lambda i,j: clean(s.diff(F,coords[i],coords[j])-sum(Gamma[k][i][j]*s.diff(F,coords[k]) for k in range(3))))
    def trace(A):
        return clean(sum(inverse[i,j]*A[i,j] for i in range(3) for j in range(3)))
    def TF(A):
        return A-metric*trace(A)/3
    def div(A):
        return [clean(sum(inverse[i,k]*(s.diff(A[i,j],coords[k])-sum(Gamma[l][k][i]*A[l,j]+Gamma[l][k][j]*A[i,l] for l in range(3))) for i in range(3) for k in range(3))) for j in range(3)]
    Y = x*s.cos(2*phi)
    F2 = x-s.Rational(1,2)
    Q = Y**2-x/2+s.Rational(1,12)
    dY = s.Matrix([s.diff(Y,z) for z in coords])
    stress = dY*dY.T
    YY = Y*hessian(Y)
    T = s.diag(1/(24*x*(1-x)), x*(s.Rational(1,6)-x/2), (1-x)*(s.Rational(1,6)-(1-x)/2))
    checks = []
    def check(name, values):
        entries = list(values) if isinstance(values, (list, tuple, s.MatrixBase)) else [values]
        for i,v in enumerate(entries):
            residual = clean(v)
            if residual != 0:
                raise ArithmeticError(f"{name}[{i}] failed: {residual}")
            checks.append(f"{name}[{i}]")
    check("Y-eigenvalue", trace(hessian(Y))+8*Y)
    check("Q-eigenvalue", trace(hessian(Q))+24*Q)
    check("F2-eigenvalue", trace(hessian(F2))+8*F2)
    check("Y-square-decomposition", Y**2-Q-F2/2-s.Rational(1,6))
    l4coef = s.Rational(1,2) if args.wrong_formula == "wrong-l4-coefficient" else s.Rational(1,3)
    tcoef = 0 if args.wrong_formula == "omit-tensor" else -s.Rational(4,3)
    check("stress-SVT-decomposition", TF(stress)-TF(hessian(F2))/2-l4coef*TF(hessian(Q))-tcoef*T)
    check("Y-Hessian-product", YY-hessian(Y**2)/2+stress)
    check("tensor-trace", trace(T))
    check("tensor-divergence", div(T))
    gradient_norm = trace(stress)
    divcoef = s.Rational(1,3) if args.wrong_formula == "wrong-vector-divergence" else s.Rational(1,6)
    expected = [s.diff(-4*Y**2+divcoef*gradient_norm,z) for z in coords]
    check("stress-TF-gradient-divergence", [a-b for a,b in zip(div(TF(stress)), expected)])
    u, v = s.symbols("u v", real=True)
    complete = 24*u**2*TF(YY)+(12*u**2-2*v**2)*TF(stress)
    decomposed = (-v**2*TF(hessian(F2))+(8*u**2-s.Rational(2,3)*v**2)*TF(hessian(Q))
                  +(16*u**2+s.Rational(8,3)*v**2)*T)
    check("complete-angular-TF-source", complete-decomposed)
    receipt = {"status":"PASS", "method":"exact symbolic direct Hopf Christoffels on x=sin(theta)^2 chart", "check_count":len(checks), "checks":checks,
               "python_optimized":not __debug__, "sympy_version":s.__version__, "wrong_formula":args.wrong_formula}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError) as error:
        print(str(error),file=sys.stderr)
        raise SystemExit(1)
