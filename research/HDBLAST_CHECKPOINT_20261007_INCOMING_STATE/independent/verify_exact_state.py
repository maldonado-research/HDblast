#!/usr/bin/env python3
"""Independent exact-symbolic checks for conditional incoming-state propagation.

This program needs SymPy 1.14.0 and writes no repository or remote state.
It deliberately uses explicit exceptions rather than Python assertions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import sympy as s


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("Fresh output required; preserve existing receipts")
    checks: list[str] = []
    mutations: list[dict[str, str]] = []

    def zero(name: str, expression: s.Expr) -> None:
        residual = s.trigsimp(s.simplify(expression))
        if residual != 0:
            raise RuntimeError(f"{name}: residual {residual}")
        checks.append(name)

    def reject(name: str, expression: s.Expr, witness: dict | None = None) -> None:
        residual = s.trigsimp(s.simplify(expression))
        if witness is not None:
            residual = s.trigsimp(s.simplify(residual.subs(witness)))
        if residual == 0 or residual.is_zero is True:
            raise RuntimeError(f"Mutation escaped rejection: {name}")
        if residual.is_zero is not False and residual.equals(0) is not False:
            raise RuntimeError(f"Mutation lacks a proven nonzero residual: {name}")
        mutations.append({"name": name, "nonzero_residual": str(residual)})

    k, L, K, Pi = s.symbols("k L K Pi", positive=True)
    t, a = s.symbols("t a", real=True)
    x, y, v, z, gr, gi = s.symbols("x y v z gr gi", real=True)
    c, ar, ai, ur, ui, wr, wi = s.symbols(
        "c ar ai ur ui wr wi", real=True
    )
    X, Y = s.symbols("X Y", real=True)
    chi, sigma = s.symbols("chi sigma", nonnegative=True)
    measure = k**2 / (2 * Pi**2)
    R = ((2*k**2 + 3*L**2)*x - k*z - L*v) / (2*k)
    P = ((2*k**2/s.Integer(3) - L**2)*x - k*z - L*v) / (2*k)

    def linear_derivative(expr: s.Expr, forced: bool = False) -> s.Expr:
        return (s.diff(expr, x)*v + s.diff(expr, y)*z
                + s.diff(expr, v)*(-2*k*z - (gr if forced else 0))
                + s.diff(expr, z)*(2*k*v - (gi if forced else 0))
                + s.diff(expr, L)*L**2)

    source = (k*gi + L*gr)/(2*k)
    zero("direct forced continuity operator", linear_derivative(R, True)
         - L*(R - 3*P) - source)
    zero("homogeneous continuity operator", linear_derivative(R)
         - L*(R - 3*P))
    zero("pressure primitive with L=-1/t",
         (R + t*linear_derivative(R) - 3*P).subs(L, -1/t))
    zero("forced pressure primitive with L=-1/t",
         (R + t*linear_derivative(R, True) - 3*P - t*source)
         .subs(L, -1/t))

    theta = 2*k*(t-a)
    harmonic_x = c + ar*s.cos(theta) - ai*s.sin(theta)
    harmonic_y = ui - ai + ar*s.sin(theta) + ai*s.cos(theta)
    harmonic_v = -2*k*(ar*s.sin(theta) + ai*s.cos(theta))
    harmonic_z = 2*k*(ar*s.cos(theta) - ai*s.sin(theta))
    zero("homogeneous U real evolution", s.diff(harmonic_x, t)-harmonic_v)
    zero("homogeneous U imaginary evolution", s.diff(harmonic_y, t)-harmonic_z)
    zero("homogeneous W real evolution", s.diff(harmonic_v, t)+2*k*harmonic_z)
    zero("homogeneous W imaginary evolution", s.diff(harmonic_z, t)-2*k*harmonic_v)
    incoming = {c: ur-wi/(2*k), ar: wi/(2*k), ai: -wr/(2*k)}
    for name, solution, value in (
        ("U real", harmonic_x, ur), ("U imaginary", harmonic_y, ui),
        ("W real", harmonic_v, wr), ("W imaginary", harmonic_z, wi),
    ):
        zero(f"canonical incoming {name}", solution.subs(t,a).subs(incoming)-value)
    zero("normalization defect invariant", linear_derivative(2*x-z/k))
    zero("canonical invariant equals c", harmonic_x-harmonic_z/(2*k)-c)

    eps = s.symbols("eps", real=True)
    # q=q0(1+eps U), q0=e^(-ikt)/sqrt(2k), U'=W.
    factor = 1+eps*(x+s.I*y)
    factor_prime = eps*(v+s.I*z)
    qnorm = s.expand(s.I*(s.conjugate(factor)*(-s.I*k*factor+factor_prime)
                         - factor*(s.I*k*s.conjugate(factor)
                                   + s.conjugate(factor_prime)))/(2*k))
    zero("exact oscillator Wronskian expression",
         qnorm - (s.expand_complex(factor*s.conjugate(factor))
                  - s.im(s.conjugate(factor)*factor_prime)/k))
    zero("first order Wronskian normalization defect",
         s.diff(qnorm,eps).subs(eps,0) - (2*x-z/k))
    zero("Wronskian quadratic remainder",
         qnorm - 1 - eps*(2*x-z/k)
         - eps**2*(x*x+y*y-(x*z-y*v)/k))

    canonical = {x:c+X, z:2*k*X, v:-2*k*Y}
    Rc = k+3*L**2/(2*k)
    Ra = 3*L**2/(2*k)*X + L*Y
    Pc = k/3-L**2/(2*k)
    Pa = -(2*k/3+L**2/(2*k))*X + L*Y
    zero("canonical direct density", R.subs(canonical)-Rc*c-Ra)
    zero("canonical direct pressure", P.subs(canonical)-Pc*c-Pa)

    Rnorm2 = L**2+9*L**4/(4*k**2)
    Pnorm2 = 4*k**2/9+5*L**2/3+L**4/(4*k**2)
    zero("density oscillatory Euclidean norm", Rnorm2
         - ((3*L**2/(2*k))**2+L**2))
    zero("pressure oscillatory Euclidean norm", Pnorm2
         - ((2*k/3+L**2/(2*k))**2+L**2))
    Pmajorant = 2*k/3+5*L**2/(4*k)
    zero("pressure positive majorant square gap",
         Pmajorant**2-Pnorm2-21*L**4/(16*k**2))
    q2 = 9*L**2/4
    zero("density positive majorant square gap",
         (k*k+q2/2)**2-k*k*(k*k+q2)-q2*q2/4)

    Rint = L*((K*K+q2)**s.Rational(3,2)-q2**s.Rational(3,2))/(6*Pi**2)
    Rbound = L*K**3/(6*Pi**2)+9*L**3*K/(16*Pi**2)
    Pbound = K**4/(12*Pi**2)+5*L**2*K**2/(16*Pi**2)
    Rcint = (K**4+3*L**2*K**2)/(8*Pi**2)
    Pcint = K**4/(24*Pi**2)+L**2*K**2/(8*Pi**2)
    zero("density exact amplitude integral derivative",
         s.diff(Rint,K) - (measure*s.sqrt(Rnorm2)).subs(k,K))
    zero("density exact amplitude integral lower endpoint", Rint.subs(K,0))
    zero("density rational bound integral derivative",
         s.diff(Rbound,K)
         - (L*(k*k+q2/2)/(2*Pi**2)).subs(k,K))
    zero("pressure rational bound integral derivative",
         s.diff(Pbound,K)-(measure*Pmajorant).subs(k,K))
    zero("density c integral derivative",
         s.diff(Rcint,K)-(measure*Rc).subs(k,K))
    zero("pressure c triangle integral derivative",
         s.diff(Pcint,K)-(measure*(k/3+L**2/(2*k))).subs(k,K))

    La, Lb = s.Rational(2,9), s.Rational(2,7)
    endpoint_c = s.simplify(Rc.subs(L,Lb)-Rc.subs(L,La))
    zero("endpoint work cancellation of leading k c",
         endpoint_c - s.Rational(64,1323)/k)
    Dcbound = 16*K**2/(1323*Pi**2)
    Dabound = (La+Lb)*K**3/(6*Pi**2) + 9*(La**3+Lb**3)*K/(16*Pi**2)
    zero("endpoint work c integral derivative",
         s.diff(Dcbound,K)-(measure*endpoint_c).subs(k,K))
    zero("endpoint A triangle sum",
         Dabound-Rbound.subs(L,La)-Rbound.subs(L,Lb))
    zero("endpoint A rational Pi floor",
         Dabound.subs(Pi,3)-16*K**3/1701-536*K/250047)
    Jcbound = K**4/(24*Pi**2)+(Lb-La)*K**2/(8*Pi**2)
    Jabound = (Rbound.subs(L,La)/La+Rbound.subs(L,Lb)/Lb)/3
    zero("time-integrated-pressure canonical signed coefficient",
         ((-Rc.subs(L,Lb)/Lb+Rc.subs(L,La)/La)/3)
         - (k/3-(Lb-La)/(2*k)))
    zero("time-integrated-pressure canonical triangle integral derivative",
         s.diff(Jcbound,K)
         - (measure*(k/3+(Lb-La)/(2*k))).subs(k,K))
    zero("time-integrated-pressure A rational Pi floor",
         Jabound.subs(Pi,3)-K**3/81-65*K/23814)
    phase_X = ar*s.cos(theta)-ai*s.sin(theta)
    phase_Y = ar*s.sin(theta)+ai*s.cos(theta)
    Rsolution = (Rc*c+Ra).subs({X:phase_X,Y:phase_Y,L:-1/t})
    Psolution = (Pc*c+Pa).subs({X:phase_X,Y:phase_Y,L:-1/t})
    zero("explicit canonical work derivative",
         s.diff(Rsolution,t)-(-1/t)*(Rsolution-3*Psolution))
    zero("explicit canonical pressure primitive",
         s.diff(t*Rsolution,t)-3*Psolution)

    # A bounded raw W error produces A=O(1/k); it remains integrable here.
    omega = s.symbols("omega", nonnegative=True)
    zero("IR bounded-c density weighted coefficient",
         s.limit(measure*Rc,k,0,dir="+"))
    zero("IR bounded-A density weighted amplitude",
         s.limit(measure*s.sqrt(Rnorm2),k,0,dir="+"))
    zero("IR bounded-A pressure weighted amplitude",
         s.limit(measure*s.sqrt(Pnorm2),k,0,dir="+"))
    zero("IR bounded-W density leading weighted coefficient",
         s.limit(measure*3*L**2/(2*k)*omega/(2*k),k,0,dir="+")
         - 3*L**2*omega/(8*Pi**2))
    for exponent in (s.Rational(0), s.Rational(1), s.Rational(3,2)):
        expected = 1/(2-exponent)
        zero(f"IR sufficient power p={exponent}",
             s.integrate(k**(1-exponent),(k,0,1))-expected)
    reject("IR p=2 is not integrable",
           s.limit(s.log(k),k,0,dir="+") - 0)

    witness = {k:1,L:s.Rational(1,4),x:1,v:2,z:3,gr:5,gi:7,c:1,X:2,Y:3}
    reject("wrong pressure geometry sign",
           P.subs(canonical)-((k/3+L**2/(2*k))*c+Pa),witness)
    reject("dropping density geometry contribution",
           R.subs(canonical)-(k*c+L*Y),witness)
    reject("wrong Im W density sign",
           R-((2*k*k+3*L*L)*x+k*z-L*v)/(2*k),witness)
    reject("wrong real W density sign",
           R-((2*k*k+3*L*L)*x-k*z+L*v)/(2*k),witness)
    reject("wrong direct-work pressure coefficient",
           linear_derivative(R)-L*(R-2*P),witness)
    wrong_R = R-3*L**2*c/(2*k)
    wrong_P = P+L**2*c/(2*k)
    zero("paired canonical geometry omission passes Ward",
         linear_derivative(wrong_R)-L*(wrong_R-3*wrong_P))
    reject("paired geometry omission changes direct density",R-wrong_R,witness)
    reject("paired geometry omission changes direct pressure",P-wrong_P,witness)
    reject("missing forced real source term",
           linear_derivative(R,True)-L*(R-3*P)-gi/2,witness)
    reject("missing canonical incoming shift",
           (ur+phase_X).subs(t,a).subs(incoming)-ur,{k:1,wi:2})
    reject("wrong canonical A imaginary sign",
           harmonic_v.subs(t,a).subs(ai,wr/(2*k))-wr,{wr:1})
    reject("wrong oscillator phase sign",
           s.diff(harmonic_x.subs(theta,-theta),t)-harmonic_v,
           {t:a,k:1,ai:1})
    reject("wrong first-order norm sign",
           s.diff(qnorm,eps).subs(eps,0)-(2*x+z/k),witness)
    reject("pressure amplitude underbound",
           (2*k/3+L**2/(2*k))**2-Pnorm2,witness)
    reject("density rational majorant wrong factor",
           (k*k+q2/4)**2-k*k*(k*k+q2),{k:1,L:1})
    reject("canonical work cancellation omitted",
           endpoint_c-(k+3*Lb**2/(2*k)),{k:1})
    reject("pressure primitive wrong factor",
           (R+t*linear_derivative(R)-2*P).subs(L,-1/t),
           {t:-4,k:1,x:1,v:2,z:3})

    rational_coefficients = {
        "R_c": Rcint.subs({L:Lb,Pi:3}),
        "R_A": Rbound.subs({L:Lb,Pi:3}),
        "P_c": Pcint.subs({L:Lb,Pi:3}),
        "P_A": Pbound.subs({L:Lb,Pi:3}),
        "work_c": Dcbound.subs(Pi,3),
        "work_A": Dabound.subs(Pi,3),
        "int_pressure_c": Jcbound.subs(Pi,3),
        "int_pressure_A": Jabound.subs(Pi,3),
    }
    expected_coefficients = {
        "R_c":K**4/72+K**2/294,
        "R_A":K**3/189+K/686,
        "P_c":K**4/216+K**2/882,
        "P_A":K**4/108+5*K**2/1764,
        "work_c":16*K**2/11907,
        "work_A":16*K**3/1701+536*K/250047,
        "int_pressure_c":K**4/216+K**2/1134,
        "int_pressure_A":K**3/81+65*K/23814,
    }
    for name in rational_coefficients:
        zero(f"rational coefficient {name}",
             rational_coefficients[name]-expected_coefficients[name])
    sensitivity = {
        str(cutoff): {name:str(s.factor(expr.subs(K,cutoff)))
                     for name,expr in rational_coefficients.items()}
        for cutoff in (64,128,256)
    }
    scientific = {
        "result":"PASS_CONDITIONAL_LINEAR_INCOMING_STATE_PROPAGATION",
        "hypotheses":[
            "same forcing and contacts; delta g=0",
            "fixed displayed linear direct-stress operators",
            "t in [-9/2,-7/2], L=-1/t, k>0, same Pi>=3",
            "sufficient spectral integrability or stated uniform chi,sigma bounds",
        ],
        "checks":checks,
        "mutation_controls":mutations,
        "rational_sensitivity":sensitivity,
        "full_12_case_pressure_contact_certificate":"UNRESOLVED",
        "historical_metric":"FAIL",
        "higher_dimensional_big_bang_origin":"NOT_ESTABLISHED",
        "incoming_state_error":"NOT_MEASURED",
        "physical_source_evaluations":0,
        "saved_array_decodes":0,
        "physical_trajectory_evaluations":0,
    }
    scientific_bytes = json.dumps(scientific,sort_keys=True,separators=(",",":")).encode()
    receipt = {
        "format":"HDBLAST independent SymPy exact-state receipt v1",
        "python":platform.python_version(),
        "python_optimization":sys.flags.optimize,
        "sympy":s.__version__,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scientific_sha256":hashlib.sha256(scientific_bytes).hexdigest(),
        "check_count":len(checks),
        "mutation_count":len(mutations),
        "scientific":scientific,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({key:receipt[key] for key in (
        "python_optimization","check_count","mutation_count","scientific_sha256"
    )},sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
