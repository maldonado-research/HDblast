#!/usr/bin/env python3
"""Pure symbolic specialization of the generic metric WKB inventory.

Produces elementwise common-subexpression code, without evaluating a physical
frequency, source, mode, or integral. The exact derivative identities are
retained in SPECIALIZATION_PROOF.json.
"""
from pathlib import Path
import hashlib
import json
import sympy as s
from sympy.printing.pycode import PythonCodePrinter


class ExtendedPrinter(PythonCodePrinter):
    def _print_Rational(self, expr):
        return f"(LD({expr.p})/LD({expr.q}))"


def main():
    here = Path(__file__).resolve().parent
    w = s.symbols("w0:5", nonzero=True)
    z = s.symbols("z0:5")
    C = s.symbols("C0:3")
    dC = s.symbols("dC0:3")
    L,dL,M,dM,k2 = s.symbols("L dL M dM k2")
    inputs = (*w,*C,L,M)
    variations = (*z,*dC,dL,dM)
    derivative_pairs = list(zip(w[:-1],w[1:]))+list(zip(C[:-1],C[1:]))
    def D(expr):
        return sum(s.diff(expr,a)*b for a,b in derivative_pairs)
    def delta(expr):
        return sum(s.diff(expr,a)*b for a,b in zip(inputs,variations))
    U = -C[0]/(2*w[0])-w[2]/(4*w[0]**2)+3*w[1]**2/(8*w[0]**3)
    U1 = -C[1]/(2*w[0])+C[0]*w[1]/(2*w[0]**2)-w[3]/(4*w[0]**2)+5*w[1]*w[2]/(4*w[0]**3)-9*w[1]**3/(8*w[0]**4)
    U2 = (-C[2]/(2*w[0])+C[1]*w[1]/w[0]**2+C[0]*w[2]/(2*w[0]**2)
          -C[0]*w[1]**2/w[0]**3-w[4]/(4*w[0]**2)
          +7*w[1]*w[3]/(4*w[0]**3)+5*w[2]**2/(4*w[0]**3)
          -57*w[1]**2*w[2]/(8*w[0]**4)+9*w[1]**4/(2*w[0]**5))
    checks = {}
    for name,residual in (("W2_prime",D(U)-U1),("W2_second",D(U1)-U2)):
        require_zero = s.expand(residual)
        if require_zero != 0:
            raise RuntimeError(name+" symbolic identity failed: "+str(s.factor(require_zero)))
        checks[name] = "exact_zero"
    V = -U**2/(2*w[0])-U2/(4*w[0]**2)+w[2]*U/(4*w[0]**3)+3*w[1]*U1/(4*w[0]**3)-3*w[1]**2*U/(4*w[0]**4)
    b = -k2/3-M
    c = 1-b/w[0]**2
    J2 = L*w[1]/w[0]**2+w[1]**2/(4*w[0]**3)
    J4 = L*U1/w[0]**2-2*L*w[1]*U/w[0]**3+w[1]*U1/(2*w[0]**3)-3*w[1]**2*U/(4*w[0]**4)
    R0,R2,R4 = w[0]/2,(L**2/w[0]+J2)/4,(U**2/w[0]-L**2*U/w[0]**2+J4)/4
    P0,P2,P4 = k2/(6*w[0]),(c*U+L**2/w[0]+J2)/4,(c*V+b*U**2/w[0]**3-L**2*U/w[0]**2+J4)/4
    S0,S2 = 1/(2*w[0]),-U/(2*w[0]**2)
    S=S0+S2
    S1=-w[1]/(2*w[0]**2)-U1/(2*w[0]**2)+U*w[1]/w[0]**3
    Ssecond=(-w[2]/(2*w[0]**2)+w[1]**2/w[0]**3-U2/(2*w[0]**2)
             +2*U1*w[1]/w[0]**3+U*w[2]/w[0]**3-3*U*w[1]**2/w[0]**4)
    for name,residual in (("S_prime",D(S)-S1),("S_second",D(S1)-Ssecond)):
        if s.expand(residual) != 0:
            raise RuntimeError(name+" symbolic identity failed")
        checks[name] = "exact_zero"
    inventory={"R0":R0,"R2":R2,"R4":R4,"P0":P0,"P2":P2,"P4":P4,
               "S0":S0,"S2":S2,"R":R0+R2+R4,"P":P0+P2+P4,"S":S,
               "W2":U,"W2_prime":U1,"W2_second":U2,"W4":V,"J2":J2,"J4":J4,
               "S_prime":S1,"S_second":Ssecond}
    outputs=dict(inventory)
    outputs.update({"delta"+key:delta(expr) for key,expr in inventory.items()})
    replacements,reduced=s.cse(list(outputs.values()),symbols=s.numbered_symbols("t"),order="canonical")
    printer=ExtendedPrinter()
    lines=["\"\"\"Generated pure directional WKB specialization; regenerate with specialize_wkb.py.\"\"\"",
           "import numpy as np", "LD = np.longdouble", "", "def evaluate_wkb(k2, L, M, omega, delta_omega, C, delta_C, delta_L, delta_M):",
           "    w0,w1,w2,w3,w4 = omega", "    z0,z1,z2,z3,z4 = delta_omega", "    C0,C1,C2 = C", "    dC0,dC1,dC2 = delta_C", "    dL,dM = delta_L,delta_M"]
    lines.extend("    "+str(symbol)+" = "+printer.doprint(expr) for symbol,expr in replacements)
    lines.append("    return {")
    lines.extend("        "+repr(key)+": "+printer.doprint(expr)+"," for key,expr in zip(outputs,reduced))
    lines.extend(["    }",""])
    code="\n".join(lines)
    path=here/"metric_wkb.py"
    path.write_text(code)
    proof={"schema_version":1,"kind":"pure_exact_symbolic_specialization","physical_evaluations":0,
           "exact_checks":checks,"directional_rule":"sum_i delta_input_i * partial_input_i(full_inventory)",
           "common_subexpressions":len(replacements),"output_inventory":list(outputs),
           "generated_sha256":hashlib.sha256(code.encode()).hexdigest(),
           "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (here/"SPECIALIZATION_PROOF.json").write_text(json.dumps(proof,indent=2)+"\n")
    print(json.dumps(proof,indent=2))


if __name__ == "__main__":
    main()
