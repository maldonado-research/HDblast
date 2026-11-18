"""Exact two-field check of the conformal trace expansion, including cross gradients."""
import json
from fractions import Fraction as F


def dot(left, right):
    return sum((x*y for x, y in zip(left, right, strict=True)), F(0))


def check():
    # Algebraic test values. These are not numerical fits to the shell.
    R0, kappa, F1, F2 = map(F, ("3/7", "2/3", "-5/9", "1/11"))
    masses_squared = [F(2), F(5)]
    d = [F(3, 4), F(-1, 3)]
    b = [-F1*d_i / (2*R0 + 6*m) for m,d_i in zip(masses_squared,d)]
    H = [[((kappa if i==j else 0)-6*b[i]*b[j]-d[i]*d[j])/6
          for j in range(2)] for i in range(2)]
    val = [F(2, 5), F(-3, 7)]
    blocks = [[[F(5), F(2)], [F(2), F(-3)]],
              [[F(-2), F(1, 7)], [F(1, 7), F(1, 3)]]]
    mass_val=[m*v for m,v in zip(masses_squared,val)]
    b0=dot(b,val)
    h_val=[dot(row,val) for row in H]
    b2=dot(val,h_val)/2
    d0=dot(d,val)
    # Potential part P(psi) of R_J through orders one and two.
    p1=-2*R0*b0-6*dot(b,mass_val)
    p2=R0*(2*b0*b0-2*b2)+12*b0*dot(b,mass_val)+2*kappa*dot(val,mass_val)-6*dot(h_val,mass_val)
    s2=(p2-F2*d0**2/2)/F1
    if p1!=F1*d0 or p2!=F1*s2+F2*d0*d0/2:
        raise ValueError("field-value terms")
    answers=[]
    for gram in blocks:
        ricci_gradient=kappa*sum(gram[i][i] for i in range(2))
        ricci_gradient-=6*sum((H[i][j]+b[i]*b[j])*gram[i][j] for i in range(2) for j in range(2))
        scalar_gradient=sum(d[i]*d[j]*gram[i][j] for i in range(2) for j in range(2))
        answers.append(ricci_gradient-scalar_gradient)
    if any(answers):raise ValueError("gradient terms")
    return {"linear_couplings":list(map(str,b)),"conformal_hessian":[list(map(str,row)) for row in H],
            "linear_potential":str(p1),"quadratic_potential":str(p2),"quadratic_scalar":str(s2),
            "gradient_residuals":list(map(str,answers))}

if __name__=="__main__":
    print(json.dumps(check(),indent=2))
