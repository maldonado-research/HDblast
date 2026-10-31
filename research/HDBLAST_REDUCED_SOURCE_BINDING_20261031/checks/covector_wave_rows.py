"""Eight exact physical covector (Box0+Ric0) modal rows. Standard library only."""
from fractions import Fraction as F


def wave_rows(z,deriv,rho,a,p,c,theta):
    """Return (Box0+Ric0)Z amplitudes; deriv(v,0)=Dy, deriv(v,1)=Dt.

    z has y,0,a rows of lengths 3,3,2 in the degree 0,2,4 scalar and
    degree 2,4 gradient bases. Ordinary derivatives apply to full coefficients.
    This function is independently checked against full-frame contraction.
    """
    potential=6/rho**2+p*p/2-6*a*a
    ric=F(2,3)*potential
    def scalar(v,eigenvalue):
        return (deriv(deriv(v,0),0)+4*a*deriv(v,0)
                -(deriv(deriv(v,1),1)+3*theta*deriv(v,1)+eigenvalue*v/c**2)/rho**2)
    out={'y':[],'0':[],'a':[]}
    for j,l in enumerate((0,2,4)):
        lam=l*(l+2)
        radial=z['y'][j];temporal=z['0'][j]
        angular=0 if j==0 else z['a'][j-1]
        out['y'].append(scalar(radial,lam)+(ric+p*p-4*a*a)*radial
                        +2*a/rho*(deriv(temporal,1)+3*theta*temporal)
                        +2*a*lam*angular/(rho*c))
        out['0'].append(scalar(temporal,lam)+(ric-a*a+3*theta*theta/rho**2)*temporal
                        +2*a*deriv(radial,1)/rho+2*lam*theta*angular/(rho**2*c))
        if j:
            out['a'].append(scalar(angular,lam-2)+(ric-a*a+theta*theta/rho**2)*angular
                            +2*a*radial/(rho*c)-2*theta*temporal/(rho**2*c))
    return out
