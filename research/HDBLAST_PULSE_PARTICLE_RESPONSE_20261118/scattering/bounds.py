"""Exact rational enclosure of noncrossing in/out mixing for every mode in a band."""
from fractions import Fraction as F
import json
from pathlib import Path


def certify(c, indices, h, mass):
    # Frequencies in units of H0. These integer bounds enclose sqrt(mass²+k²).
    lo,hi=F(50),F(54)
    if not lo**2<mass**2+5**2<mass**2+18**2<hi**2:
        raise ValueError('Frequency interval')
    if mass**2+min([F(0),*c])<=0:
        raise ValueError('Squared mass may cross zero')
    center=F(12,100)
    sign=1 if sum(c)>0 else -1
    projection=F(0)
    for j,coef in zip(indices,c,strict=True):
        radius=abs(2*hi*(h*j-center))
        if radius>=3:
            raise ValueError('Cosine bound domain')
        cosine_lower=1-radius**2/2+radius**4/24-radius**6/720
        aligned=sign*coef
        projection+=aligned*(cosine_lower if aligned>=0 else 1)
    sinc_lower=1-(hi*h)**2/6
    if projection<=0 or sinc_lower<=0:
        raise ValueError('No positive Born projection')
    born_lower=h*sinc_lower**4*projection/(2*hi)
    L=h*sum(map(abs,c))/lo
    if L>=3:
        raise ValueError('Geometric bound on exponential is unavailable')
    remainder=L**2/(4*(1-L/3))
    lower=born_lower-remainder
    upper=L/2+remainder
    if not 0<lower<=upper:
        raise ValueError('No certified signal')
    # Round reported lower bounds down and upper bounds up in exact decimal arithmetic.
    def bound(x, upper=False):
        scale=10**12
        i=-((-x*scale)//1) if upper else x*scale//1
        return f'{i//scale}.{i%scale:012d}'
    return dict(born_lower=bound(born_lower),remainder_upper=bound(remainder,True),
                beta_abs_lower=bound(lower),beta_abs_upper=bound(upper,True),
                occupation_lower=bound(lower**2),occupation_upper=bound(upper**2,True))


def main():
    data=json.loads(Path(__file__).with_name('PULSES.json').read_text())
    mass=100*(F('0.5')-F(data['phi0']))
    h=F(data['h'])
    report={n:certify(list(map(F,c)),data['centers_in_h'],h,mass) for n,c in data['protocols'].items()}
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
