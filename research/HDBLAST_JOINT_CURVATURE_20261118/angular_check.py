"""Evaluate the full-cut null-test coefficients by exact sphere moments."""
from fractions import Fraction as F
import json
import math


def sphere_moment(exponents):
    """Normalized monomial integral on S^3 from Gaussian radial factorization."""
    if any(power % 2 for power in exponents):
        return F(0)
    numerator = math.prod(factor for power in exponents for factor in range(1, power, 2))
    denominator = math.prod(range(4, 4 + sum(exponents), 2))
    return F(numerator, denominator)


def integral(a, b):
    return sum((ca*cb*sphere_moment(tuple(i+j for i,j in zip(pa,pb)))
                for ca,pa in a for cb,pb in b), F(0))


one = [(1,(0,0,0,0))]
y = [(1,(2,0,0,0)), (-1,(0,2,0,0))]
z = [(1,(0,0,2,0)), (-1,(0,0,0,2))]
e_y = [(2,(1,0,1,0)), (-2,(0,1,0,1))]
hess = [(2*c,p) for c,p in z] + [(-2*c,p) for c,p in y]
rows = []
for sign in (1,-1):
    weight = one + [(sign*c,p) for c,p in z]
    vals = [integral(weight,target) for target in (one,y,e_y,hess)]
    if vals != [F(1),F(0),F(0),F(sign,3)]:
        raise ValueError('angular moment contractions')
    rows.append({'shutter': '1+Z' if sign==1 else '1-Z',
                 'total_mass':str(vals[0]),'weighted_Y':str(vals[1]),
                 'weighted_eY':str(vals[2]),'weighted_Hessian':str(vals[3]),
                 'null_Ricci_per_U_over_a_squared':str(-2*vals[3])})
print(json.dumps({'dimension_of_sphere':3,'moments':rows},indent=2))
