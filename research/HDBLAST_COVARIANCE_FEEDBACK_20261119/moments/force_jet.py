"""Check the initial force moments of a closed Gaussian scalar on a shared smooth shell."""
import json
from fractions import Fraction as F
from math import comb


def mode_jets(order, potential, hubble):
    # Polynomial coefficients in lambda=q^2, with Delta C factored out.
    zero = [F(0)]*(order+1)
    rows = [zero.copy(), zero.copy(), zero.copy()]
    rows[1][0], rows[2][0] = 2, 4*hubble
    for n in range(order-2):
        row = zero.copy()
        for j in range(order):
            row[j+1] -= 4 * rows[n+1][j]
        for k in range(n+1):
            weight = comb(n,k)
            for j in range(order+1):
                row[j] -= weight*(4*potential[k]*rows[n-k+1][j]
                                  + 2*potential[k+1]*rows[n-k][j])
        rows.append(row)
    return rows


def discrete_recovery(rows):
    """Recover one covariance difference per positive-weight momentum node."""
    lambdas = [F(1,9), F(4,9), F(1), F(16,9), F(25,9)]
    weights = [F(1,8), F(1,8), F(1,4), F(1,4), F(1,4)]
    initial = [F(-1,5), F(2,9), F(1,8), F(-2,7), F(1,3)]
    moments = []
    for r in range(len(lambdas)):
        coefficients = rows[2*r+1][:r+1]
        observed = sum(w*c*sum(a*x**j for j,a in enumerate(coefficients))
                       for x,w,c in zip(lambdas,weights,initial))
        known = sum(a*m for a,m in zip(coefficients,moments))
        moments.append((observed-known)/coefficients[r])
    recovered = []
    for target,weight in zip(lambdas,weights):
        poly = [F(1)]
        for x in lambdas:
            if x == target:
                continue
            scaled = [a/(target-x) for a in poly]
            poly = [a-x*b for a,b in zip([0]+scaled,scaled+[0])]
        recovered.append(sum(a*m for a,m in zip(poly,moments))/weight)
    if recovered != initial:
        raise ValueError('Covariance interpolation')
    return [str(v) for v in recovered]


def check():
    # Nonconstant exact rational background derivatives U^(k) at the fork.
    order = 23
    potential = [F((-1)**k*(k+2),k+3) for k in range(order+1)]
    rows = mode_jets(order, potential, F(2,7))
    for r in range((order+1)//2):
        odd = rows[2*r+1]
        if any(odd[r+1:]) or odd[r] != 2*(-4)**r:
            raise ValueError('Odd-derivative polynomial identity')
    # S(x) in units q/Q: normalized initial correlation difference, up to
    # an arbitrary common scale. Integrate unit-band weight x^2 exactly.
    moments = [1/F(2*j+3) - F(14,3*(2*j+5)) + F(21,5*(2*j+7))
               for j in range(order+1)]
    force_jets = [sum(a*b for a,b in zip(row,moments)) for row in rows]
    if any(force_jets[:5]) or force_jets[5]!=32*moments[2] or not force_jets[5]:
        raise ValueError('Moment cancellation example')
    return dict(verified_odd_derivatives=order//2+1,highest_derivative=order,
                recovered_five_node_covariance=discrete_recovery(rows),
                normalized_moments=[str(v) for v in moments[:3]],
                reduced_force_derivatives_through_fifth=[str(v) for v in force_jets[:6]])

if __name__=='__main__':
    print(json.dumps(check(),indent=2))
