"""Exact rational witness for the endpoint controllability determinant at the shell."""
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def witness():
    path = Path(__file__).resolve().parents[1] / 'B1_summary.json'
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != 'bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f':
        raise ValueError('Invalid source fixture')
    model = json.loads(data, parse_float=Decimal)
    field = Fraction(model['phi_b_static'])
    delta, c, d = (Fraction(model['params'][i]) for i in ('delta','c','d'))
    w = 1-field+field**3/3
    wp = field**2-1
    potential = wp**2/2-2*w**2/3
    potential_prime = 2*field*wp-4*w*wp/3
    tension = 2*w+delta*(1+c*field+d*field**2/2)
    tension_prime = 2*wp+delta*(c+d*field)
    # Divide the controllability determinant by the positive cube of the wall radius.
    determinant = tension*tension_prime**2/24+tension_prime*potential_prime/3
    inverse_radius_squared = tension**2/36-tension_prime**2/48+potential/6
    if determinant >= 0 or inverse_radius_squared <= 0:
        raise ValueError('Nonzero-rank or regular-background condition failed')
    return {
        'determinant_over_radius_cubed': str(determinant),
        'inverse_radius_squared': str(inverse_radius_squared),
        'determinant_over_radius_cubed_approx': float(determinant),
    }


if __name__ == '__main__':
    print(json.dumps(witness(),indent=2))
