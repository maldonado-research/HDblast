"""Exact cap moments for a bulk-volume deficit proportional to (v1²-v2²)²."""
import json
from fractions import Fraction as F


def moments(c):
    # On a cap the cosine t of the angle to its center is uniform on [c,1].
    t2 = (1 + c + c*c) / 3
    t4 = sum(c**i for i in range(5)) / 5
    axis = F(19,8)*t4 - F(7,4)*t2 + F(3,8)
    diagonal = 2*(t2 - t4)
    pole = (1 - 2*t2 + t4) / 2
    return axis, diagonal, pole


def run():
    cutoff = F(4,5)
    axis, diagonal, pole = moments(cutoff)
    sky = F(4,15)
    # Independently integrate the Legendre polynomial P4 over a spherical cap.
    fourth = cutoff*(1+cutoff)*(7*cutoff**2-3)/8
    if (axis-diagonal != fourth or (axis-diagonal)/sky != F(999,1000)
            or not pole < diagonal < axis):
        raise ValueError('angular contrast')
    return {k:str(v) for k,v in {
        'cap_cosine_cutoff':cutoff,'tangent_sky_fraction_per_cap':(1-cutoff)/2,
        'axes':axis, 'equatorial_diagonals':diagonal, 'polar_axis':pole,
        'full_sky':sky, 'axial_minus_diagonal':fourth, 'contrast_over_sky':fourth/sky
    }.items()}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
