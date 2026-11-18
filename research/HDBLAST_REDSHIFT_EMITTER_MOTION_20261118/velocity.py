"""Response of matched-redshift cap contrasts to additional emitter velocities."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json


def compute(path):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != "40b61502f6877d6c91614d2eded6db2c043fd38a69a84928ccc7309bbeb08a09":
        raise ValueError("registered distance coefficients")
    source = json.loads(raw)
    a = F(source['receiver_x'])
    b = F(source['source_background_x'])
    if (a,b)!=(F(49,32),F(99,64)):
        raise ValueError('registered clock')
    ell = F(1,10**6)
    cap = F(18,25)
    mean_forward = F(9,10)
    K = 1+F(1,10**11)
    q0,q1 = a*ell,b*ell
    if not 1/(1-q1**2) < K:
        raise ValueError('cosh envelope')
    # |G|=(sinh(q1)-sinh(q0))*sinh(q1)/cosh(q1)^2.
    leading = (q1-q0)*q1
    jacobian = [leading/K**2,leading*K**2]
    root2 = [F(14142135,10**7),F(14142136,10**7)]
    if not root2[0]**2 < 2 < root2[1]**2:
        raise ValueError('square root')
    # Model physical curvature units, using the wider inherited 129<rho_b<131.
    radius = [129*(q1-q0)/K,131*(q1-q0)*K]
    out={}
    for name, result in source['comparisons'].items():
        low,high = map(F,result['fractional_distance_cap_half_difference'])
        if low*high<=0 or low>high:
            raise ValueError('fiducial sign')
        h=sorted((abs(low),abs(high)))
        bias=[g*y for g,y in zip(jacobian,h)]
        dipole=[s*x/mean_forward for s,x in zip(root2,bias)]
        axis=[x/cap for x in bias]
        gradient=[bias[0]/radius[1],bias[1]/radius[0]]
        out[name] = {'minimum_uniform_radial_speed_for_two_caps':bias,
                     'radial_quadrupole_coefficient_for_full_sky_cancellation':axis,
                     'minimum_common_dipole_speed_for_two_one_sided_caps':dipole,
                     'antipodal_shear_necessary_threshold':gradient}
    return {'registered_distance_sha256':hashlib.sha256(raw).hexdigest(),
            'absolute_background_log_redshift_slope_times_length':jacobian,
            'cap_area_factor':cap,'cap_dipole_factor':mean_forward,
            'motion_thresholds_per_unit_canonical_amplitude':out}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('redshift',nargs='?',type=Path,default=Path(__file__).resolve().parent/'REDSHIFT.json')
    print(json.dumps(compute(parser.parse_args().redshift),indent=2,default=str))
