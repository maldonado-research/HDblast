"""Verify the numerical choices in the spectral cutoff argument."""
import json
from fractions import Fraction as F

ell = F(1, 10**6)
cutoff = 10**16
h = F(1, 10**8)
hat_error = F(6, 100000)
energy_error = F(101, 1000000)
form_bound = 1/h + 103000*ell**2*h/3 + 300*ell
if h/3 >= hat_error**2 or form_bound >= cutoff * energy_error**2:
    raise ValueError('endpoint projection error')
output_error = F(3, 250) * 417 / ell * (hat_error + energy_error)
if output_error != F(201411, 250) or output_error >= 806:
    raise ValueError('physical readout error bound')
# Since pi > 3, 6*sqrt(cutoff)/pi is strictly below the integer below.
max_modes = 2 * 10**8
if max_modes**2 != 4*cutoff:
    raise ValueError('Dirichlet eigenvalue count')
print(json.dumps({'scaled_spectral_cutoff': cutoff, 'maximum_number_of_fields': max_modes,
                  'strict_uniform_error_ceiling': str(output_error)}, indent=2))
