"""Uniform real-frequency Gram bounds for early and late clock exposures."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

PINS = {
    "SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "SHAPE.json": "247d0a9cd51aada5dac103732278429981aeade6d1f6943cb3ee179d629eec1a",
}


def report(directory):
    def load(name):
        blob = (directory / name).read_bytes()
        if hashlib.sha256(blob).hexdigest() != PINS[name]:
            raise ValueError(name)
        return json.loads(blob)
    source = load("SHELL_BOX.json")
    late_lower = F(load("SHAPE.json")["scaled_window_means"][0][0])
    hubble_lower = F(source["boxes"]["H_squared"]["exact"][0])
    ell, horizon, h = F(1,10**6), F(11,8), F(19,40)
    center = h/2
    midpoint = F(21,16)
    p_cut = F(1,2)
    cos_floor = 1-p_cut**2/12
    scaled_sin_floor = F(1,12)-p_cut**2/240
    sinc_upper = 1-p_cut**2/6+p_cut**4/120
    high_gram_floor = (1-sinc_upper)/2
    def low_ratio_squared(time):
        return 1/cos_floor + ((time-center)/h)**2/scaled_sin_floor
    mean2 = max(low_ratio_squared(midpoint),1/high_gram_floor)
    sup2 = max(low_ratio_squared(horizon),1/high_gram_floor)
    k = F(35,8)*ell**4*horizon**4
    alpha = 9*(k/(1-k)+3*ell**2*horizon**2/4)
    frozen_ratio = F('7.954')
    if not (mean2 < frozen_ratio**2 and sup2 < 9**2 and
            alpha < F(1,10**10) and (frozen_ratio+alpha)/(1-alpha) < 8):
        raise ValueError("single-mode transfer constants")
    data = {}
    for abs_error in (F(0),F(1,2000)):
        field_power = (late_lower-abs_error)**2/64
        trace = F(2,3)*hubble_lower*(1-(ell*h)**2)*field_power
        if field_power <= F('0.2515') or trace <= F('9.9363e-6'):
            raise ValueError("reported gravitational floor")
        data[str(abs_error)] = {"D_squared_times_early_mean_sum_amplitude_squared_lower": field_power,
                                "D_squared_times_QE_over_kappa_lower": trace}
    return {"input_sha256": PINS, "early_window": [F(0),h], "late_window": [F(5,4),horizon],
            "low_cosine_gram_lower":cos_floor, "low_scaled_sine_gram_lower":scaled_sin_floor,
            "high_gram_lower":high_gram_floor, "flat_late_ratio_squared_upper":mean2,
            "flat_sup_ratio_squared_upper":sup2, "distortion_per_early_flat_rms":alpha,
            "physical_transfer_constant":8, "model_lower_bounds_by_late_scalar_error":data}


if __name__ == "__main__":
    print(json.dumps(report(Path(__file__).resolve().parent/"inputs"),indent=2,default=str))
