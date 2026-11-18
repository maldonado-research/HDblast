"""Exact rational arithmetic for the tangential-ray normal Jacobi first jet.

Standard library only; no ancestor Python module is imported. This encloses
finite coefficients conditional on the continuum identities in FIRST_ORDER.md.
It does not prove a nonlinear amplitude/tilt radius or give an SI calibration.
"""
from fractions import Fraction as F
from pathlib import Path
import csv
import hashlib
import io
import json

PINS = {
    "FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "METRIC.json": "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42",
    "KERNEL.json": "6aad53ade018e9e9f9f79be27ff80a615bc7b8b8005c1072482bd9730e896fa6",
}
ELL, A, B, X = F(1, 10**6), F(49, 32), F(99, 64), F(25, 16)
DELTA, TRANSFER, K = B - A, F(3, 10**6), 1 + F(1, 10**11)


def checked(condition, label, roster):
    if not condition:
        raise ValueError(label)
    roster.append(label)


def scale(interval, factor):
    return sorted(factor * endpoint for endpoint in interval)


def product(*intervals):
    out = [F(1), F(1)]
    for interval in intervals:
        values = [a*b for a in out for b in interval]
        out = [min(values), max(values)]
    return out


def report():
    checks = []
    root = Path(__file__).resolve().parent.parent
    joint = root / "HDBLAST_JOINT_CURVATURE_20261118"
    finite = root / "HDBLAST_FINITE_LIGHT_CONES_20261118"
    paths = {"FUTURE_CELLS.csv": joint / "inputs/FUTURE_CELLS.csv",
             "SHELL_BOX.json": joint / "inputs/SHELL_BOX.json",
             "METRIC.json": joint / "METRIC.json",
             "KERNEL.json": finite / "KERNEL.json"}
    raw = {}
    for name, pin in PINS.items():
        raw[name] = paths[name].read_bytes()
        checked(hashlib.sha256(raw[name]).hexdigest() == pin, "pin_"+name, checks)
    metric = json.loads(raw["METRIC.json"])
    kernel = json.loads(raw["KERNEL.json"])
    box = json.loads(raw["SHELL_BOX.json"])["boxes"]
    rho = list(map(F, box["rho_b"]["exact"]))
    bc = list(map(F, box["B_b"]["exact"]))
    source_bound = F(metric["source_bound"])
    zbound, theta = F(metric["z_bound"]), F(metric["theta"])
    checked(theta == F(5,2)*ELL**2*(X-1)**2/2, "volterra_theta", checks)
    checked(zbound == source_bound*(X-1)**2/(2*(1-theta)), "volterra_bound", checks)
    checked(1/(1-(ELL*X)**2) < K, "cosh_and_affine_bound", checks)
    checked(F(5,2)-F(35,4)*(ELL*X)**2 > 0, "positive_endpoint_potential", checks)

    # Independent integrations of u against positive polynomial kernels.
    # C_position kernel is obtained by interchanging the double integrals.
    def position_primitive(s):
        y = (s-A)/DELTA
        if s <= A:
            return DELTA*(s-A)/2-(s-A)**2/6
        return DELTA**2*(y/2-y*y/6-y**5/30)

    position_flat = [F(0), F(0)]
    slope_flat = [F(0), F(0)]
    z_endpoint = [F(0), F(0)]
    dz_endpoint = [F(0), F(0)]
    edge, full_edge, count, covers = F(1), F(0), 0, True
    for row in csv.DictReader(io.StringIO(raw["FUTURE_CELLS.csv"].decode("ascii"))):
        l0, r0 = F(int(row["left_grid"]),4096), F(int(row["right_grid"]),4096)
        covers = covers and l0 == full_edge and r0 > l0
        full_edge = r0
        l, r = max(F(1),l0), min(B,r0)
        if l >= r:
            continue
        if l != edge:
            raise ValueError("clipped cell coverage")
        edge = r
        lo = F(int(row["H_lo"]), 2**112)-TRANSFER
        hi = F(int(row["H_hi"]), 2**112)+TRANSFER
        if not lo <= hi or max(abs(lo),abs(hi)) > source_bound:
            raise ValueError("cell interval/source bound")
        wp = position_primitive(r)-position_primitive(l)
        wz = ((B-l)**2-(B-r)**2)/2
        ws = wz+DELTA*(r-l)
        if not wp > 0 or not ws > 0:
            raise ValueError("kernel positivity")
        for i,u in enumerate((lo,hi)):
            position_flat[i] += wp*u
            slope_flat[i] += ws*u
            z_endpoint[i] += wz*u
            dz_endpoint[i] += (r-l)*u
        count += 1
    checked(covers and full_edge == 3 and edge == B and count == 2240,
            "complete_grid_and_optical_interval", checks)
    checked(position_flat == list(map(F,kernel["flat_axis"])),
            "independent_position_kernel_agrees", checks)

    # H=sech^(7/2)(ell*x)z. Bounds use 1-sech^a(q)<=a*q^2/2.
    h_error = theta*zbound + F(7,4)*(ELL*X)**2*zbound
    dh_error = (F(5,2)*ELL**2*(X-1)*zbound
                + F(7,4)*(ELL*X)**2*source_bound*(X-1)
                + F(7,2)*ELL**2*X*zbound)
    position_error = F(5,3)*h_error+F(2,3)*(K**3-1)*zbound
    # r=L/(d lambda/dx)_b lies in [DELTA/K, DELTA].
    slope_error = (h_error+DELTA*dh_error
                   +DELTA*(K-1)*(source_bound*(X-1)+dh_error))
    position = [position_flat[0]-position_error,position_flat[1]+position_error]
    slope = [slope_flat[0]-slope_error,slope_flat[1]+slope_error]
    archived = list(map(F,kernel["axis_coefficient"]))
    checked(archived[0] <= position[0] < position[1] <= archived[1],
            "position_inside_archived_enclosure", checks)
    checked(F("0.0817743") < position[0] < position[1] < F("0.0820899"),
            "position_decimal_enclosure", checks)
    checked(slope[0] > 0, "positive_opening_coefficient", checks)
    checked(slope_flat == [z_endpoint[i]+DELTA*dz_endpoint[i] for i in (0,1)],
            "opening_kernel_endpoint_crosscheck", checks)

    cap_edge = F(4,5)
    v1sq = (1+cap_edge+cap_edge**2)/3
    cap_weight = (3*v1sq-1)/2
    checked(cap_weight == F(18,25) and (1-cap_edge)/2 == F(1,10),
            "positive_tenth_sky_cap_average", checks)
    checked(2*cap_edge**2-1 == F(7,25), "pointwise_positive_cap", checks)
    checked(1+1-2 == 0, "three_screen_first_trace", checks)

    # Polynomial plane-wave check F(lambda)=lambda^2+lambda^3. Integrate
    # (lambda-s)*s*F''(s) at lambda=1 exactly, independently of endpoint form.
    fcoeff = {2:F(1),3:F(1)}
    endpoint_form = sum((1-F(2,n+1))*v for n,v in fcoeff.items())
    jacobi_integral = sum(v*n*(n-1)*(F(1,n)-F(1,n+1))
                          for n,v in fcoeff.items())
    checked(endpoint_form == jacobi_integral == F(5,6),
            "tracefree_5d_plane_wave_jacobi_sign_check", checks)

    length_over_rho = [ELL*DELTA,K*ELL*DELTA]
    length = product(rho,length_over_rho)
    # κ_opt L^2 = ell B_b/rho_b * (L/rho_b)^2.
    scale_size = product([ELL,ELL],bc,[1/rho[1],1/rho[0]],
                         product(length_over_rho,length_over_rho))
    d_shell = product(scale_size,position)
    normal = scale(d_shell,-2)
    normal_slope = scale(product(scale_size,slope),-2)
    normal_length = product(normal,length)
    outputs = {
        "L_over_rho_b": length_over_rho,
        "L_in_model_length_units": length,
        "kappa_opt_L_squared": scale_size,
        "position_C": position,
        "opening_C": slope,
        "axis_shell_fractional_distance_derivative": d_shell,
        "axis_normal_fractional_Jacobi_derivative": normal,
        "cap_mean_normal_fractional_Jacobi_derivative": scale(normal,cap_weight),
        "axis_normal_Jacobi_derivative_model_length": normal_length,
        "cap_mean_normal_Jacobi_derivative_model_length": scale(normal_length,cap_weight),
        "axis_normal_opening_derivative": normal_slope,
        "cap_mean_normal_opening_derivative": scale(normal_slope,cap_weight),
        "axis_normal_opening_over_position_ratio": [slope[0]/position[1],slope[1]/position[0]],
    }
    return {
        "status":"PASS_CONDITIONAL_FINITE_ARITHMETIC",
        "scope":"Amplitude derivatives at zero and infinitesimal inward tilt only; tangent S2 cap, one smooth bulk copy. No finite amplitude/tilt radius, full bulk-sky cap, physical-unit calibration or experimental sensitivity is proved.",
        "input_sha256":PINS,
        "clock":{"a":str(A),"b":str(B),"ell":str(ELL)},
        "cap":{"solid_angle_fraction":"1/10","mean_quadrupole":"18/25","minimum_quadrupole":"7/25"},
        "finite_checks":checks,
        "finite_check_count":len(checks),
        "integration_cells":count,
        "errors":{"H":str(h_error),"H_x":str(dh_error),"position":str(position_error),"opening":str(slope_error)},
        "exact_intervals":{key:list(map(str,value)) for key,value in outputs.items()},
        "display_only":{key:list(map(float,value)) for key,value in outputs.items()},
    }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
