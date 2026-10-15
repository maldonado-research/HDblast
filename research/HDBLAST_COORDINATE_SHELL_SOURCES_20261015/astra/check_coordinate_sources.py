#!/usr/bin/env python3
"""Independent exact boundary-source audit; standard library only.

No production/planning/upstream scientific Python module is read or imported.
The retained analytic waveform certificate is an input, not re-proved here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from fractions import Fraction as F

PINS = {
    "inputs/prior/INERTIAL_WAVEFORM_CELLS.csv": (1043169, "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7"),
    "inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json": (10082, "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8"),
    "inputs/prior/INERTIAL_SHELL_WAVEFORM.md": (14813, "7897073c52ee87514b26d9377ac9b3cd88a4b082c47f35bf36164e6f97ad2458"),
    "inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json": (19535, "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8"),
    "inputs/reference/PHYSICAL_SHELL_RESPONSE.md": (16284, "a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b"),
    "inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md": (15786, "8b4d08bf3c4b1704be99e291049a61c902feb39439f8edbda54a5d7d3f9735c8"),
    "inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json": (31129, "fd06ac6eee239834f052c5bbcf7e29531cc127d1fa6621b3069ab4dcf1a77831"),
    "inputs/reference/FULL_CONTACT_MANIFEST.json": (14212, "0ee47e71ac4371a0f118e2c9a6bff835181a8300e2805423a444ffe700eac843"),
    "inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json": (13270, "9b42295d819f0ed5caf94d0a46af38274a27fb11648dd6bff91f5d992385f140"),
    "inputs/reference/COMPLETE_L4_SCALAR_BULK.md": (17865, "e7e23e4576fe89ba072241b981a250ec781f90c0b315a2616ce420f80c589e36"),
}
CONTROLS = {
    "normal_quadratic": "unit_normal_series",
    "one_copy_doubled": "one_copy_scalar_source",
    "taylor_half": "derivative_taylor_coefficient",
    "sk_metric_sign": "metric_coordinate_source",
    "angular_norm": "angular_total_norm",
    "volume_power": "proper_volume_power",
    "scaled_force": "scaled_boundary_force",
    "drop_resolvent": "volterra_resolvent",
    "central_mu": "positive_central_inverse_mu",
    "tail_drop_square": "strict_tail_scalar_sign",
    "proper_derivative": "proper_clock_chain_rule",
    "use_comparator": "actual_trace_envelope",
    "uncancelled_gamma": "dynamic_boundary_gamma_cancellation",
    "metric_radius": "metric_reconstruction_radius",
    "endpoint_drop_square": "central_endpoint_source_signs:central_left",
    "erase_central_metric": "both_central_source_signs",
    "ratio_denominator": "norm_ratio_interval_upper",
}


class Failure(Exception):
    def __init__(self, gate, detail=""):
        self.gate, self.detail = gate, detail
        super().__init__(gate + (": " + detail if detail else ""))


class Poly:
    """Small exact polynomial ring, eight variables, no CAS dependency."""
    def __init__(self, terms=0):
        if isinstance(terms, Poly):
            terms = terms.t
        if not isinstance(terms, dict):
            terms = {(0,) * 8: F(terms)}
        self.t = {tuple(k): F(v) for k, v in terms.items() if v}

    @classmethod
    def var(cls, i):
        k = [0] * 8
        k[i] = 1
        return cls({tuple(k): 1})

    def __add__(self, other):
        out = dict(self.t)
        for k, v in Poly(other).t.items():
            out[k] = out.get(k, 0) + v
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) - self

    def __mul__(self, other):
        out = {}
        for k, v in self.t.items():
            for l, w in Poly(other).t.items():
                m = tuple(a + b for a, b in zip(k, l))
                out[m] = out.get(m, 0) + v * w
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, n):
        out = Poly(1)
        for _ in range(n):
            out = out * self
        return out

    def __eq__(self, other):
        return self.t == Poly(other).t

    def diff(self, i):
        out = {}
        for k, v in self.t.items():
            if k[i]:
                l = list(k)
                l[i] -= 1
                out[tuple(l)] = k[i] * v
        return Poly(out)

    def subs(self, mapping):
        out = Poly()
        for k, v in self.t.items():
            term = Poly(v)
            for i, degree in enumerate(k):
                term *= mapping.get(i, Poly.var(i)) ** degree
            out += term
        return out


def jet_mul(x, y):
    return [sum((x[i] * y[n-i] for i in range(n+1)), Poly()) for n in range(3)]


def sphere_mean(poly):
    # S^3 coordinate moments: E prod n_i^(2k_i) = prod (1/2)_ki / (2)_sumki.
    total = F(0)
    for powers, coeff in poly.t.items():
        if any(powers[4:]):
            raise Failure("moment_dimension")
        if any(p % 2 for p in powers[:4]):
            continue
        ks = [p // 2 for p in powers[:4]]
        term = coeff
        for k in ks:
            for j in range(k):
                term *= F(2*j+1, 2)
        for j in range(sum(ks)):
            term /= 2+j
        total += term
    return total


def box(v):
    return tuple(F(x) for x in v)


def packed(v):
    if isinstance(v, F):
        return {"exact": str(v), "display": float(v)}
    if isinstance(v, dict):
        return {k: packed(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [packed(x) for x in v]
    return v


def audit(root, inject):
    groups = []

    def gate(name, condition, detail=""):
        if not condition:
            raise Failure(name, detail)
        groups.append(name)

    authenticated = []
    registered = json.loads((root / "inputs/INPUT_PINS.json").read_text())
    registry = {x["path"]: x for x in registered["files"]}
    for rel, (size, sha) in PINS.items():
        raw = (root / rel).read_bytes()
        gate("pin:"+rel, len(raw) == size and hashlib.sha256(raw).hexdigest() == sha)
        entry = registry[rel]
        gate("registry:"+rel, entry["bytes"] == size and entry["sha256"] == sha)
        authenticated.append({"path": rel, "sha256": sha, "bytes": size,
                              "public_source": entry["original_public_source"]})
    for manifest_rel, member, local in [
        ("inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json", "PHYSICAL_SHELL_RESPONSE.md", "inputs/reference/PHYSICAL_SHELL_RESPONSE.md"),
        ("inputs/reference/FULL_CONTACT_MANIFEST.json", "boundary/SCALAR_SHELL_BOUNDARY_SYSTEM.md", "inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md"),
        ("inputs/reference/FULL_CONTACT_MANIFEST.json", "astra-bulk/COMPLETE_L4_SCALAR_BULK.md", "inputs/reference/COMPLETE_L4_SCALAR_BULK.md"),
        ("inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json", "INERTIAL_SHELL_WAVEFORM.md", "inputs/prior/INERTIAL_SHELL_WAVEFORM.md"),
        ("inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json", "computed/INERTIAL_WAVEFORM_CELLS.csv", "inputs/prior/INERTIAL_WAVEFORM_CELLS.csv"),
        ("inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json", "computed/INERTIAL_WAVEFORM_ENCLOSURE.json", "inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json"),
        ("inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json", "archive/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json", "inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json"),
    ]:
        manifest = json.loads((root / manifest_rel).read_text())
        matching = [x for x in manifest["files"] if x["path"] == member]
        gate("manifest_membership:"+member, len(matching) == 1 and
             matching[0]["sha256"] == PINS[local][1] and matching[0]["bytes"] == PINS[local][0])

    a,p,s,u,v,uy,vy,e = [Poly.var(i) for i in range(8)]
    normal_second = 5 if inject == "normal_quadratic" else 6
    ninv = [Poly(1), 2*u, normal_second*u*u]
    squared_normal = jet_mul(jet_mul(ninv, ninv), [Poly(1), -4*u, Poly()])
    gate("unit_normal_series", squared_normal == [Poly(1), Poly(), Poly()],
         "(1-4 epsilon u) times the proposed inverse-normal square must equal one through order two")
    kjet = jet_mul([a, uy+2*a*u, Poly()], ninv)
    pjet = jet_mul([p, vy, Poly()], ninv)
    d2k = 2*kjet[2]
    d2p = (1 if inject == "taylor_half" else 2)*pjet[2]
    gate("derivative_taylor_coefficient", d2p == 4*u*vy+12*p*u*u)
    gate("extrinsic_derivative_coefficient", d2k == 4*u*uy+20*a*u*u)
    first_rows = {5: -2*a*u-p*v*F(1,3), 6: -2*p*u-s*v*F(1,2)}
    metric_material = -F(4,3)*p*u*v + F(1,6)*s*v*v
    scalar_material = (-4 if inject == "one_copy_doubled" else -2)*v*v
    scalar_coordinate = scalar_material-d2p.subs(first_rows)
    metric_coordinate = metric_material-d2k.subs(first_rows)
    gate("one_copy_scalar_source", scalar_coordinate == -4*p*u*u+2*s*u*v-2*v*v)
    metric_expected = s*v*v*F(1,6) + (12 if inject == "sk_metric_sign" else -12)*a*u*u
    gate("metric_coordinate_source", metric_coordinate == metric_expected)
    c_y = 32*u*first_rows[5]-F(8,3)*v*first_rows[6]
    en = c_y-4*metric_coordinate
    gate("robin_join", en == -16*a*u*u-F(16,3)*p*u*v+F(2,3)*s*v*v)
    # Original bulk note (21),(23): derivative TF source equals
    # 12u² Hess(Y²) + (-12u²-2v²) (dY tensor dY)^TF.
    source_tf_l2 = 12*u*u*F(1,2)+(-12*u*u-2*v*v)*F(1,2)
    c2 = 2*source_tf_l2
    gate("degree_two_algebraic_row", c2 == -2*v*v)
    en2 = -4*v*first_rows[6]-2*metric_coordinate
    gate("degree_two_robin_join", en2 == 24*a*u*u+8*p*u*v+F(5,3)*s*v*v and en2 != en*F(1,2))

    n = [Poly.var(i) for i in range(4)]
    r2 = sum((x*x for x in n), Poly())
    angular_y = n[0]*n[0]-n[1]*n[1]
    j = n[0]*n[0]+n[1]*n[1]
    f2 = j-F(1,2)
    f4 = angular_y**2-j*F(1,2)+F(1,12)
    gate("harmonic_decomposition", angular_y**2 == F(1,6)+f2*F(1,2)+f4)
    homogeneous = [j-r2*F(1,2), angular_y**2-j*r2*F(1,2)+r2*r2*F(1,12)]
    gate("harmonic_extensions", all(sum((z.diff(i).diff(i) for i in range(4)), Poly()) == 0 for z in homogeneous))
    pieces = [Poly(F(1,6)), f2*F(1,2), f4]
    gate("harmonic_orthogonality", all(sphere_mean(pieces[i]*pieces[k]) == 0 for i in range(3) for k in range(i)))
    piece_means = [sphere_mean(z*z) for z in pieces]
    total_mean = sphere_mean(angular_y**4)
    gate("harmonic_piece_norms", piece_means == [F(1,36), F(1,48), F(19,720)])
    angular_coefficient = F(3,40) if inject == "angular_norm" else F(3,20)
    gate("angular_total_norm", angular_coefficient == 2*total_mean and sum(piece_means) == total_mean,
         "ordinary S3 measure has area 2 pi squared; mean-normalized measure does not")
    volume_rho_power = 3 if inject == "volume_power" else 4
    gate("proper_volume_power", volume_rho_power == 3+1,
         "spatial volume supplies rho cubed and proper time supplies one more rho")
    ep, mup, gammap, hp, yp, yxxp, tanhsq, unused = [Poly.var(i) for i in range(8)]
    reduced_beta = (2 if inject == "uncancelled_gamma" else 1)*gammap + mup*(F(5,2)-F(35,4)*tanhsq)
    boundary_residual = -ep*hp-gammap*ep**2*yp+mup*yxxp+ep**2*reduced_beta*yp
    target_residual = mup*yxxp-ep*hp+mup*ep**2*(F(5,2)-F(35,4)*tanhsq)*yp
    gate("dynamic_boundary_gamma_cancellation", boundary_residual == target_residual)

    report = json.loads((root / "inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json").read_text())
    shell = json.loads((root / "inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json").read_text())
    ell, L = F(report["state"]["ell"]), F(3)
    mu = box(report["shell_parameter_intervals"]["mu"])
    sb = lambda key: box(shell["boxes"][key]["exact"])
    invr2, rho, ab, pb, sigma2 = [sb(k) for k in ("H_squared", "rho_b", "a_b", "p_b", "sigma_second_b")]
    gate("registered_state", ell == F(1,10**6) and report["state"]["domain"] == [0,3])
    gate("shell_signs", 0 < mu[0] <= mu[1] and 0 < ab[0] <= ab[1] and pb[0] <= pb[1] < 0 and 0 < sigma2[0] <= sigma2[1])
    exact_mu = sb("mu")
    gate("rounded_mu_encloses_exact_box", mu[0] <= exact_mu[0] <= exact_mu[1] <= mu[1])
    dv, dt = report["arithmetic"]["value_denominator"], report["cell_data"]["time_grid_denominator"]
    cells = []
    with (root / "inputs/prior/INERTIAL_WAVEFORM_CELLS.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            cells.append((F(row["left_grid"])/dt, F(row["right_grid"])/dt,
                          F(row["H_lo"])/dv, F(row["H_hi"])/dv,
                          F(row["ell_S_lo"])/dv/ell, F(row["ell_S_hi"])/dv/ell,
                          F(row["H0_lo"])/dv, F(row["H0_hi"])/dv))
    gate("whole_cell_partition", len(cells) == 4353 and cells[0][0] == 0 and cells[-1][1] == L and
         all(z[0] < z[1] and z[2] <= z[3] and z[4] <= z[5] for z in cells) and
         all(cells[i][1] == cells[i+1][0] for i in range(len(cells)-1)))
    gate("zero_prefront", cells[0][:6] == (F(0), F(1), F(0), F(0), F(0), F(0)))
    weight = lambda left,right,x: x*(right-left)-(right*right-left*left)/2
    actual_abs_integral = sum(max(abs(lo),abs(hi))*weight(left,right,L) for left,right,lo,hi,*_ in cells)
    selected = (6,7) if inject == "use_comparator" else (2,3)
    abs_integral = sum(max(abs(z[selected[0]]),abs(z[selected[1]]))*weight(z[0],z[1],L) for z in cells)
    gate("actual_trace_envelope", abs_integral >= actual_abs_integral,
         "the comparator alone does not retain the actual-H whole-cell error enlargement")
    force = ell*ell/mu[0] if inject == "scaled_force" else ell/mu[0]
    gate("scaled_boundary_force", force*mu[0] == ell,
         "y_xx=ell squared times w_tt and H=ell(w_z-gamma*w) leave one ell")
    M = F(5,2)*ell*ell
    q = M*L*L/2
    gate("window_operator_bound", 9*ell*ell < F(4,7) and q == F(9,800000000000))
    yfree = force*abs_integral
    by = yfree if inject == "drop_resolvent" else yfree/(1-q)
    gate("volterra_resolvent", by >= yfree+q*by,
         "the proposed uniform bound must absorb the Volterra feedback term")
    error = q*by
    bu = invr2[1]*by
    radius_scale = 1/rho[0] if inject == "metric_radius" else 1/rho[0]**2
    gate("metric_reconstruction_radius", radius_scale*rho[0]**2 == 1,
         "u=V/rho squared and alpha(0)=1 in the retained representative")
    abs_h_integral = sum(max(abs(lo),abs(hi))*(right-left) for left,right,lo,hi,*_ in cells)
    byprime = ell/mu[0]*abs_h_integral + M*L*by
    chain = 1/ell if inject == "proper_derivative" else 1/(rho[0]*ell)
    gate("proper_clock_chain_rule", chain*rho[0]*ell == 1)
    butau = invr2[1]*chain*(byprime+F(3,2)*ell*ell*L*by)

    # Independent central envelope: split the whole history by signs and the
    # two bracket ends; no polynomial-prefix implementation is used here.
    A,B = box(report["zero_existence_brackets"][0])
    lower_integral = F(0)
    for left,right,lo,hi,*_ in cells:
        before_right = min(right,A)
        if left < before_right:
            lower_integral += lo*weight(left,before_right,A if lo >= 0 else B)
        middle_left, middle_right = max(left,A), min(right,B)
        if middle_left < middle_right:
            lower_integral += min(F(0),lo)*weight(middle_left,middle_right,B)
    gate("central_positive_convolution", lower_integral > 0)
    central_scale = ell/mu[0] if inject == "central_mu" else ell/mu[1]
    gate("positive_central_inverse_mu", central_scale*lower_integral <= ell/mu[1]*lower_integral,
         "a positive lower convolution bound must use mu maximum")
    ymin = central_scale*lower_integral-error
    umin = invr2[0]*(1-(ell*B)**2)*ymin
    gate("central_metric_positive", ymin > F(843,10**7) and umin > 0)
    witnesses = report["sign_witnesses"]
    gate("central_actual_scalar_zero", F(witnesses["central_left"]["x"]) == A and F(witnesses["central_right"]["x"]) == B and
         F(witnesses["central_left"]["H"][0]) > 0 and F(witnesses["central_right"]["H"][1]) < 0)
    sp_at_zero_min = (0 if inject == "erase_central_metric" else -4*pb[1])*umin**2
    sk_at_zero_max = -12*ab[0]*umin**2
    gate("both_central_source_signs", sp_at_zero_min > 0 and sk_at_zero_max < 0)
    endpoint_bounds = {}
    for name in ("central_left", "central_right"):
        witness = witnesses[name]
        vinterval = tuple(F(z)/ell for z in witness["ell_S"])
        gate("central_endpoint_nonzero:"+name, vinterval[0]*vinterval[1] > 0)
        vabsmin, vabsmax = min(map(abs,vinterval)), max(map(abs,vinterval))
        dominant_square = 0 if inject == "endpoint_drop_square" else -2*vabsmin*vabsmin
        sp_upper = -4*pb[0]*bu*bu+2*sigma2[1]*bu*vabsmax+dominant_square
        sk_lower = sigma2[0]*vabsmin*vabsmin/6-12*ab[1]*bu*bu
        gate("central_endpoint_source_signs:"+name, sp_upper < 0 and sk_lower > 0)
        endpoint_bounds[name] = {"x":F(witness["x"]),"actual_H":box(witness["H"]),
                                 "retained_ell_S":box(witness["ell_S"]),"v":vinterval,
                                 "Sp_upper":sp_upper,"SK_lower":sk_lower}

    tail = box(report["tail_theorem"]["actual_H_interval"])
    gain = box(report["physical_gain_interval"])
    gate("tail_cells_in_retained_theorem", 0 < tail[0] <= tail[1] and gain[0] <= gain[1] < 0 and
         all(tail[0] <= z[2] <= z[3] <= tail[1] for z in cells if z[0] >= 2))
    vmin, vmax = -gain[1]/ell*(1-(ell*L)**2)*tail[0], -gain[0]/ell*tail[1]
    negative_square = 0 if inject == "tail_drop_square" else -2*vmin*vmin
    sp_tail_upper = -4*pb[0]*bu**2+2*sigma2[1]*bu*vmax+negative_square
    sk_tail_lower = sigma2[0]*vmin*vmin/6-12*ab[1]*bu**2
    gate("strict_tail_scalar_sign", sp_tail_upper < F("-2.0579e-7"))
    gate("strict_tail_metric_sign", sk_tail_lower > F("6.8724e-8"))
    bv = max(max(abs(z[4]),abs(z[5])) for z in cells)
    bsp = -4*pb[0]*bu**2+2*sigma2[1]*bu*bv+2*bv*bv
    bsk = sigma2[1]*bv*bv/6+12*ab[1]*bu**2
    det_min = 24*ab[0]*sigma2[0]-16*max(abs(pb[0]),abs(pb[1]))**2
    gate("registered_robin_symmetrizer_positive", 6*ab[0] > 0 and det_min > 0)

    # This join checks the new report's norm factors and ratio inference.
    # It does not independently regenerate its complete-cell quadratures;
    # the separate Sol checker does that independent task.
    new_report = json.loads((root / "computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json").read_text())
    new = new_report["exact"]
    source_table = (root / "computed/COORDINATE_SHELL_SOURCE_CELLS.csv").read_bytes()
    source_sha = hashlib.sha256(source_table).hexdigest()
    gate("new_source_table_identity", source_sha == "d7ac199c89109c360c5da52b6d8693f4d34950f314fc42a5bde35a39d71859c5" and
         source_sha == new_report["cell_data"]["sha256"] and len(source_table) == new_report["cell_data"]["bytes"])
    def arctan_enclosure(denominator, terms):
        partial = sum((F((-1)**k,(2*k+1)*denominator**(2*k+1)) for k in range(terms)), F(0))
        next_partial = partial+F((-1)**terms,(2*terms+1)*denominator**(2*terms+1))
        return min(partial,next_partial),max(partial,next_partial)
    atan5, atan239 = arctan_enclosure(5,70), arctan_enclosure(239,20)
    independent_pi = (16*atan5[0]-4*atan239[1],16*atan5[1]-4*atan239[0])
    pi = box(new["pi"])
    gate("norm_pi_interval", pi[0] <= independent_pi[0] <= independent_pi[1] <= pi[1])
    gate("norm_background_identity", box(new["rho"]) == rho and F(new["ell"]) == ell)
    pref = tuple(F(3,20)*pi[i]**2*rho[i]**4*ell for i in (0,1))
    enpref = tuple(F(19,360)*pi[i]**2*rho[i]**4*ell for i in (0,1))
    gate("norm_full_Y_squared_prefactor", box(new["L2_prefactor"]) == pref)
    gate("norm_EN_degree_four_prefactor", box(new["EN_L2_prefactor"]) == enpref)
    for name, record in new["norms"].items():
        factor = enpref if name in ("EN","dEN","EN_material") else pref
        field = "H4" if name in ("EN","dEN","EN_material") else "Y^2"
        integral = box(record["weighted_x_squared_integral"])
        norm_square = box(record["squared_L2"])
        gate("norm_factor_join:"+name, record["angular_field"] == field and
             norm_square == tuple(factor[i]*integral[i] for i in (0,1)))
    dnorm = box(new["norms"]["dSp"]["squared_L2"])
    matnorm = box(new["norms"]["Sp_material"]["squared_L2"])
    ratio = dnorm[1]/matnorm[1 if inject == "ratio_denominator" else 0]
    gate("norm_ratio_interval_upper", matnorm[0] > 0 and ratio >= dnorm[1]/matnorm[0])
    gate("norm_ratio_report_identity", ratio == F(new["metric_correction_squared_norm_ratio_upper"]["Sp"]))
    gate("norm_ratio_strict_nine_e_minus_thirteen", ratio < F(9,10**13)**2)
    norm_join = {"source_table_sha256":source_sha,"dSp_squared_norm":dnorm,"material_squared_norm":matnorm,
                 "squared_ratio_upper":ratio,"strict_unsquared_ratio_bound":F(9,10**13),
                 "scope":"Exact report factor and ratio inference; full quadratures independently reconstructed by the separate Sol checker."}

    return {"schema": "hdblast-independent-astra-coordinate-sources-v1", "status": "PASS",
            "positive_groups": groups, "positive_group_count": len(groups),
            "inputs": authenticated, "independent_code": "stdlib Fraction polynomial and interval sums; no imported producer or upstream scientific code",
            "source_coefficients": {"Sp": "-4*p*u^2+2*sigma_second*u*v-2*v^2", "SK": "sigma_second*v^2/6-12*a*u^2", "EN_degree_four_algebra": "-16*a*u^2-(16/3)*p*u*v+(2/3)*sigma_second*v^2", "EN_degree_two_algebra_only": "24*a*u^2+8*p*u*v+(5/3)*sigma_second*v^2"},
            "angular_mean_square_factors": packed({"degree_0":piece_means[0], "degree_2":piece_means[1], "degree_4":piece_means[2], "total":total_mean}),
            "ordinary_angular_integral_over_pi_squared": packed(2*total_mean),
            "proper_volume_radial_power": volume_rho_power,
            "global_bounds": packed({"abs_H_triangular_integral":abs_integral,"q":q,"y":by,"y_minus_F":error,"y_x":byprime,"u":bu,"d_tau_u":butau,"v":bv,"abs_Sp":bsp,"abs_SK":bsk}),
            "central_zero": packed({"bracket":(A,B),"signed_convolution_lower":lower_integral,"y_lower":ymin,"u_lower":umin,"Sp_lower_at_any_scalar_zero":sp_at_zero_min,"SK_upper_at_any_scalar_zero":sk_at_zero_max}),
            "central_endpoint_source_signs": packed(endpoint_bounds),
            "new_report_norm_join": packed(norm_join),
            "tail": packed({"domain":(F(2),F(3)),"abs_v":(vmin,vmax),"Sp_upper":sp_tail_upper,"SK_lower":sk_tail_lower}),
            "robin_symmetrizer_determinant_lower": packed(det_min),
            "claims": ["At least one actual central scalar zero is inherited from the pinned continuous waveform theorem.", "At every such zero Sp>0 and SK<0; at both central bracket endpoints and throughout the tail Sp<0 and SK>0.", "Continuity gives at least two distinct zeros of EACH coordinate coefficient inside the central scalar-zero bracket, one on either side of any actual scalar zero there. No uniqueness or count of all zeros is claimed.", "Y squared has angular zeros; a strict coefficient sign is not strict everywhere on the sphere."],
            "scope": "Full known one-copy coordinate quadratic boundary-source coefficients in the specified first longitudinal representative and derivative Taylor convention. These are not U2, physical forces, gauge-invariant sources, tensor energy norms, or a nonlinear remainder. Original analytic/PDE certificates are retained hypotheses."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inject-error", choices=sorted(CONTROLS))
    parser.add_argument("--run-controls", action="store_true")
    args = parser.parse_args()
    if args.run_controls:
        results = []
        for optimized in (False, True):
            for control, expected in CONTROLS.items():
                command = [sys.executable, "-I", "-B"] + (["-O"] if optimized else []) + [str(Path(__file__).resolve()), "--root", str(args.root.resolve()), "--inject-error", control]
                run = subprocess.run(command, text=True, capture_output=True, check=False)
                try:
                    receipt = json.loads(run.stdout)
                except json.JSONDecodeError as exc:
                    raise Failure("control_receipt", control+": "+run.stderr) from exc
                if run.returncode != 1 or receipt.get("failed_gate") != expected:
                    raise Failure("control_wrong_failure", control+": "+run.stdout+run.stderr)
                results.append({"control": control, "optimized": optimized, "exit_code": run.returncode, "failed_gate": expected, "status": "EXPECTED_REJECTION"})
        result = {"schema": "hdblast-astra-coordinate-controls-v1", "status": "PASS", "controls": results,
                  "rejected_processes": len(results), "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        code = 0
    else:
        try:
            result = audit(args.root.resolve(), args.inject_error)
            result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            code = 0
        except Failure as exc:
            result = {"status": "FAIL", "failed_gate": exc.gate, "detail": exc.detail, "injected_error": args.inject_error}
            code = 1
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
