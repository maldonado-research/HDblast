"""Test local energy accounting of an adiabatic-target bath on pinned B1 wall histories."""
import argparse
from bisect import bisect_left
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

WALL = Path("research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip")
CLASSICAL = Path("research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs")
WALL_HASH = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
EVENT_HASH = {
    "fine": "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810",
    "coarse": "91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c",
}
KEYS = ("H0tau", "A_b", "phi_b", "H_over_H0", "v_over_H0")
EDGES = (8, 12, 14)
CUT = 2.10
G = 100.0


def checked(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"Input checksum mismatch for {path}")
    return raw


def wall_rows(archive, name):
    with archive.open(name) as stream:
        rows = csv.DictReader(io.TextIOWrapper(stream, encoding="ascii"))
        return [{key: float(row[key]) for key in KEYS} for row in rows]


def cut_at(rows, cut, count=2):
    times = [row["H0tau"] for row in rows]
    if any(a >= b for a, b in zip(times, times[1:])) or not times[0] < cut <= times[-1]:
        raise ValueError("Invalid proper clock")
    right = bisect_left(times, cut)
    if count == 2:
        lo, hi = rows[right-1], rows[right]
        w = (cut-lo["H0tau"])/(hi["H0tau"]-lo["H0tau"])
        return {key: lo[key] + w*(hi[key]-lo[key]) for key in KEYS}
    if count != 4 or right < 2 or right+1 >= len(rows):
        raise ValueError("Invalid four-point interpolation")
    near = rows[right-2:right+2]
    weights = [math.prod((cut-other["H0tau"])/(r["H0tau"]-other["H0tau"])
                         for other in near if other is not r) for r in near]
    return {key: math.fsum(w*r[key] for r, w in zip(near, weights)) for key in KEYS}


def simpson(fn, end, steps):
    ds = end/steps
    odd = math.fsum(fn(i*ds) for i in range(1, steps, 2))
    even = math.fsum(fn(i*ds) for i in range(2, steps, 2))
    return ds*(fn(0) + 4*odd + 2*even + fn(end))/3


def modes(k, a, h, signed_mass, signed_rate):
    """Independent canonical-mode construction. Evaluate the two signed daughter rates."""
    omega = math.hypot(k, a*signed_mass)  # conformal frequency
    energy = omega/a
    omega_s = a*a*signed_mass*(h*signed_mass + signed_rate)/omega
    # One arbitrary common instantaneous phase suffices.
    u = complex(1/math.sqrt(2*omega))
    u_eta = complex(-a*omega_s/(2*omega), -omega)*u
    chi = u/a
    momentum = a*u_eta - a*a*h*u
    g = a**3*energy
    if abs((chi*momentum.conjugate()).imag - .5) > 2e-14:
        raise ValueError("Canonical Wronskian check failed")
    # Exact B-b overlap from the negative-frequency component; robust near the vacuum.
    beta = (math.sqrt(g/2)*chi.conjugate()
            + 1j*momentum.conjugate()/math.sqrt(2*g))
    excited = abs(beta)**2
    # b-vacuum mode has chi_b=1/sqrt(2g), momentum_b=-i sqrt(g/2).
    # Physical parent energy gained per unit damping from the covariance difference.
    d_kinetic = .5/a**3*(abs(momentum)**2-g/2)
    d_potential = .5*a**3*energy**2*(abs(chi)**2-1/(2*g))
    gamma = signed_mass**2/(8*math.pi*energy)
    measure = k*k/(2*math.pi**2*a**3)
    return (-measure*gamma*(d_kinetic+d_potential), measure*gamma*energy*excited)


def point(row, astar, verify=True):
    a = math.exp(row["A_b"]-astar)
    h = row["H_over_H0"]
    signed_mass = G*(row["phi_b"]-.5)
    signed_rate = G*row["v_over_H0"]
    mass = abs(signed_mass)
    if mass < 1e-6:
        raise ValueError("This comparison requires a massive endpoint")
    mix = (h*mass**2 + signed_mass*signed_rate)/2
    coef = mass**2/(64*math.pi**3)

    def density(p):
        e2 = p*p + mass*mass
        return coef*p*p*(h + mix/e2)**2/e2

    def primitive(P):
        q = P*P+mass*mass
        t = math.atan(P/mass)
        first = P-mass*t
        second = (t/mass-P/q)/2
        third = (t/mass**3+P/(mass**2*q))/8-P/(4*q*q)
        return coef*math.fsum((h*h*first, 2*h*mix*second, mix*mix*third))

    result = {}
    for edge in EDGES:
        pmax = edge/a
        exact = primitive(pmax)
        if not verify:
            result[str(edge)] = exact
            continue
        coarse, fine = (simpson(density, pmax, n) for n in (128, 256))
        if exact <= 0 or not math.isclose(fine, exact, rel_tol=3e-10, abs_tol=1e-14):
            raise ValueError("Power or analytic integration check failed")
        signed = [simpson(lambda k: modes(k, a, h, signed_mass, signed_rate)[i],
                          edge, 256) for i in (0, 1)]
        for calculated, expected in zip(signed, (-exact,exact)):
            if not math.isclose(calculated,expected,rel_tol=2e-9):
                raise ValueError("Canonical covariance or volume check failed")
        result[str(edge)] = dict(physical_p_max=pmax, physical_daughter_power=-exact,
                                 assigned_daughter_power=exact, gap=2*exact,
                                 simpson128_error=coarse-exact, simpson256_error=fine-exact,
                                 canonical_trace_error=signed[0]+exact,
                                 canonical_count_error=signed[1]-exact,
                                 max_physical_frequency=math.hypot(pmax,mass))
    if not verify:
        return result

    # Massive IR: check the leading p³ integral at a scale where the elementary
    # primitive loses precision from cancellation.
    ir_p = mass*1e-4
    ir_lead = coef*(h + mix/mass**2)**2/mass**2 * ir_p**3/3
    ir_checked = simpson(density,ir_p,32)
    if not math.isclose(ir_lead,ir_checked,rel_tol=1e-6):
        raise ValueError("Infrared scaling check failed")
    return dict(proper_time=row["H0tau"], a=a, phi=row["phi_b"], phi_dot=row["v_over_H0"],
                h=h, mass=mass, mass_squared=mass**2, mass_squared_dot=2*signed_mass*signed_rate,
                delta_ir_limit=h+mix/mass**2, ir_leading=ir_lead, ir_integrated=ir_checked,
                bands=result)


def compute(repo):
    reference = {}
    for grid, pin in EVENT_HASH.items():
        path = CLASSICAL/grid/"quiet/rolling_summary.json"
        reference[grid] = json.loads(checked(repo/path,pin))["matter"]["events"][0]["A_star"]
    out = {}
    raw = checked(repo/WALL,WALL_HASH)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        for grid, astar in reference.items():
            for model in ("y0_b0","y0_b1","y1_b0","y1_b1"):
                name = f"{grid}/{model}"
                prep = point(wall_rows(archive,name+"/pre.csv")[0],astar)
                arms = {}
                for arm in ("quiet","plus","minus","far"):
                    rows = wall_rows(archive,f"{name}/{arm}.csv")
                    reading = point(cut_at(rows,CUT),astar)
                    cubic = point(cut_at(rows,CUT,count=4),astar,verify=False)
                    for edge, band in reading["bands"].items():
                        band["cubic_minus_linear"] = cubic[edge] - band["assigned_daughter_power"]
                    arms[arm] = reading
                out[name] = dict(preparation=prep,cut=arms)
    return dict(wall_sha256=WALL_HASH, event_sha256=EVENT_HASH, log_a_star=reference,
                bath_y=1.0, G=G, phi_star=.5, upper_comoving_edges=EDGES, results=out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository",type=Path)
    args = parser.parse_args()
    print(json.dumps(compute(args.repository),indent=2,allow_nan=False))
