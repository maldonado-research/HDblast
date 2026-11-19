"""Retarded response of the published Gaussian spectral continuum to two positive probes."""
import argparse
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicHermiteSpline
from scipy.special import roots_legendre

PI = math.pi
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
PHASES = (4, 8, 16, 32, 64)
CUTOFFS = (60, 120, 240, 480)


def native_probe(repo):
    native = repo / NATIVE
    fork = json.loads((native / "FORK.json").read_text())["fine"]
    with zipfile.ZipFile(native / "WALLS.zip") as archive:
        history = np.genfromtxt(io.BytesIO(archive.read("fine/retain/quiet/wall.csv")), delimiter=",", names=True)
    loga = CubicHermiteSpline(history["H0tau"], history["A_b"] - fork["reference_A"], history["H_over_H0"])
    knots = loga.x[(loga.x > 1.35) & (loga.x < 2.05)]
    width = quad(lambda s: math.exp(-float(loga(s))), 1.35, 2.05, points=knots,
                 limit=len(knots) + 50, epsabs=2e-13)[0]
    return float(min(fork["q"], key=lambda q: abs(q - 2.0))), width


def transform(w, width):
    """Real centered Fourier amplitude of cos^4(pi t/T) on [-T/2,T/2]."""
    x = w * width / 2
    if abs(x) < 3 * PI:
        return width * (3/8*np.sinc(x/PI) + (np.sinc(x/PI-1)+np.sinc(x/PI+1))/4
                        + (np.sinc(x/PI-2)+np.sinc(x/PI+2))/16)
    return width * (1.5*PI**4) * math.sin(x) / (x*(x*x-PI*PI)*(x*x-4*PI*PI))


def continuum(q, width, v):
    # The Abel-regulated retarded kernel is ordinary and smooth at every lag in this integral.
    def integrand(t):
        lag = v - t
        kernel = -(q*lag*math.sin(q*lag)+math.cos(q*lag))/(2*PI**2*lag**3)
        return math.cos(PI*t/width)**4 * kernel
    return quad(integrand, -width/2, width/2, epsabs=1e-16, epsrel=3e-13)[0]


def band(q, width, v, cutoff, smooth):
    def spectral(w):
        return (w*w-q*q)/(4*PI**2)*transform(w,width)
    def oscillatory(fn, left, right):
        value, error = quad(fn,left,right,weight="sin",wvar=v,epsabs=3e-15,limit=400)
        return value
    if not smooth:
        return oscillatory(spectral,q,cutoff)
    edge = .8*cutoff
    def tapered(w):
        x = (w-edge)/(cutoff-edge)
        return spectral(w)*(1-x)**3*(1+3*x+6*x*x)
    return oscillatory(spectral,q,edge) + oscillatory(tapered,edge,cutoff)


def discrete(q,width,v,cutoff,order):
    nodes, weights = roots_legendre(order)
    w = q + (cutoff-q)*(nodes+1)/2
    coefficients = np.array([transform(float(k),width) for k in w])
    return float((cutoff-q)/2 * np.dot(weights, (w*w-q*q)/(4*PI**2)*coefficients*np.sin(w*v)))


def combine(fn,q,v):
    last = fn(v)
    return {"single":last, "half":last+fn(v+PI/q), "full":last+fn(v+2*PI/q)}


def experiment(repo):
    q,width = native_probe(repo)
    if PI/q <= width:
        raise ValueError("The declared probe pulses overlap")
    threshold = transform(q,width)
    normalization = -q*threshold/(2*PI**2)
    rows = []
    for n in PHASES:
        v = (2*PI*n+PI/2)/q
        raw = combine(lambda t:continuum(q,width,t),q,v)
        scaled = dict(single=raw["single"]*v**2/normalization,
                      half=raw["half"]*v**3/(2*PI/q*normalization),
                      full=raw["full"]*v**2/(2*normalization))
        cutoffs = {}
        for cutoff in CUTOFFS:
            cutoffs[str(cutoff)] = {}
            for name, smooth in (("sharp",False),("tapered",True)):
                vals = combine(lambda t:band(q,width,t,cutoff,smooth),q,v)
                cutoffs[str(cutoff)][name] = {k:vals[k]/raw[k]-1 for k in raw}
        rows.append({"n":n,"v":v,"continuum":raw,"scaled":scaled,
                     "half_to_full":raw["half"]/raw["full"],"cutoff_fractional_errors":cutoffs})
    # Compare the published discrete quadrature to the evaluated sharp integral on an early and late phase.
    convergence = {}
    for n in (1,8):
        v = (2*PI*n+PI/2)/q
        target = band(q,width,v,120,False)
        convergence[str(n)] = {str(order):discrete(q,width,v,120,order)/target-1
                               for order in (64,256,1024,4096)}
    if max(abs(x-1) for x in rows[-1]["scaled"].values()) > .025:
        raise ValueError("Phase-resolved asymptotic amplitudes")
    if not all(a["half_to_full"] > b["half_to_full"] for a,b in zip(rows,rows[1:])):
        raise ValueError("Loss of the planned half-period cancellation")
    for row in rows:
        if max(abs(v) for v in row["cutoff_fractional_errors"]["480"]["tapered"].values()) > 3e-5:
            raise ValueError("Independent finite-band integral disagrees with causal force")
    if max(abs(case["4096"]) for case in convergence.values()) > 1e-5:
        raise ValueError("Unresolved quadrature control")
    return dict(q=q,conformal_pulse_width=width,separation_half=PI/q,threshold_form_factor=threshold,
                phase_clock="time from last pulse center",rows=rows,discrete_sharp_120_fractional_errors=convergence)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(experiment(args.repo),indent=2,allow_nan=False)+"\n")
