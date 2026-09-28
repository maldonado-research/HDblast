"""Shared definitions for the mechanism screen (read-only use of earlier packages).

Registered model in dimensionless units (kappa_5^2 = 1, model length unit L):
  W = 1 - phi + phi^3/3,  U = W_phi^2/2 - (2/3) W^2,  sigma = 2W + delta(1 + c phi).
AdS radii of the two bulk vacua: l(+1) = 3/W(1) = 9, l(-1) = 3/W(-1) = 9/5.
"""
from pathlib import Path
import hashlib, json, math

C = 2/1.0357712571566784 - 4/3          # registered c
DELTA = 1e-3                            # registered detuning

ROOT = Path('/home/user/unified-theory-maldonado')
PKG = ROOT/'new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922'
CHAT14_ZIP = PKG/'source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip'
PLUS_JSON = PKG/'static_branch/PLUS_BRANCH_RESULTS.json'
CHAT14_PREFIX = 'HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/'
HERE = Path(__file__).resolve().parent

def W(p): return 1 - p + p**3/3
def Wp(p): return p*p - 1
def U(p): return 0.5*(p*p-1)**2 - (2/3)*W(p)**2
def Up(p): return (p*p-1)*2*p - (4/3)*W(p)*(p*p-1)
def sigma(p, d=DELTA, c=C): return 2*W(p) + d*(1 + c*p)
def sigma1(p, d=DELTA, c=C): return 2*(p*p-1) + d*c

# ---- observational / physical constants (values quoted in README with sources) ----
HBARC_EV_M = 1.973269804e-7           # eV m
MPL_RED_EV = 2.435e27                 # reduced Planck mass, eV
H0_TODAY_EV = 2.1332e-33*0.674        # 100 h km/s/Mpc = 2.1332e-33 h eV, h = 0.674
OMEGA_L = 0.685
SEC_PER_EV_INV = 6.582119569e-16      # hbar in eV s
GYR_S = 3.15576e16

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(name, obj):
    def clean(v):
        if isinstance(v, dict): return {k: clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, float) and not math.isfinite(v): return None
        if hasattr(v, 'item'): return clean(v.item())
        return v
    (HERE/name).write_text(json.dumps(clean(obj), indent=2) + '\n')
