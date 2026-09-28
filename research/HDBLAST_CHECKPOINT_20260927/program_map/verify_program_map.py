#!/usr/bin/env python3
"""Source-consistency verifier for PROGRAM_MAP.md (HDBLAST program map, 27 Sept 2026).

This script does NOT produce new physics. It
  1. recomputes the small set of numbers that PROGRAM_MAP.md quotes and that can be
     recomputed cheaply from the model definition or from saved JSON outputs of earlier
     checkpoints (closed forms, expansions, ratios, interval containments);
  2. checks that every dated checkpoint in the chronology is anchored by a date string
     (or dated folder name) inside a named source file;
  3. rebuilds the Zenodo version chronology of the HDBLAST concept record from the
     29 Aug 2026 provenance audit table and locates the later records in named files;
  4. checks that every absolute path quoted in PROGRAM_MAP.md exists;
  5. runs deliberate controls (wrong formulas, perturbed parameters, a fabricated record
     number, a fabricated path) that MUST fail.

Read-only with respect to all source folders. Writes only program_map_checks.json next
to this script. Runtime: a few seconds, one core.
"""
from pathlib import Path
import hashlib, json, math, re, sys, time

import numpy as np
from scipy.optimize import brentq
from scipy.special import kv

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/user/unified-theory-maldonado/new-files')
DB3 = ROOT / 'D-Blast 3'
LATEST = ROOT / 'latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922'
MAP = HERE / 'PROGRAM_MAP.md'
OUT = HERE / 'program_map_checks.json'

checks = []


def check(name, ok, **detail):
    checks.append(dict(name=name, passed=bool(ok), **detail))
    return bool(ok)


def control(name, should_fail_condition, **detail):
    """A control passes when the deliberately wrong input is DETECTED (condition True)."""
    checks.append(dict(name='CONTROL: ' + name, passed=bool(should_fail_condition), control=True, **detail))


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def loglog_slope(x, y):
    x = np.log(np.asarray(x, float)); y = np.log(np.abs(np.asarray(y, float)))
    return float(np.polyfit(x, y, 1)[0])


t0 = time.time()
result = dict(script='verify_program_map.py', started_unix=t0,
              python=sys.version.split()[0], numpy=np.__version__)

# ---------------------------------------------------------------------------
# 1. Model constants (registered model, frozen solver)
# ---------------------------------------------------------------------------
solver_src = (LATEST / 'frozen/registered_solver.py').read_text()
IP = float(re.search(r'^IP=([0-9.eE+-]+)', solver_src, re.M).group(1))
c = 2 / IP - 4 / 3
c_quoted = 0.5975949350280132
check('c = 2/I_+ - 4/3 reproduces registered c', abs(c - c_quoted) < 1e-15,
      I_plus=IP, c_recomputed=c, c_quoted=c_quoted, abs_diff=abs(c - c_quoted),
      source=str(LATEST / 'frozen/registered_solver.py'))

# Chat 9 closed-form tachyon mass and O(t) slope; certified interval; growth rate
mu2_0 = -4 * (3 * c * c - 4 * c + 8) / (c * (3 * c + 4))
slope = 1.9243896  # Chat 9 analytic slope, quoted
mu2_t = mu2_0 + slope * 1e-3
cert_lo, cert_hi = -7.71788, -7.71786  # Chat 9 interval certificate (quoted)
mu2_numeric = -7.7178716  # Chat 9 shooting value (quoted)
check('Chat 9 closed-form mu^2(t->0) = -7.719796', abs(mu2_0 - (-7.719796)) < 5e-7, mu2_closed_form=mu2_0)
check('Chat 9 closed form + O(t) slope lies inside certified interval at t=1e-3',
      cert_lo < mu2_t < cert_hi, mu2_linear_in_t=mu2_t, certified_interval=[cert_lo, cert_hi])


def growth(mu2):
    """Growth rate s (units of H) of a scalar mode of mass^2 mu2*H^2 on de Sitter 4: s^2+3s+mu2=0."""
    return (-3 + math.sqrt(9 - 4 * mu2)) / 2


s = growth(mu2_numeric)
check('growth rate from mu^2=-7.7178716 equals 1.65719', abs(s - 1.65719) < 5e-6, growth_rate=s)
control('wrong growth formula s=sqrt(-mu^2) is rejected', abs(math.sqrt(-mu2_numeric) - 1.65719) > 0.1,
        wrong_value=math.sqrt(-mu2_numeric))
# e-folding time quoted as 0.60 Hubble times
check('e-folding time 1/s ~ 0.60 H^-1', abs(1 / s - 0.603) < 0.005, efold_time=1 / s)

# Registered-shell controls (folder 152): quoted rates and their difference
r_prev, r_eig, r_evo = 1.657193631245, 1.657193663767, 1.657194158956
check('folder-152 eigenvalue vs previous rate differ by ~3.25e-8', abs((r_eig - r_prev) - 3.25e-8) < 1e-10,
      diff=r_eig - r_prev, source=str(DB3 / 'untitled folder 152/00_READ_FIRST.md'))

# ---------------------------------------------------------------------------
# 2. New static +1 branch (22 Sept 2026 package)
# ---------------------------------------------------------------------------
pb = json.loads((LATEST / 'static_branch/PLUS_BRANCH_RESULTS.json').read_text())
rows = pb['rows']
deltas = np.array([r['delta'] for r in rows])
phib = np.array([r['phi_b'] for r in rows])
rhob = np.array([r['rho_b'] for r in rows])
H2 = np.array([r['H2'] for r in rows])
H2_metric = np.array([r['metric_only_H2'] for r in rows])
cc = pb['c']
reg = [r for r in rows if abs(r['delta'] - 1e-3) < 1e-15][0]
check('registered branch values match quoted phi_b, rho_b, H^2',
      abs(reg['phi_b'] - 0.9999159473169134) < 1e-15 and abs(reg['rho_b'] - 129.9247628497) < 1e-9
      and abs(reg['H2'] - 5.92401479433e-5) < 1e-15,
      phi_b=reg['phi_b'], rho_b=reg['rho_b'], H2=reg['H2'], eta_h=reg['eta_h'],
      source=str(LATEST / 'static_branch/PLUS_BRANCH_RESULTS.json'))
check('de Sitter identity H^2 * rho_b^2 = 1 on all six detunings',
      float(np.max(np.abs(H2 * rhob ** 2 - 1))) < 1e-10, max_abs_dev=float(np.max(np.abs(H2 * rhob ** 2 - 1))))

# metric-only (constant phi=1) Randall-Sundrum benchmark: H^2 = (sigma/6)^2 - (W(1)/3)^2
H2_rs = deltas * (1 + cc) / 27 + deltas ** 2 * (1 + cc) ** 2 / 36
check('metric-only H^2 equals exact RS formula delta(1+c)/27 + delta^2(1+c)^2/36',
      float(np.max(np.abs(H2_rs / H2_metric - 1))) < 1e-12,
      max_rel_dev=float(np.max(np.abs(H2_rs / H2_metric - 1))))

# small-detuning expansions quoted in the map
res_phi = phib - (1 - 9 * cc / 64 * deltas)
res_H2 = H2 - (deltas * (1 + cc) / 27 + deltas ** 2 * ((1 + cc) ** 2 / 36 - cc ** 2 / 384))
sl_phi = loglog_slope(deltas[:4], res_phi[:4])
sl_H2 = loglog_slope(deltas[:4], res_H2[:4])
check('phi_b = 1 - (9c/64) delta + O(delta^2): residual slope ~2', 1.8 < sl_phi < 2.2,
      loglog_slope=sl_phi, residuals=res_phi.tolist(), fit_points='delta<=0.01')
check('H^2 two-term expansion: residual slope ~3', 2.7 < sl_H2 < 3.3,
      loglog_slope=sl_H2, residuals=res_H2.tolist(), fit_points='delta<=0.01')
# controls: wrong coefficients must degrade the residual scaling
res_phi_bad = phib - (1 - 9 * cc / 32 * deltas)
res_H2_bad = H2 - (deltas * (1 + cc) / 27 + deltas ** 2 * ((1 + cc) ** 2 / 36 - cc ** 2 / 192))
control('wrong scalar coefficient 9c/32 gives O(delta) residual', loglog_slope(deltas[:4], res_phi_bad[:4]) < 1.3,
        slope=loglog_slope(deltas[:4], res_phi_bad[:4]))
control('wrong H^2 coefficient -c^2/192 gives O(delta^2) residual', loglog_slope(deltas[:4], res_H2_bad[:4]) < 2.3,
        slope=loglog_slope(deltas[:4], res_H2_bad[:4]))

# H/H0 relative to the original unstable shell; rho0 quoted in the independent review
review_txt = (LATEST / 'static_branch/INDEPENDENT_BRANCH_REVIEW.md').read_text()
rho0 = float(re.search(r'ρ₀=([0-9.]+)', review_txt).group(1))
ratio = rho0 / reg['rho_b']
ratio_metric_formula = math.sqrt(1e-3 * (1 + cc) / 27 + 1e-6 * (1 + cc) ** 2 / 36) * rho0
ppm = (ratio / 0.606726506 - 1) * 1e6
check('H/H0 = rho0/rho_b = 0.606721732', abs(ratio - 0.606721732) < 5e-10, rho0=rho0, H_over_H0=ratio)
check('RS metric benchmark H/H0 ~ 0.6067265', abs(ratio_metric_formula - 0.606726506) < 5e-8,
      H_over_H0_metric=ratio_metric_formula)
check('scalar-profile correction to H is about -7.87 ppm', abs(ppm - (-7.868927)) < 0.01, ppm=ppm)
control('perturbed c (+1%) moves the RS benchmark by >> 7.87 ppm',
        abs(math.sqrt(1e-3 * (1 + 1.01 * cc) / 27 + 1e-6 * (1 + 1.01 * cc) ** 2 / 36) * rho0 / 0.606726506 - 1) * 1e6 > 100)

# ---------------------------------------------------------------------------
# 3. Chat 14 registered-detuning roll-off (from the 22 Sept independent recompute)
# ---------------------------------------------------------------------------
c14 = json.loads((LATEST / 'source_audit/CHAT14_RECOMPUTED_RESULTS.json').read_text())
turn = {k: v['events']['turnaround'] for k, v in c14['runs'].items()
        if k.startswith('runs/v2_t1e3_minus_c4') and (v.get('events') or {}).get('turnaround')}
tt = [v['H0_tau'] for v in turn.values()]; tp = [v['phi_b'] for v in turn.values()]
check('Chat 14 turnaround H0*tau in [5.8944, 5.8947] across three c4 grids',
      len(tt) == 3 and 5.8943 < min(tt) and max(tt) < 5.8948, runs=sorted(turn), H0_tau=tt, phi_b=tp)
check('Chat 14 turnaround phi_b in [-1.9417, -1.9412]', len(tp) == 3 and -1.94175 < min(tp) and max(tp) < -1.94115)
plus69 = {k: (v.get('at_proper_time') or {}).get('6.9') for k, v in c14['runs'].items() if k.startswith('runs/v2_t1e3_plus_c4')}
hv = {k: round(v['H_over_H0'], 6) for k, v in plus69.items() if v}
check('Chat 14 refined-far-grid +1 run: H/H0 = 0.638642 at H0*tau = 6.9',
      abs(hv.get('runs/v2_t1e3_plus_c4_finecoarse', 0) - 0.638642) < 1e-6, values_at_6p9=hv)
check('Chat 14 archive integrity recheck passed (84 payloads)', c14['passed'] and c14['manifest_payload_count'] == 84)

# ---------------------------------------------------------------------------
# 4. Linear-response determinant intervals (v23 / v24 / numerics)
# ---------------------------------------------------------------------------
v23 = (2.19870426928051, 2.46151588838361)
v24 = (2.12278154373842557755059, 2.53743323308364125340774)
adj, fwd = 2.330107384842, 2.330107321657
check('v23 and v24 determinant intervals exclude zero', v23[0] > 0 and v24[0] > 0, v23=v23, v24=v24)
check('numerical adjoint and forward determinants lie in both intervals',
      all(v23[0] < x < v23[1] and v24[0] < x < v24[1] for x in (adj, fwd)),
      adjoint=adj, forward=fwd, rel_diff=abs(adj - fwd) / adj)
# Chat 8 global forward certificate (stored inside folder 144 inputs zip)
import zipfile
fz = DB3 / 'untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/inputs/PREVIOUS_GLOBAL_FORWARD.zip'
ftxt = zipfile.ZipFile(fz).read('00_READ_FIRST.md').decode()
m = re.search(r'\\boxed\{([0-9.]+)\\le\\det\\mathcal R\\le([0-9.]+)\.\}', ftxt)
gf = (float(m.group(1)), float(m.group(2))) if m else None
check('Chat 8 global forward interval [2.3290685, 2.3311463] parsed; positive; inside v23 and v24; contains numerics',
      gf is not None and gf[0] > 0 and v23[0] < gf[0] and gf[1] < v23[1] and v24[0] < gf[0] and gf[1] < v24[1]
      and gf[0] < adj < gf[1] and gf[0] < fwd < gf[1], forward_certificate=gf, source=str(fz) + '::00_READ_FIRST.md')

# ---------------------------------------------------------------------------
# 5. PTA knee convention (Contract v3, 22 Jan 2026): K2 transfer x50
# ---------------------------------------------------------------------------
def F(x, nu=2):
    return ((x ** 2 * kv(nu, x)) / 2) ** 2 if nu == 2 else None


x50 = brentq(lambda x: F(x) - 0.5, 0.1, 5, xtol=1e-14)
check('RS2 K2 filter F(x)=((x^2 K2(x))/2)^2 has x50 = 1.339139, ln x50 = 0.292027',
      abs(x50 - 1.339139) < 5e-7 and abs(math.log(x50) - 0.292027) < 5e-7, x50=x50, ln_x50=math.log(x50),
      source=str(DB3 / 'DBLAST FILES PURGE/Dblast 181/HDBLAST_fscale_vs_f50_clarification_20260122.md'))
# control: a K1-type filter (x K1(x))^2 must give a different half-suppression point
x50_k1 = brentq(lambda x: (x * kv(1, x)) ** 2 - 0.5, 0.05, 5, xtol=1e-14)
control('K1-type filter gives a different x50 (convention matters)', abs(x50_k1 - x50) > 0.05, x50_K1=x50_k1)

# ---------------------------------------------------------------------------
# 6. Chronology anchors: each dated entry must be backed by a string in a named file
# ---------------------------------------------------------------------------
anchors = [
    ('Brane-world two-link line starts (Zenodo 16866884)', DB3 / 'VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md', '16866884](https://zenodo.org/records/16866884) | 2025-08-14'),
    ('vLASTMILE+ canonical SBPL snapshot', DB3 / 'VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md', '17918095](https://zenodo.org/records/17918095) | 2025-12-12'),
    ('Kernel-free NANOGrav 15-yr certificate result', DB3 / 'ZENODO-PUBLICATIONS/Zenodo-18158842-jan-06-2026/HDBLAST_RESULTS_NANOGrav15yr_HD_20260106/README.md', 'Generated: 2026-01-06'),
    ('Theory Contract v3 (knee = f50, K2 filter)', DB3 / 'ZENODO-PUBLICATIONS/ZENODO-18341882-Version-01-22-26/ZENODO_DESCRIPTION_HDBLAST_DBLAST_20260122.md', '2026.01.22'),
    ('v7plus go/no-go memo: pilots FAIL/FAIL', DB3 / 'DBLAST FILES PURGE/Dblast 571/HDBLAST_GO_NO_GO_MEMO_20260307.md', 'All current pilots are **FAIL / FAIL**'),
    ('v7plus frozen gate criteria', DB3 / 'D-BLAST-PROGRESS-3-FOLDERS/untitled folder 89/HDBLAST_v7plus_Authoritative_LastMile_20260306/README_CURRENT_STATE_20260306.md', 'κ50 sanity: κ50_band_med < 4.0'),
    ('Status snapshot v8', DB3 / 'PURGE ARCHIVE APRIL 2026 FILES/untitled folder 160 ZENODO/HDBLAST_ZENODO_V8_UPLOAD_KIT_20260410/HDBLAST_ZENODO_V8_DESCRIPTION_20260410.md', 'v8 (2026-04-10)'),
    ('PHYS-M51 dated 2026-06-25 (earliest PHYS-M file name located)', DB3 / 'VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS.md', 'HDBLAST_PHYS_M51_COVARIANT_BULK_BRANE_TENSOR_TRANSFER_THEORY_ONLY_NO_V432_20260625'),
    ('Bessel-K v7 record 18807647 dated 2026-02-27', DB3 / 'VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md', '18807647](https://zenodo.org/records/18807647) | 2026-02-27'),
    ('Canonical SBPL knee parameters (secondary compilation)', DB3 / 'VECTOR STORE/HDBLAST VECTOR STORE/DBlast 2 vector store/HDBLAST_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md', '- \\(f_k=3.16\\times10^{-8}\\,{\\rm Hz}\\),'),
    ('PHYS-M309-M312 radiation/energy branch', DB3 / 'untitled folder 68/HDBLAST_NEW_CHAT_HANDOFF_PHYS_M285_M312_20260710.md', 'PHYS-M309'),
    ('Zenodo v20 correction (withdrawal)', DB3 / 'untitled folder 89 ZENODO/HDBLAST_ZENODO_V20_METADATA.txt', '2026-07-14'),
    ('v20 withdrawal text', DB3 / 'untitled folder 89 ZENODO/HDBLAST_ZENODO_V20_CORRECTION_AND_CLAIM_BOUNDARY.md', 'radiative Fourier-support certification withdrawn'),
    ('M462R1 registered shell root', DB3 / 'untitled folder 120/M462_MAIN_REPORT_RESEND.md', 'Date: 2026-08-17'),
    ('M468R1 / Zenodo v21 kit', DB3 / 'untitled folder 126 ZENODO/HDBLAST_ZENODO_V21_METADATA.txt', '2026-08-18'),
    ('M489G-A (v22 content)', DB3 / 'untitled folder 136/00_READ_FIRST_PHYS_M489GA_20260903.md', 'Date: 2026-09-03'),
    ('M489G-C rank-two determinant (v23 content)', DB3 / 'untitled folder 138 ZENODO/M489GC_ZENODO_RECOMMENDATION_AND_PLAIN_LANGUAGE_SUMMARY.md', '2.19870426928051 to 2.46151588838361'),
    ('Five-field forward/adjoint comparison', DB3 / 'untitled folder 138/00_READ_FIRST.md', 'Prepared 4 September 2026'),
    ('Continuum certification components', DB3 / 'untitled folder 139/00_READ_FIRST.md', 'Prepared 4 September 2026'),
    ('Conditional global adjoint (v24 content)', DB3 / 'untitled folder 140/00_READ_FIRST.md', 'Research checkpoint: 2026-09-04'),
    ('Forward error control components', DB3 / 'untitled folder 141/00_READ_FIRST.md', 'Research checkpoint, 5 September 2026'),
    ('Chat 8 fixed startup', DB3 / 'untitled folder 142/00_READ_FIRST.md', '6 September 2026'),
    ('Chat 8 scalar closure', DB3 / 'untitled folder 143/00_READ_FIRST.md', '6 September 2026'),
    ('Chat 8 information and quantum (dated folder name)', DB3 / 'untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906.zip', None),
    ('Chat 7 transcript', DB3 / 'untitled folder 142 new chat/HDBLAST_CHAT7_COMPLETE_TRANSCRIPT_AND_PROGRESS_20260906.md', '6 September 2026'),
    ('APS author-review package (not submitted)', DB3 / 'untitled folder 145 APS JOURNAL/HDBLAST-v24-APS-review-package/README.txt', 'Prepared 9 September 2026. Not submitted.'),
    ('Chat 9 shell dynamics', DB3 / 'untitled folder 146/00_READ_FIRST.md', '16 September 2026'),
    ('Chat 10 stable shell', DB3 / 'untitled folder 147/00_READ_FIRST.md', '18 September 2026'),
    ('Chat 11 symbolic verification', DB3 / 'untitled folder 148/00_READ_FIRST.md', '18 September 2026'),
    ('Chat 12 all sectors', DB3 / 'untitled folder 149/00_READ_FIRST.md', '19 September 2026'),
    ('Chat 13 nonlinear roll-off', DB3 / 'untitled folder 150/00_READ_FIRST.md', '21 September 2026'),
    ('Chat 14 registered-detuning roll-off', DB3 / 'untitled folder 151/00_READ_FIRST.md', '22 September 2026'),
    ('Registered-shell controls', DB3 / 'untitled folder 152/00_READ_FIRST.md', '22 September 2026'),
    ('Scalar-profile branch and matter', LATEST / '00_READ_FIRST.md', '22 September 2026'),
]
anchor_rows = []
for label, path, needle in anchors:
    ok = path.exists() and (needle is None or needle in path.read_text(errors='replace'))
    anchor_rows.append(dict(label=label, path=str(path), needle=needle, found=ok))
check('all chronology anchors found in their source files', all(a['found'] for a in anchor_rows),
      n=len(anchor_rows), missing=[a['label'] for a in anchor_rows if not a['found']])
# the latest-work copy and folder 153 carry the same overview
same = sha256(LATEST / '00_READ_FIRST.md') == sha256(DB3 / 'untitled folder 153/00_READ_FIRST.md')
check('latest-work 00_READ_FIRST.md identical to D-Blast 3/untitled folder 153 copy', same)

# ---------------------------------------------------------------------------
# 7. Zenodo chronology of the HDBLAST concept record
# ---------------------------------------------------------------------------
audit = (DB3 / 'VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md').read_text()
sect = audit.split('### HDBLAST / D-Blast primary', 1)[1].split('###', 1)[0]
row_re = re.compile(r'^\| (\d+) \| (historical|current) \| \[(\d+)\]\(https://zenodo.org/records/\d+\) \| (\d{4}-\d{2}-\d{2}) \| ([^|]*) \| (\w+) \| [^|]* \| ([^|]*) \| (.*) \|$', re.M)
chain = [dict(index=int(m[0]), status_at_audit=m[1], record=m[2], date=m[3], version_label=m[4].strip(),
              resource_type=m[5], amendment=m[6].strip(), title=m[7].strip()[:160]) for m in row_re.findall(sect)]
dates = [r['date'] for r in chain]
check('concept 17088132 chain parsed: 21 versions, dates non-decreasing',
      len(chain) == 21 and dates == sorted(dates), n=len(chain), first=chain[0]['date'] if chain else None,
      last=chain[-1]['date'] if chain else None)
later = [
    ('22285737', 'v22 = M489G-A, 2026-09-03', DB3 / 'untitled folder 137/M489GB_ZENODO_RECOMMENDATION_AND_PLAIN_LANGUAGE_SUMMARY.md', 'retain v22, DOI 10.5281/zenodo.22285737'),
    ('22287013', 'v23 = M489G-C, published 3 Sept 2026 per transcript', DB3 / 'untitled folder 142 new chat/HDBLAST_CHAT7_COMPLETE_TRANSCRIPT_AND_PROGRESS_20260906.md', 'Your HDBLAST v23 (M489G-C) is now published!'),
    ('22347452', 'v24 = conditional global adjoint, 2026-09-05', DB3 / 'untitled folder 142 new chat/HDBLAST_CHAT7_COMPLETE_TRANSCRIPT_AND_PROGRESS_20260906.md', 'version **2026.09.05-v24**'),
    ('22922928', 'scalar junction consistency / corrected static branch, 2026-09-23', DB3 / 'GPD-site-deploy/gpd5-articles.json', 'Newest HD-Blast record, 23 Sept 2026'),
]
later_rows = []
for rec, label, path, needle in later:
    txt = path.read_text(errors='replace')
    later_rows.append(dict(record=rec, label=label, path=str(path), record_found=rec in txt, phrase_found=needle in txt))
check('post-audit records v22, v23, v24 and 22922928 located with their identifying phrases',
      all(r['record_found'] and r['phrase_found'] for r in later_rows), rows=later_rows)
fake = '22922929'
control('fabricated record number 22922929 is absent from the same sources',
        not any(fake in Path(p).read_text(errors='replace') for _, _, p, _ in later))
result['zenodo_concept_17088132_chain_from_20260829_audit'] = chain
result['zenodo_post_audit_records'] = later_rows

# ---------------------------------------------------------------------------
# 8. Paths quoted in PROGRAM_MAP.md must exist
# ---------------------------------------------------------------------------
if MAP.exists():
    paths = sorted(set(re.findall(r'`(/home/user/[^`]+)`', MAP.read_text())))
    missing = [p for p in paths if not Path(p).exists()]
    check('every absolute path quoted in PROGRAM_MAP.md exists', not missing, n_paths=len(paths), missing=missing)
    control('a fabricated path is reported missing',
            not Path('/home/user/unified-theory-maldonado/new-files/D-Blast 3/untitled folder 999/00_READ_FIRST.md').exists())
    result['program_map_sha256'] = sha256(MAP)
else:
    check('PROGRAM_MAP.md present', False)

result['anchors'] = anchor_rows
result['checks'] = checks
result['n_checks'] = len(checks)
result['n_failed'] = sum(not c_['passed'] for c_ in checks)
result['all_passed'] = result['n_failed'] == 0
result['runtime_s'] = round(time.time() - t0, 3)
result['tolerances_note'] = ('Numerical tolerances are stated per check. Quoted source values are compared at the '
                             'precision printed in their source; expansions are tested by log-log residual slopes over '
                             'delta in {3e-4,1e-3,3e-3,1e-2}. Nothing here is an interval certificate.')
OUT.write_text(json.dumps(result, indent=2, default=float))
print(json.dumps(dict(n_checks=result['n_checks'], n_failed=result['n_failed'], all_passed=result['all_passed'],
                      runtime_s=result['runtime_s'])))
for c_ in checks:
    if not c_['passed']:
        print('FAILED:', c_['name'])
