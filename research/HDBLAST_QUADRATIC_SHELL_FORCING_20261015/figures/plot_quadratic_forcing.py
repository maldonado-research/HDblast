"""Plot continuous coefficient bounds; plotting floats are not proof arithmetic."""
import csv
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='hdblast-matplotlib-'))
os.environ.setdefault('XDG_CACHE_HOME', tempfile.mkdtemp(prefix='hdblast-font-cache-'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
raw = (ROOT/'computed/QUADRATIC_FORCING_CELLS.csv').read_bytes()
rows = list(csv.DictReader(raw.decode().splitlines()))
D = 2**112
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={'width_ratios': [1.6, 1]})
for ax, domain, scale, label in [(axes[0], (0, 3), 1, r'$\ell^2\langle Q_\phi\rangle$'),
                                (axes[1], (2, 3), 10**19, r'$10^{19}\ell^2\langle Q_\phi\rangle$')]:
    selected = [r for r in rows if int(r['left_grid']) >= domain[0]*4096
                and int(r['right_grid']) <= domain[1]*4096]
    x = [int(r['left_grid'])/4096 for r in selected]+[int(selected[-1]['right_grid'])/4096]
    lower = [float(Fraction(int(r['ell2_Q0_lo']), D)*scale) for r in selected]
    upper = [float(Fraction(int(r['ell2_Q0_hi']), D)*scale) for r in selected]
    lower.append(lower[-1]); upper.append(upper[-1])
    ax.fill_between(x, lower, upper, step='post', color='#156b89', alpha=.35, label='Whole-cell enclosure')
    ax.step(x, lower, where='post', color='#156b89', lw=.8)
    ax.step(x, upper, where='post', color='#156b89', lw=.8)
    ax.set(xlim=domain, xlabel=r'$x=\tau/(\rho_b\ell)$', ylabel=label)
    ax.grid(alpha=.18)
axes[0].set_title('Angular mean: complete time window')
axes[1].set_title('Continuous postpulse tail')
axes[1].axhline(-2.3, color='#a94f22', ls='--', lw=1, label='Strict theorem bounds')
axes[1].axhline(-.67, color='#a94f22', ls='--', lw=1)
axes[1].set_ylim(-2.5, -.45)
axes[1].legend(loc='best', frameon=False, fontsize=9)
fig.suptitle('Explicit quadratic material tension contribution — model units', fontsize=14)
fig.text(.5, .015, r'Outward-normal SUM: $Q_\phi=-4\psi_1^2$. Full second-order response remains to be evaluated.',
         ha='center', fontsize=10)
fig.tight_layout(rect=(0, .05, 1, .94))
target = ROOT/'computed/QUADRATIC_MATERIAL_FORCING.png'
fig.savefig(target, dpi=170, metadata={'Software': 'Matplotlib 3.10.8; deterministic source plot'})
pdf = ROOT/'computed/QUADRATIC_MATERIAL_FORCING.pdf'
fig.savefig(pdf, metadata={'Creator': 'Matplotlib 3.10.8', 'CreationDate': None, 'ModDate': None})
receipt = {'status': 'PASS', 'matplotlib_version': matplotlib.__version__,
           'input_csv_sha256': sha256(raw).hexdigest(), 'rows': len(rows),
           'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
           'artifacts': [{'path': p.relative_to(ROOT).as_posix(), 'bytes': p.stat().st_size,
                          'sha256': sha256(p.read_bytes()).hexdigest()} for p in [target, pdf]],
           'scope': 'Rendering of all continuous cell coefficient bounds. Float conversion is presentation only; exact rational CSV and proofs are authoritative.'}
(ROOT/'figures/FIGURE_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
