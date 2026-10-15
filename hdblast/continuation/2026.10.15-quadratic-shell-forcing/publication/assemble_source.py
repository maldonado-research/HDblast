"""Assemble only named reviewed scientific outputs; never copy scratch or raw papers."""
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT/'assembled-source'
FILES = [
    'README.md', 'REPRODUCE.md', 'requirements.txt', 'QUADRATIC_MATERIAL_TENSION_FORCING.md',
    'compute_quadratic_forcing.py', 'check_quadratic_forcing.py', 'verify_manifest.py',
    'computed/QUADRATIC_FORCING_ENCLOSURE.json', 'computed/QUADRATIC_FORCING_CELLS.csv',
    'computed/CHECK_QUADRATIC_FORCING.json', 'computed/QUADRATIC_MATERIAL_FORCING.png',
    'computed/QUADRATIC_MATERIAL_FORCING.pdf', 'figures/plot_quadratic_forcing.py', 'figures/FIGURE_RECEIPT.json',
    'archive/SOURCE_AUTHENTICATION.json', 'archive/REGISTERED_SHELL_COEFFICIENT_AUDIT.md',
    'archive/check_registered_shell_coefficients.py',
    'archive/REGISTERED_SHELL_COEFFICIENT_AUDIT_NORMAL.json',
    'archive/REGISTERED_SHELL_COEFFICIENT_AUDIT_OPTIMIZED.json',
    'astra/ASTRA_INPUT_PINS.json', 'astra/ASTRA_QUADRATIC_MATERIAL_BOUNDARY_FORCING.md',
    'astra/check_astra_material_forcing.py', 'astra/ASTRA_NORMAL.json', 'astra/ASTRA_OPTIMIZED.json',
    'sol/SOURCE_PINS.json', 'sol/QUADRATIC_MATERIAL_FORCING.md',
    'sol/check_quadratic_material_forcing.py', 'sol/QUADRATIC_MATERIAL_NORMAL.json',
    'sol/QUADRATIC_MATERIAL_OPTIMIZED.json',
    'referee/check_material_source.py', 'referee/MATERIAL_SOURCE_AUDIT_NORMAL.json',
    'referee/MATERIAL_SOURCE_AUDIT_OPTIMIZED.json', 'referee/MATERIAL_SCALAR_SOURCE_REVIEW.md',
    'referee/ROOT_CHECK_NORMAL.json', 'referee/ROOT_CHECK_OPTIMIZED.json',
    'referee/FINAL_QUADRATIC_SOURCE_REVIEW.json', 'referee/FINAL_QUADRATIC_SOURCE_REVIEW.md',
    'literature/LITERATURE_CHECK.md', 'literature/LITERATURE_LEDGER.json',
    'inputs/INPUT_PINS.json',
]
FILES += ['inputs/prior/'+n for n in [
    'INERTIAL_WAVEFORM_CELLS.csv', 'INERTIAL_WAVEFORM_ENCLOSURE.json', 'INERTIAL_SHELL_WAVEFORM.md',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json', 'PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
    'FINAL_INERTIAL_WAVEFORM_REVIEW.json', 'FINAL_INERTIAL_WAVEFORM_REVIEW.md']]
FILES += ['inputs/reference/'+n for n in ['SCALAR_SHELL_BOUNDARY_SYSTEM.md', 'MATERIAL_SECOND_ORDER_JUNCTIONS.md']]
if DEST.exists():
    raise RuntimeError('New assembly directory required')
records = []
for name in FILES:
    source = ROOT/name
    if not source.is_file() or source.is_symlink():
        raise RuntimeError('Missing reviewed scientific source: '+name)
    raw = source.read_bytes()
    target = DEST/name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(raw)
    records.append({'path': name, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})
for name in ['reproduce.py', 'REPLAY_ROSTER.json']:
    raw = (ROOT/'publication/replay-tools'/name).read_bytes()
    (DEST/name).write_bytes(raw)
    records.append({'path': name, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})
receipt = {'status': 'PASS', 'source_files': len(records), 'files': records,
           'scope': 'Explicit named reviewed source inventory. Raw responses, private references, subprocess copies and scratch outputs excluded.'}
(ROOT/'publication/SOURCE_ASSEMBLY_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'status': 'PASS', 'source_files': len(records)}))
