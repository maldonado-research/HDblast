"""Select reviewed support and planning files for the gated final evidence builder."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
names = [
    'NEXT_SCIENTIFIC_TARGET.md',
    'publication/INDEPENDENT_FROZEN_SCIENCE_AND_REPLAY_REVIEW.json',
    'publication/INDEPENDENT_FRESH_REPLAY_AND_INITIAL_EVIDENCE_REVIEW.json',
    'publication/FINAL_EVIDENCE_HELPER_SOURCE_REVIEW.json',
    'publication/WIDER_STATE_PLANNING_REVIEW.json',
    'publication/LIVE_DISPLAY_REVIEW.json',
    'publication/INDEX_SOURCE_PREFLIGHT.json',
]
for name in ['freeze_source.py', 'assemble_source.py', 'authenticate_fresh_zip.py',
             'install_checkpoint.py', 'install_initial_evidence.py', 'integrate_summaries.py',
             'record_merges.py', 'stage_completed_evidence.py', 'validate_replay_integrity.py',
             'REPLAY_TOOL_VALIDATION.json', 'FINAL_EVIDENCE_TOOL_PREPARATION.json']:
    names.append('publication/' + name)
for path in sorted((ROOT / 'planning').rglob('*')):
    if path.is_file() and '__pycache__' not in path.parts:
        if path.suffix not in {'.py', '.json', '.md', '.txt'}:
            raise RuntimeError('Unexpected planning artifact')
        names.append(path.relative_to(ROOT).as_posix())
if len(names) != len(set(names)):
    raise RuntimeError('Duplicate evidence selection')
for name in names:
    path = ROOT / name
    if not path.is_file() or path.is_symlink():
        raise RuntimeError('Missing reviewed support artifact: ' + name)
args = [sys.executable, '-I', '-B', str(ROOT / 'publication/build_final_evidence.py')]
for name in names:
    args += ['--include', name]
subprocess.run(args, check=True)
