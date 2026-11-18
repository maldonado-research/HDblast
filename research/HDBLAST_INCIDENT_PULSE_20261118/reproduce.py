"""Replay the initial constraints, wall response and interface work from the public files."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parent


def close(expected, value):
    """Cross-platform rounding tolerance for reported numerical values."""
    if isinstance(expected,dict) and isinstance(value,dict):
        return expected.keys()==value.keys() and all(close(expected[k],value[k]) for k in expected)
    if isinstance(expected,list) and isinstance(value,list):
        return len(expected)==len(value) and all(close(a,b) for a,b in zip(expected,value))
    if type(expected) is float and isinstance(value,(float,int)):
        return math.isclose(expected,value,rel_tol=1e-10,abs_tol=1e-12)
    return expected==value


def checked(script, saved, *args):
    flags=['-O'] if sys.flags.optimize else []
    run=subprocess.check_output([sys.executable,'-I','-B',*flags,ROOT/script,*args])
    if not close(json.loads((ROOT/saved).read_text()),json.loads(run)):
        raise ValueError('Reported values differ: '+script)
    print(f'Checked {saved}')


def main():
    with tempfile.TemporaryDirectory() as tmp:
        data=Path(tmp)
        with zipfile.ZipFile(ROOT/'evolution/wall_histories.zip') as archive:
            if archive.testzip():
                raise ValueError('Invalid data archive')
            archive.extractall(data)
        for meta in sorted(data.rglob('metadata.json')):
            d=json.loads(meta.read_text())
            for name,check in d['cases'].items():
                rows=(meta.parent/(name+'.csv')).read_bytes()
                if hashlib.sha256(rows).hexdigest()!=check['sha256']:
                    raise ValueError('Corrupt wall history: '+name)
        summary=ROOT/'B1_summary.json'
        checked('theory/rank_witness.py','theory/WITNESS.json')
        checked('scalar/rank_certificate.py','scalar/RANK.json',summary)
        checked('scalar/fixed_grid.py','scalar/FIXED_GRID.json',summary)
        checked('controlled/replay_matching.py','controlled/MATCHING.json')
        checked('evolution/report_evolution.py','evolution/REPORT.json',data/'runs')
        checked('interface.py','INTERFACE.json',data/'runs','--mass',data/'mass_runs','--half',data/'half_runs')
        checked('original/compare_histories.py','original/COMPARISON.json',data/'runs')


if __name__=='__main__':
    main()
