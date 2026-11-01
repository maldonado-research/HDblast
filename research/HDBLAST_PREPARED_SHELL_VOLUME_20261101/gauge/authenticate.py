"""Portable standard-library authentication, before any science imports."""
import hashlib
import json
from pathlib import Path

COMPONENT=Path(__file__).resolve().parent.parent


class CheckFailure(ValueError):
    def __init__(self,witness):
        self.witness=witness
        super().__init__(witness)


def failure_payload(exc,mutation=None):
    codes={
        'GAUGE_RESIDUAL_NONZERO':'INITIAL_GAUGE_RESIDUAL',
        'same initial lapse-square value':'CONDITIONED_LAPSE_VALUE',
        'mixed velocity vector-wave join':'CONDITIONED_MIXED_VELOCITY',
        'pullback displacement':'CONDITIONED_PULLBACK',
        'canonical_momentum_factor: 3*Qm(y)':'PREPARATION_MOMENTUM',
        'output must be fresh':'OUTPUT_NOT_FRESH',
        'output must be outside source and component':'OUTPUT_IN_SOURCE',
    }
    witness=getattr(exc,'witness',codes.get(str(exc),'UNEXPECTED_FAILURE'))
    if witness=='INITIAL_GAUGE_RESIDUAL' and mutation in (
        'omit_radial_second_profile','omit_trace_inverse_cross','omit_metric_mixed_cross',
        'omit_lapse_square_value','halve_second_gauge_velocity'):
        witness+='_'+mutation.upper()
    return json.dumps({'status':'FAIL','witness':witness},sort_keys=True)+'\n'


def authenticate(repo_root):
    repo_root=Path(repo_root).resolve()
    essential_bytes=(COMPONENT/'gauge/ESSENTIAL_PINS.json').read_bytes()
    essential=json.loads(essential_bytes)
    checked=[]
    for row in essential['files']:
        relative=Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('invalid essential path')
        raw=(COMPONENT/relative).read_bytes()
        if len(raw)!=row['bytes'] or hashlib.sha256(raw).hexdigest()!=row['sha256']:
            raise CheckFailure('LOCAL_AUTHENTICATION')
        checked.append(row['path'])
    pinbytes=(COMPONENT/'gauge/SOURCE_PINS.json').read_bytes()
    pins=json.loads(pinbytes)
    for row in pins['sources']:
        relative=Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('invalid source path')
        path=(repo_root/relative).resolve()
        if not path.is_relative_to(repo_root):raise ValueError('source outside repository')
        raw=path.read_bytes()
        if len(raw)!=row['bytes'] or hashlib.sha256(raw).hexdigest()!=row['sha256']:
            raise CheckFailure('SOURCE_AUTHENTICATION')
    return {'essential_pin_sha256':hashlib.sha256(essential_bytes).hexdigest(),
            'source_pin_sha256':hashlib.sha256(pinbytes).hexdigest(),
            'local_files':checked,'upstream_keys':[r['key'] for r in pins['sources']]},pins


def output_path(value,repo_root):
    if value=='-':return None
    path=Path(value).expanduser().resolve()
    if path.exists():raise ValueError('output must be fresh')
    if path.is_relative_to(Path(repo_root).resolve()) or path.is_relative_to(COMPONENT):
        raise ValueError('output must be outside source and component')
    if not path.parent.is_dir():raise ValueError('output parent must exist')
    return path
