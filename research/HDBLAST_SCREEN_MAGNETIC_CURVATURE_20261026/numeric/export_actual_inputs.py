#!/usr/bin/env python3
"""One authenticated apex export of actual first-field physical two-jets.

No source campaign, finite differences, or new magnetic curvature formula.
The frozen oracle's exact rational outward arithmetic and every C3 error
remain in force. Reports must be written to fresh external scratch paths.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from output_safety import guard_output, source_snapshot, verify_source_snapshot, source_labels, require_hash
EXTERIOR = 'HDBLAST_EXTERIOR_TENSOR_JETS_20261026'
EXTERIOR_PIN = '0bc16f65598ae4977a7b8ab70a66c13ec85d844c871d6fee7739c4ec472d1c70'
ORACLE_PIN = 'd36f8aca0b257ed8bc9de6b2cbacbdcde24ebf3b447823e54d8af006638e08cb'
CAPSULES = {
 'exterior': (EXTERIOR, EXTERIOR_PIN),
 'current': ('HDBLAST_COUPLED_MEAN_PILOT_20261025', 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de'),
 'I': ('HDBLAST_INPUT_PRECISION_20261023', 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61'),
 'P': ('HDBLAST_COMMON_GAUGE_FORCING_20261016', '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6'),
 'N': ('HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44'),
 'G': ('HDBLAST_GEOMETRIC_BULK_SOURCES_20261016', '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),
 'C': ('HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'),
 'T': ('HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023', 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'),
 'F': ('HDBLAST_FULL_ANGULAR_RESPONSE_20261026', 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac'),
 'precision35': ('HDBLAST_SHELL_FEEDBACK_PRECISION_20261025', '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4'),
}
SUFFIXES = ('', 'y', 't', 'yy', 'yt', 'tt')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(repo):
    """Independent tree validation, before importing any project Python."""
    result = {}
    for key, (label, pin) in CAPSULES.items():
        root = repo/'research'/label
        need(root.is_dir() and not root.is_symlink(), 'ROOT_TYPE:'+key)
        manifest = root/'MANIFEST.json'
        need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST:'+key)
        seen = {'MANIFEST.json'}
        for row in json.loads(manifest.read_text())['files']:
            rel = Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'MEMBER_PATH:'+key)
            seen.add(str(rel))
            path = root/rel
            need(path.is_file() and not path.is_symlink(), 'MEMBER_TYPE:'+key+'/'+str(rel))
            need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'MEMBER_BYTES:'+key+'/'+str(rel))
        for path in root.rglob('*'):
            need(not path.is_symlink(), 'TREE_SYMLINK:'+key)
            if path.is_file():
                need(str(path.relative_to(root)) in seen, 'UNLISTED_FILE:'+key+'/'+str(path.relative_to(root)))
        result[key] = {'label':label, 'manifest_sha256':pin, 'files_including_manifest':len(seen)}
    return result


def encode(value):
    if hasattr(value, 'lo') and hasattr(value, 'hi'):
        return [str(value.lo), str(value.hi)]
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key:encode(item) for key,item in value.items()}
    if isinstance(value, (tuple,list)):
        return [encode(item) for item in value]
    return value


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo, output = args.repo_root.resolve(), args.output
    guard_output(output,source_file=__file__,protected_roots=[
        ('DEPENDENCY:'+key,repo/'research'/label) for key,(label,pin) in CAPSULES.items()])
    sources_before = source_snapshot([__file__,HERE/'output_safety.py'])
    own_pin = sources_before[str(Path(__file__).resolve())]
    before = authenticate(repo)
    root = repo/'research'/EXTERIOR
    oracle_path = root/'inputs/actual_third_jets.py'
    require_hash(oracle_path,ORACLE_PIN,'ORACLE_EXTERNAL_PIN')
    oracle = load('screen_magnetic_authenticated_actual_jets', oracle_path)
    ctx = oracle.Context(repo)
    # This is the sole actual apex call; all quantities below reuse its jets.
    amps, scalar, bg = oracle.actual_frame_jets(ctx, F(1,4), F(5,4))
    symbols, audits = {}, {}
    for name, jet in {**amps, 's':scalar}.items():
        order = 0 if name in 'LM' else 1
        audits[name] = {'reconstruction_order':order, 'max_physical_derivative':2,
                        'max_actual_input_order':order+2, 'certified_actual_input_order':3}
        need(order+2 <= 3, 'ACTUAL_C3_SUFFICIENCY:'+name)
        for suffix in SUFFIXES:
            symbols[name+suffix] = oracle.physical_expression(jet,suffix,bg['rho'],bg['a'],ctx.ell).val()
    background_symbols = {}
    for name, jet in bg.items():
        for suffix in SUFFIXES:
            key = name+suffix
            # ay and py already name derivatives; retain exact original names
            # and their further derivatives in an unambiguous nested table.
            value = oracle.physical_expression(jet,suffix,bg['rho'],bg['a'],ctx.ell).val()
            background_symbols.setdefault(name,{})[suffix] = value
            if key not in symbols:
                symbols[key] = value
    # Explicit aliases use ordinary physical time and y, not e0-frame time.
    symbols['Ht'] = bg['theta'].val()
    symbols['R'] = (bg['rho']*bg['c']).val()
    symbols['h'] = (bg['theta']/bg['rho']).val()
    symbols['inverse_rho'] = (1/bg['rho']).val()
    symbols['angular_scale_inverse'] = (1/(bg['rho']*bg['c'])).val()
    symbols['radial_conversion'] = (1/(ctx.ell*bg['rho'])).val()
    ew, eu = ctx.jets.exterior_caps()
    summary = json.loads((root/'RESULTS_SUMMARY.json').read_text())
    tt_report = json.loads((root/'compute/TT_TWOJET_NORMAL.json').read_text())
    need(sha(root/'compute/TT_TWOJET_NORMAL.json') == summary['source_report_sha256'], 'TWOJET_SUMMARY_REPORT_PIN')
    need(summary['apex'] == ['1/4','5/4'] and tt_report['domain']['apex'] == summary['apex'], 'TWOJET_APEX')
    need(summary['primary_acceptance'] is True and tt_report['status'] == 'CERTIFIED_FINITE_ACTUAL_TOTAL_TT_TWOJET', 'TWOJET_ACCEPTANCE')
    total_tt = {key:value['exact'] for key,value in summary['twojet'].items()}
    # Read the maximum retained campaign result; do not execute that campaign.
    retained = tt_report['campaign_runs'][-1]['total_apex_twojet']
    for key, row in total_tt.items():
        need(row == retained[key], 'TWOJET_EXACT_ENDPOINT_IDENTITY:'+key)
    for key, row in total_tt.items():
        symbols[key] = ctx.I(*row)
    symbols['ut'] = symbols['u_x']/ctx.ell
    ctx.finish()
    after = authenticate(repo)
    need(after == before and sha(oracle_path) == ORACLE_PIN and sha(__file__) == own_pin, 'FROZEN_CLOSURE_OR_EVALUATOR_CHANGED')
    verify_source_snapshot(sources_before)
    result = {
      'schema':'HDBLAST_SCREEN_MAGNETIC_ACTUAL_APEX_INPUTS_V1',
      'status':'CERTIFIED_CONDITIONAL_ACTUAL_FIRST_FIELD_PHYSICAL_TWOJETS',
      'apex':{'x':'1/4','r':'5/4','t':'1/4000000'},
      'actual_apex_calls':1, 'source_campaign_runs':0,
      'evaluator_sha256':own_pin, 'oracle_path':str(oracle_path.relative_to(repo)), 'oracle_sha256':ORACLE_PIN,
      'source_integrity_before_and_after':source_labels(sources_before),
      'authentication_before':before, 'authentication_after':after,
      'all_ten_complete_capsules_authenticated_before_project_import_and_after_run':True,
      'oracle_certificate':ctx.certificate,
      'derivative_order_audit':audits,
      'error_caps_retained':{'master_by_total_order':ew,'gauge_by_mixed_order':eu},
      'conventions':{'derivatives':'Suffix y,t,yy,yt,tt are ordinary physical coordinate derivatives, not covariant derivatives.',
        'conversion':'t=ell*x; partial_y=-(ell*rho)^(-1)partial_r; yy includes a/(ell*rho)partial_r.',
        'first_metric':'A,B,C,D,E,L,M are the registered frame amplitudes, with retained common-gauge reconstruction.',
        'scalar':'s is retained as an optional auxiliary input; no curvature formula is assumed.',
        'theta':'theta=Ht=tanh(t); h=Ht/rho.',
        'observer':'Conditional fixed slicing observer and spatial y normal are supplied by the eventual map, with no observer-invariance claim.'},
      'symbols':symbols, 'background_physical_jets':background_symbols,
      'total_tt_import':{'path':str((root/'compute/TT_TWOJET_NORMAL.json').relative_to(repo)),
        'sha256':sha(root/'compute/TT_TWOJET_NORMAL.json'),
        'summary_sha256':sha(root/'RESULTS_SUMMARY.json'), 'exact_scaled_twojet':total_tt,
        'retained_campaign_index':len(tt_report['campaign_runs'])-1,
        'no_campaign_replayed':True, 'u_x_exact_interval':total_tt['u_x'],
        'ordinary_physical_ut':'ut=u_x/ell; e0 u=u_x/(ell*rho)'},
      'arithmetic':{'frozen_input_ring':'exact Fraction endpoints rounded outward on grid 2^112',
                    'decisions':'exact rational endpoints only; display floats are not proof'},
      'limits':['Retained smooth actual solution, differentiated integral identities, and C3 certification remain conditional premises.',
                'A new complete quadratic magnetic Weyl map must be supplied and authenticated separately.',
                'This input export alone makes no curvature or observer-invariance claim.']}
    with output.open('x') as handle:
        json.dump(encode(result),handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'output':str(output),'sha256':sha(output),
        'display_only':{name:[float(symbols[name].lo),float(symbols[name].hi)] for name in ['rho','c','a','theta','A','B','C','D','E','L','M','u_x','ut']}}),flush=True)


if __name__ == '__main__':
    main()
