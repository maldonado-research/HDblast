#!/usr/bin/env python3
"""Validate the completed independent receipts and bind the accepted bytes."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]

def need(ok,label):
    if not ok:raise ValueError(label)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(name):return json.loads((ROOT/name).read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    need(not args.output.exists(),'OUTPUT_ALREADY_EXISTS')
    reports={};checks=[]
    for stem,status in (
        ('FRAME_COORDINATE_FINAL','PASS_INDEPENDENT_COMPLETE_FRAME_SOURCE_REVIEW'),
        ('SOURCE_INTERVALS_FINAL','PASS_INDEPENDENT_EXACT_SOURCE_INTERVAL_EVALUATION'),
        ('SHELL_ZEROS_INDEPENDENT','PASS_INDEPENDENT_FOUR_SIMPLE_SHELL_MEAN_ZEROS')):
        names=[f'review/{stem}_{mode}.json' for mode in ('NORMAL','OPTIMIZED')]
        need((ROOT/names[0]).read_bytes()==(ROOT/names[1]).read_bytes(),'NORMAL_OPTIMIZED_BYTE_EQUAL_'+stem)
        report=read(names[0]);need(report['status']==status,'REVIEW_STATUS_'+stem)
        for name,sha in report.get('sha256',report.get('input_sha256',{})).items():
            need(digest(ROOT/name)==sha,'CURRENT_REVIEW_DEPENDENCY_'+name)
        reports.update({name:digest(ROOT/name) for name in names})
        checks.append(stem+'_NORMAL_OPTIMIZED_IDENTICAL_AND_CURRENT')
    current=read('projection/GEOMETRIC_COEFFICIENT_MAPS.json')
    previous=read('review/ALL_POINT_PROJECTION_REPLAY.json')
    need({k:v for k,v in current.items() if k!='engine_pins'}==
         {k:v for k,v in previous.items() if k!='engine_pins'},'REVIEWED_UNIVERSAL_MAP_PAYLOAD_UNCHANGED')
    need((ROOT/'projection/GEOMETRIC_COEFFICIENT_MAPS.json').read_bytes()==
         (ROOT/'projection/GEOMETRIC_COEFFICIENT_MAPS_OPTIMIZED.json').read_bytes(),'FINAL_UNIVERSAL_MAP_NORMAL_OPTIMIZED_EQUAL')
    for name,sha in current['engine_pins'].items():need(digest(ROOT/'geometry'/name)==sha,'CURRENT_FORMAL_ENGINE_'+name)
    need({key:value['rank'] for key,value in current['certificate']['basis'].items()}==
         {'scalar':3,'vector':2,'tensor':7},'EXACT_COEFFICIENT_BASIS_RANKS')
    need(read('review/SOURCE_INTERVALS_FINAL_NORMAL.json')['maps_sha256']==digest(ROOT/'projection/GEOMETRIC_COEFFICIENT_MAPS.json'),'FINAL_INTERVAL_MAP_BINDING')
    checks.extend(['REVIEWER_REGENERATED_UNIVERSAL_MAP_PAYLOAD_IDENTICAL',
                   'FINAL_UNIVERSAL_MAP_NORMAL_OPTIMIZED_IDENTICAL','ALL_FINAL_ENGINE_PINS_CURRENT'])
    paths=[
      'geometry/indexed_geometry.py','geometry/geometric_sources.py','geometry/GEOMETRIC_SOURCE_DERIVATION.md',
      'projection/angular_generic.py','projection/angular_project.py','projection/formal_coefficients.py',
      'projection/source_coefficients.py','projection/check_complete_source_projection.py',
      'projection/check_projection.py','projection/check_source_maps.py',
      'projection/GEOMETRIC_COEFFICIENT_MAPS.json','projection/GEOMETRIC_COEFFICIENT_MAPS_OPTIMIZED.json',
      'projection/EXACT_GEOMETRIC_PROJECTION.md','numeric/actual_first_jets.py','numeric/exterior_error_caps.py',
      'numeric/interval_core.py','numeric/ORACLE_FREEZE.json','numeric/INPUT_PINS.json','numeric/PREDECESSOR_PINS.json',
      'shell-invariant/intrinsic_scalar.py','shell-invariant/check_shell_invariant.py','shell-invariant/SHELL_INTRINSIC_SCALAR.md',
      'GAUGE_COVECTOR_BOOKKEEPING.md','review/INDEPENDENT_DERIVATION.md','review/INDEPENDENT_GEOMETRIC_REVIEW.md',
      'review/independent_coordinate_oracle.py','review/check_producer_against_coordinate.py',
      'review/check_frame_against_coordinate.py','review/check_independent_source_intervals.py',
      'review/check_independent_shell_zeros.py','review/seal_review.py','review/ALL_POINT_PROJECTION_REPLAY.json',
      'review/COORDINATE_ORACLE_NORMAL.json','review/PRODUCER_CORE_NORMAL.json',*reports]
    files={name:{'sha256':digest(ROOT/name),'bytes':(ROOT/name).stat().st_size} for name in sorted(set(paths))}
    out={'status':'ACCEPTED_INTERNAL_INDEPENDENT_GEOMETRIC_SOURCE_REVIEW','blocking_findings':[],
         'review':'Internal mathematical review, not external peer review.',
         'accepted_scope':['Universal quadratic Ricci/Einstein/scalar/gauge identities and complete frame implementation',
             'Exact all-point selected-state angular maps with retained TT sectors and 23-slot forcing convention',
             'Independent continuous-cell actual-source interval evaluation, conditional on inherited actual C3 jet premise',
             'Exactly two simple intrinsic quadratic-mean zeros in each of two disjoint stated shell time brackets'],
         'checks':checks,'authoritative_replay_sha256':reports,'accepted_files':files,
         'analytic_proof_not_test_count':'Universal identities follow from inverse/connection difference algebra and complete formal polynomial reconstruction; coordinate families test implementation.',
         'map_rebind':'Final exact integer handling leaves every formal coefficient and all-point certificate unchanged. Independent prior export equals final export except the explicitly refreshed engine hashes; final normal/optimized maps and final coordinate replays are bound.',
         'separate_scope':'Independent complete-domain numeric jet/global-cap audit belongs to numeric-review; no second-response solve, nonlinear remainder, finite-amplitude force, or exclusion of mean zeros outside the two brackets.'}
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])

if __name__=='__main__':main()
