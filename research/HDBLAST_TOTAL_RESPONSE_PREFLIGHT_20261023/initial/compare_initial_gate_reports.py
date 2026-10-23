#!/usr/bin/env python3
"""Compare numerical replay exactly, permitting only dependency-location metadata."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--replay',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rawa=a.baseline.read_bytes();rawb=a.replay.read_bytes();one=json.loads(rawa);two=json.loads(rawb)
if set(one)!=set(two):raise ValueError('report field inventory mismatch')
for key in one:
 if key!='manifest_closures' and one[key]!=two[key]:raise ValueError('numerical report mismatch: '+key)
if set(one['manifest_closures'])!=set(two['manifest_closures']):raise ValueError('closure role mismatch')
for role in one['manifest_closures']:
 omit={'root','external_pin_source','external_pin_source_sha256'}
 strip=lambda d:{k:v for k,v in d.items() if k not in omit}
 if strip(one['manifest_closures'][role])!=strip(two['manifest_closures'][role]):raise ValueError('closure identity mismatch: '+role)
r={'status':'PASS_EXACT_INITIAL_GATE_NUMERICAL_REPLAY_EQUALITY','baseline_sha256':hashlib.sha256(rawa).hexdigest(),'replay_sha256':hashlib.sha256(rawb).hexdigest(),'allowed_differences':['dependency root path','location/bytes of external-to-capsule pin source'],'all_numerical_fields_exactly_equal':True,'all_manifest_hashes_counts_and_bytes_equal':True}
a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
