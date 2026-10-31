#!/usr/bin/env python3
"""Portable replay and receipt comparison; stdlib authentication precedes imports."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise ValueError(msg)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--receipt',type=Path,default=HERE/'BACKGROUND_NORMAL.json');a=p.parse_args()
 pins=json.loads((HERE/'LOCAL_PINS.json').read_text())
 for rel,pin in pins.items():
  raw=(HERE/rel).read_bytes();need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'local preimport authentication: '+rel)
 spec=importlib.util.spec_from_file_location('proper_clock_certificate',HERE/'background_certificate.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 expected=module.canonical_bytes(module.certify(a.repo_root));got=a.receipt.read_bytes()
 need(got==expected,'background receipt differs from exact recomputation')
 print(json.dumps({'status':'PASS_BACKGROUND_RECEIPT_REPLAY','sha256':hashlib.sha256(got).hexdigest(),'source_pins_before_repository_import':True},sort_keys=True))
if __name__=='__main__':main()
