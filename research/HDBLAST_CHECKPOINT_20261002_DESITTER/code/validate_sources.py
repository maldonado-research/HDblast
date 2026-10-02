#!/usr/bin/env python3
"""Fail-closed validator for the frozen numerical source benchmark."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('--registration',type=Path)
    a=ap.parse_args();root=Path(__file__).resolve().parents[1]
    registration=a.registration or root/'PREREGISTRATION.json'
    pin=hashlib.sha256(registration.read_bytes()).hexdigest()
    protocol=json.loads(registration.read_text())
    started=json.loads((a.run/'started.json').read_text())
    result=json.loads((a.run/'results.json').read_text())
    if (a.run/'failure.json').exists():raise ValueError('Run contains preserved failure artifact')
    if started['registration_sha256']!=pin or result['registration_sha256']!=pin:raise ValueError('Registration hash mismatch')
    if started['producer_sha256']!=protocol['local_sha256']['code/desitter_sources.py']:raise ValueError('Producer hash mismatch')
    expected={(x,y) for x in ['1','2'] for y in ['0.5','1']}
    observed={(p['x_over_r'],p['H2_over_r']) for p in result['points']}
    if observed!=expected or len(result['points'])!=4:raise ValueError('Missing/duplicate/unregistered points')
    if len(result['gates'])!=57 or len({g['name'] for g in result['gates']})!=57:raise ValueError('Missing or duplicate gates')
    for g in result['gates']:
        if not g['pass'] or Decimal(g['absolute_residual'])>Decimal(g['threshold']):raise ValueError('Failed gate: '+g['name'])
    if len(result['negative_controls'])!=3 or not all(g['detected'] for g in result['negative_controls']):raise ValueError('Missing or undetected negative control')
    if result['status']!='PASS':raise ValueError('Terminal status not PASS')
    print(json.dumps({'status':'PASS','grid_points':4,'gates':57,'negative_controls':3,'registration_sha256':pin},indent=2))

if __name__=='__main__':main()
