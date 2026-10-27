#!/usr/bin/env python3
"""Check inherited background enclosure and domain arithmetic; not a solver."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)

research = Path(__file__).resolve().parents[2]
path = research / 'HDBLAST_COMMON_GAUGE_FORCING_20261016' / 'bulk-producer' / 'inputs' / 'ACTUAL_POTENTIAL_JET_BOUNDS.json'
raw = path.read_bytes()
expected = 'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0'
require(hashlib.sha256(raw).hexdigest() == expected, 'inherited enclosure hash mismatch')
d = json.loads(raw)
require(d['status'] == 'PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ', 'inherited status mismatch')
ell, T, Z = (F(d['domain'][key]) for key in ('ell', 'T', 'Z'))
pmin,pmax = map(F,d['background_intervals']['p']['exact'])
rmin,rmax = map(F,d['background_intervals']['rho']['exact'])
require((ell,T,Z)==(F(1,1000000),F(3,1000000),F(1,125000)), 'domain mismatch')
require((pmin,pmax)==(F(-7,50000),F(-3,25000)) and pmin<=pmax<0, 'scalar gradient bound')
require((rmin,rmax)==(F(129),F(131)) and rmin>0, 'rho bound')
t=ell/4
z_target=5*ell/4
require(0<=t<=T and 0<=z_target and t+z_target<=Z, 'target/segment outside triangle')
require(d['domain']['triangle']=='0<=t<=T, z>=0, z+t<=Z', 'triangle convention mismatch')
print(json.dumps({'status':'PASS','inherited_sha256':expected,'conditional':'retained regular-branch theorem; checks inputs and containment, does not re-prove bounds','target_t':str(t),'target_z':str(z_target),'radial_segment':'0 <= z <= target_z at fixed target_t','p_bounds':[str(pmin),str(pmax)],'absolute_p_lower_bound':str(-pmax),'rho_bounds':[str(rmin),str(rmax)]},sort_keys=True))
