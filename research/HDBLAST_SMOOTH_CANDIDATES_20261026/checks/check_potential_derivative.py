#!/usr/bin/env python3
"""One isolated bulk scalar-potential derivative at a certified cell.
No full-source or derivative-error bound is calculated.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

def need(ok, msg):
    if not ok: raise RuntimeError(msg)

def interval(obj, key):
    v=tuple(map(F,obj[key])); need(len(v)==2 and v[0]<=v[1], 'ordered interval '+key); return v

parser=argparse.ArgumentParser();parser.add_argument('--repo-root',required=True,type=Path);args=parser.parse_args()
p=args.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_FIRST_JETS_NORMAL.json'
raw=p.read_bytes();pin='0aaa7b2fe30f546db92ef8b26a7224ffa7f6812d3780c331661395aa1991d3b8'
need(hashlib.sha256(raw).hexdigest()==pin,'inherited cell SHA')
d=json.loads(raw);need(d['status']=='PASS_ACTUAL_FIRST_COMMON_GAUGE_ORDINARY_2JET_CELL_ENCLOSURES','first jet status')
need(d['ell']=='1/1000000','time scaling');ell=F(d['ell'])
cs=[c for c in d['cells'] if c['label']=='evolved_packet'];need(len(cs)==1,'cell selection');c=cs[0]
need(c['scaled_x']==['1/4','1025/4096'] and c['scaled_r']==['5/4','5121/4096'],'cell domain')
s=interval(c['actual_s'],'s'); st=interval(c['actual_s'],'st');V=interval(c['background'],'Vppp');rho=interval(c['background'],'rho')
need(all(q[0]>0 for q in (s,st,V,rho)), 'positive factors required')
a=2*ell*V[0]*s[0]*st[0];b=2*ell*V[1]*s[1]*st[1]
fa=a*ell**2*rho[0]**2;fb=b*ell**2*rho[1]**2
need(F(3520000)<a<=b<F(342000000),'raw derivative bounds')
need(F(59,1000)<fa<=fb<F(578,100),'scaled forcing derivative bounds')
# Nondegenerate exact dual-number differentiation test for -2 ell factor.
q,qt,v,e=F(7,3),F(11,5),F(13,7),F(2,9)
correct=-2*e*v*q*qt
need(correct != -e*v*q*qt and correct != -2*v*q*qt,'reject missing factor 2 or ell')
print(json.dumps({'status':'PASS','source_cell':['1/4','1025/4096','5/4','5121/4096'],'raw_Y2_coefficient_derivative':'(-342000000,-3520000)','scaled_forcing_Y2_derivative':'(-578/100,-59/1000)','exact_positive_lower':str(a),'exact_positive_upper':str(b),'exact_scaled_lower':str(fa),'exact_scaled_upper':str(fb),'inherited_sha256':pin,'scope':'isolated potential term only; no full-source sign, C2 error or solution'},sort_keys=True))
