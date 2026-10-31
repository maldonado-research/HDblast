#!/usr/bin/env python3
"""Source-bound exact rational join to the selected exact-family slab.
No repository scientific code is imported. Frozen background receipt unchanged.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
CAPSULE=Path('research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031')
BACKGROUND_SHA='abd05a8a93fd2f543a112c344ca7f92cf95ce6d56e379cb80ba1e7d9dfb94882'

def need(ok,msg):
 if not ok:raise ValueError(msg)
def pin(path):
 raw=path.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def make_report(repo_root):
 repo_root=Path(repo_root);pins=json.loads((HERE/'FAMILY_DOMAIN_SOURCE_PINS.json').read_text())['files']
 for rel,want in pins.items():
  need(not Path(rel).is_absolute() and '..' not in Path(rel).parts,'relative source path required')
  need(pin(repo_root/rel)==want,'family domain source authentication: '+rel)
 manifest=json.loads((repo_root/CAPSULE/'MANIFEST.json').read_text());index={r['path']:r for r in manifest['files']}
 for rel,want in pins.items():
  sub=str(Path(rel).relative_to(CAPSULE))
  if sub=='MANIFEST.json':continue
  need(sub in index and index[sub]['bytes']==want['bytes'] and index[sub]['sha256']==want['sha256'],'retained exact-family manifest mismatch: '+rel)
 need(pin(HERE/'BACKGROUND_NORMAL.json')['sha256']==BACKGROUND_SHA,'frozen background receipt authentication')
 r=json.loads((HERE/'BACKGROUND_NORMAL.json').read_text());domain=json.loads((repo_root/CAPSULE/'domain/NORMAL.json').read_text())
 need(domain['status']=='PASS_PINNED_DOMAIN_AND_COMMON_SEED_ALGEBRA','source domain status')
 dd=domain['domain'];ell=F(r['ell']);need(F(dd['ell'])==ell and F(dd['slab_t_over_ell'])==F(3,2) and list(map(F,dd['annulus_z_over_ell']))==[F(1,2),F(5,2)],'selected exact-family slab and correction annulus')
 tauh=F(1,5500);H=F(8,5)*ell;C=1/(1-(15*H)**2/2);time_upper=tauh*C/129;Tphys=F(dd['slab_t_over_ell'])*ell
 tauhi=F(r['shoot']['tau_dagger'][1]);need(tauhi<131*5*ell/4<tauh,'target proper clock lies inside restricted family tube')
 need(time_upper<Tphys and time_upper/ell<F(10,7)<F(3,2),'restricted proper-time tube lies strictly in retained slab')
 need(F(r['actual_background_local']['rho'][0])>129,'retained rho lower bound')
 chi_min=F(r['compact_tube']['rapidity'][0]);need(abs(chi_min)<15*H,'auxiliary rapidity bound implies conservative cosh ceiling')
 label_hi=F(3,16)*ell;annulus_lo=F(1,2)*ell;need(label_hi<annulus_lo,'whole initial-label band inside exact initial zero collar')
 need(F(r['gaussian']['coordinate_determinant_dt_dz_over_dtau_dzi_lower'])>0 and F(r['gaussian']['radial_gram_lower'])>0,'inherit inverse-map and radial Gram margins on sub-tube')
 return {'schema':'proper-clock-selected-family-domain-join-v1','status':'PASS_SELECTED_EXACT_FAMILY_PROPER_CLOCK_DOMAIN_JOIN','frozen_background_sha256':BACKGROUND_SHA,'frozen_background_report_unchanged':True,
 'restricted_family_tube':{'proper_tau':['0',str(tauh)],'scaled_initial_z':['1/16','3/16'],'whole_sphere':True,'enclosing_scaled_t':['0','10/7'],'enclosing_scaled_z':r['compact_tube']['enclosing_scaled_r']},
 'exact_inequalities':{'target_tau_upper':str(tauhi),'target_tau_headroom':str(tauh-tauhi),'rough_target_tau_headroom':str(tauh-131*5*ell/4),'conservative_cosh_rapidity_upper':str(C),'coordinate_time_upper':str(time_upper),'scaled_coordinate_time_upper':str(time_upper/ell),'selected_coordinate_time_upper':str(Tphys),'physical_time_margin':str(Tphys-time_upper),'scaled_enclosing_box_time_margin':str(F(3,2)-F(10,7)),'initial_zero_collar_gap':str(annulus_lo-label_hi)},
 'scope_clarification':{'auxiliary_background_only':'The frozen tau<=1/5000, x<=8/5 tube is a valid background overcover. No selected exact-family existence beyond coordinate time3ell/2 is inferred.','family_persistence_application':'Apply the frozen report Gaussian stability argument only to the restricted tau<=1/5500 sub-tube, whose strict coordinate-time margin is proved here. The report premise requiring family definition near the chosen compact tube is satisfied here under the retained selected-family theorem.','parameter_radius':'Qualitative common small parameter neighborhood; no effective radius. The positive clock/inverse-map/Gram margins and strict selected-slab margin persist after parameter shrinkage.'},
 'initial_observer_data':{'selected_correction_annulus_over_ell':['1/2','5/2'],'initial_seed_support_over_ell':['1','2'],'exact_canonical_data_equal_background_near_all_launch_labels':True,'initial_lapse':'N_e=rho(1+e uY)=rho on the label neighborhood','initial_shift':'beta_e=0','initial_future_unit_normal':'rho^-1 partial_t, independent of parameter near all launch labels','initial_first_position_variation':'zero','initial_first_velocity_variation':'zero','initial_first_position_label_derivatives':'zero','initial_first_velocity_label_derivatives':'zero','first_position_system_W0':'0','reason':'Exact correction is supported in closure of ell/2<z<5ell/2; seed support is ell<=z<=2ell; harmonic initial completion is local and preserves the exact open zero collar.'},
 'source_sha256':{k:v['sha256'] for k,v in pins.items()}}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args();r=make_report(a.repo_root);raw=(json.dumps(r,indent=2,sort_keys=True)+'\n').encode()
 if a.output:need(not a.output.exists(),'refuse overwrite');a.output.write_bytes(raw)
 print(json.dumps({'status':r['status'],'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
