#!/usr/bin/env python3
"""Conditional production-approximation screen; does not produce particles.

Read the frozen external Chat14 ZIP without importing its solver. Derivatives
come from proper-time local polynomial fits to phi, with stored velocity used
only as a consistency comparison. All exact physical exclusions remain open.
"""
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ARCHIVE = ROOT/'source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip'
PREFIX = 'HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/runs/'
TAGS = ['v2_t1e3_plus_c4', 'v2_t1e3_plus_c4_fine', 'v2_t1e3_plus_c4_finecoarse']
THRESHOLD = .1
TARGETS = [.1,.25,.5,.75,.9]
WINDOWS = [.05,.10]

def sha(blob): return hashlib.sha256(blob).hexdigest()

def first_cross(rec, phi):
    yy = rec[:,1]-phi
    ix = np.flatnonzero((yy[:-1]<=0)&(yy[1:]>0))
    if not len(ix): raise RuntimeError('No crossing')
    i=int(ix[0]); a=-yy[i]/(yy[i+1]-yy[i])
    return (1-a)*rec[i]+a*rec[i+1]

def fit(rec, tau, halfwidth):
    ok = np.abs(rec[:,9]-tau)<=halfwidth
    if ok.sum()<9: raise RuntimeError('Insufficient fit points')
    x=(rec[ok,9]-tau)/halfwidth
    cp=np.polynomial.polynomial.polyfit(x,rec[ok,1],4)
    ch=np.polynomial.polynomial.polyfit(x,rec[ok,2],4)
    v=cp[1]/halfwidth; vd=2*cp[2]/halfwidth**2
    hh=ch[0]; hd=ch[1]/halfwidth
    stored=np.interp(tau,rec[:,9],rec[:,5])
    gh=hh*hh/(THRESHOLD**2*abs(v))
    gv=vd*vd/(THRESHOLD**2*abs(v)**3)
    ghd=abs(hd)/(THRESHOLD*abs(v))
    return {'halfwidth_H0_tau':halfwidth,'n_samples':int(ok.sum()),
        'v_dphi_dH0tau':float(v),'vd_d2phi_dH0tau2':float(vd),
        'H_over_H0':float(hh),'dHratio_dH0tau':float(hd),
        'stored_velocity':float(stored),
        'stored_vs_fitted_velocity_relative':float(abs(v-stored)/abs(v)),
        'Gmin_expansion':float(gh),'Gmin_scalar_acceleration':float(gv),
        'Gmin_curvature':float(ghd),'Gmin_all_three':float(max(gh,gv,ghd))}

def main():
    out={'status':'CONDITIONAL_APPROXIMATION_SCREEN','archive_sha256':sha(ARCHIVE.read_bytes()),
         'epsilon':THRESHOLD,'free_coupling_definition':'G=abs(gb/(kappa5*H0))=abs(gbar/H0)',
         'criteria':['abs(Hratio)/sqrt(G*abs(v)) <= epsilon',
                     'abs(vd)/(abs(v)*sqrt(G*abs(v))) <= epsilon',
                     'abs(dHratio_dH0tau)/(G*abs(v)) <= epsilon'],
         'scope':'Local fixed-trajectory approximation eligibility only; no quantum state, EFT cutoff, actual coupling, backreaction, thermalization or heating is solved',
         'runs':{},'comparisons':{},'checks':{}}
    with zipfile.ZipFile(ARCHIVE) as z:
        for tag in TAGS:
            meta=json.loads(z.read(PREFIX+tag+'_summary.json'))
            rec=np.load(io.BytesIO(z.read(PREFIX+tag+'_timeseries.npz')))['rec']
            causal=meta['L']*(1-meta['taper'])
            runs={}
            for phi in TARGETS:
                row=first_cross(rec,phi)
                fits=[fit(rec,float(row[9]),window) for window in WINDOWS]
                runs[str(phi)]={'coordinate_t':float(row[0]),'H0_tau':float(row[9]),
                    'before_taper_signal':bool(row[0]<causal),
                    'fits':fits,'window_relative_Gmin_change':float(abs(fits[0]['Gmin_all_three']/fits[1]['Gmin_all_three']-1))}
            out['runs'][tag]={'dz_fine':meta['dz_fine'],'dz_coarse':meta['dz_coarse'],
                              'causal_taper_threshold':causal,'crossings':runs}
        meta0=json.loads(z.read(PREFIX+TAGS[0]+'_summary.json'))
    for phi in TARGETS:
        vals=[out['runs'][tag]['crossings'][str(phi)]['fits'][0]['Gmin_all_three'] for tag in TAGS]
        out['comparisons'][str(phi)]={'Gmin_by_grid':dict(zip(TAGS,vals)),
              'relative_grid_spread':float((max(vals)-min(vals))/min(vals))}

    branchfile=ROOT/'static_branch/PLUS_BRANCH_RESULTS.json'
    branch=json.loads(branchfile.read_text())
    row=next(r for r in branch['rows'] if abs(r['delta']-.001)<1e-14)['refined']
    H2=row['H2']; rb0=meta0['rho_b']; alpha=H2*rb0**2
    ratio=(1-alpha)/alpha
    nmax=np.log(ratio)/4
    c=branch['c']
    W=lambda x:1-x+x**3/3
    sigma0=2*W(meta0['phi_b_static'])+.001*(1+c*meta0['phi_b_static'])
    sigmaf=2*W(row['phi_b'])+.001*(1+c*row['phi_b'])
    out['closed_4D_toy_budget']={'branch_source_sha256':sha(branchfile.read_bytes()),
        'initial_rho_b':rb0,'candidate_H2':H2,'alpha_residual_fraction':alpha,
        'maximum_radiation_to_vacuum_ratio':ratio,'maximum_radiation_dominated_efolds':float(nmax),
        'sigma_initial':sigma0,'sigma_candidate':sigmaf,'sigma_candidate_over_initial':sigmaf/sigma0,
        'premises':['Closed low-energy 4D energy budget with initial total 3*M4^2*H0^2',
                   'Same constant M4 at initial and final states',
                   'No subsequent injection or changing residual vacuum',
                   '100 percent of all nonvacuum energy converted instantly to radiation'],
        'qualification':'Not a bound on the registered 5D theory. The local low-energy identification kappa4^2=kappa5^2*sigma/6 changes about threefold across these configurations, so even the fixed-M4 premise requires a new effective-theory demonstration.'}
    out['checks']={
        'all_crossings_precede_taper':all(v['before_taper_signal'] for run in out['runs'].values() for v in run['crossings'].values()),
        'all_estimates_finite_positive':all(np.isfinite(f['Gmin_all_three']) and f['Gmin_all_three']>0 for run in out['runs'].values() for v in run['crossings'].values() for f in v['fits']),
        'closed_toy_fraction_between_zero_and_one':bool(0<alpha<1),
        'all_fit_windows_have_at_least_9_samples':all(f['n_samples']>=9 for run in out['runs'].values() for v in run['crossings'].values() for f in v['fits'])}
    if not all(out['checks'].values()): raise RuntimeError('Extraction check failed')
    (HERE/'CROSSING_ELIGIBILITY_AND_TOY_BUDGET.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':out['checks'],'comparisons':out['comparisons'],
          'toy_alpha':alpha,'toy_Nmax':float(nmax),'sigma_ratio':sigmaf/sigma0},indent=2))

if __name__=='__main__': main()
