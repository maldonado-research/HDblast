"""Stdlib-only nonlinear frozen-state B1 reservoir bound.

This verifies pinned numerical snapshots and algebraic consequences, not PDE
error bounds, future evolution, or a replacement for the registered plateau rule.
"""
from pathlib import Path
from fractions import Fraction as F
import csv, hashlib, json, math

ROOT=Path(__file__).resolve().parent
PINS={
 'STATES.csv':'6863286bc42448a251339b1d736d925daeb4f98e0978507134ed389053c42d08',
 'CONSTRAINTS.csv':'f333e95f53e7190b2c3189fab0866b6602d747cca763d8fc81bc7d6cf3132a33',
 'INPUTS.json':'161b6a7f4eaab80f7198643bfddd2d19cfadb3993836a171aef76ec6cdc44940',
}

def require(ok,message):
    if not ok: raise ValueError(message)

def input_bytes(name):
    data=(ROOT/name).read_bytes()
    require(hashlib.sha256(data).hexdigest()==PINS[name],'source pin '+name)
    return data

def read_rows(name):
    lines=input_bytes(name).decode().splitlines()
    out=[]
    for row in csv.DictReader(lines):
        out.append({k:(v if k=='tag' else int(v) if k=='index' else float(v)) for k,v in row.items()})
    return out

def g(sigma,rho,rb): return rb*rb*(sigma*rho/18+rho*rho/36)

def root_required(target,sigma,rb):
    require(target>=0 and sigma>0 and rb>0,'root domain')
    if target==0:return 0.0
    q=36*target/(rb*rb)
    return q/(math.sqrt(sigma*sigma+q)+sigma)

def density_ratio(required,available):
    require(required>=0 and available>=0,'density ratio domain')
    if available==0:
        return (None,'infinite_requirement') if required>0 else (None,'zero_over_zero_no_energy_required')
    return required/available,'finite'

def derive(row,rb):
    sg,R,X,h,W=row['sigma'],row['R'],row['X'],row['H_over_H0'],row['Wy']
    require(all(math.isfinite(v) for k,v in row.items() if k not in ('tag','index')),'finite state')
    require(sg>0 and R>=0 and X>=0 and h>0,'admitted expanding positive-density state')
    rho=R+X; env=g(sg,rho,rb); H2=h*h
    require(math.isclose(g(sg,R,rb),row['rad'],rel_tol=2e-12,abs_tol=1e-25),'radiation normalization')
    require(math.isclose(env,row['rad']+row['chi'],rel_tol=2e-12,abs_tol=1e-25),'quadratic cross term and total-density identity')
    closure=sum(row[k] for k in ('vac','rad','chi','kin','fric','Wy'))-H2
    scale=sum(abs(row[k]) for k in ('vac','rad','chi','kin','fric','Wy'))+H2
    require(abs(closure)<=5e-13*max(scale,1e-30),'recorded nonlinear Friedmann identity')
    target_w=abs(W)/.1; target_share=.9*H2; target=max(target_w,target_share)
    req=root_required(target,sg,rb);factor,status=density_ratio(req,rho)
    require(math.isclose(g(sg,req,rb),target,rel_tol=2e-14,abs_tol=1e-25),'stable root substitution')
    require(env>0 and status=='finite','positive selected reservoir')
    return dict(tag=row['tag'],index=row['index'],H0tau=row['H0tau'],rho=rho,
      radiation_envelope=env,r_abs_min=abs(W)/env,Omega_max=env/H2,
      target=target,target_weyl=target_w,target_share=target_share,
      contribution_multiplier=target/env,required_density=req,density_multiplier=factor,
      share_density_multiplier=root_required(target_share,sg,rb)/rho,
      additional_density=max(0.0,req-rho),remaining_massive_fraction=X/rho,
      H2=H2,W=W,sigma=sg,closure_abs=abs(closure),
      target_weyl_over_share=target_w/target_share,
      actual_r_abs=abs(W)/row['rad'] if row['rad']>0 else None)

def extrema(rows):
    output=dict(count=len(rows),first_H0tau=min(r['H0tau'] for r in rows),last_H0tau=max(r['H0tau'] for r in rows))
    for key,fn in [('r_abs_min',min),('Omega_max',max),('density_multiplier',min),('contribution_multiplier',min),('share_density_multiplier',min)]:
        r=fn(rows,key=lambda x:x[key]);output[key]=dict(value=r[key],tag=r['tag'],index=r['index'],H0tau=r['H0tau'])
    return output

def compute():
    meta=json.loads(input_bytes('INPUTS.json'));states=read_rows('STATES.csv');constraints=read_rows('CONSTRAINTS.csv')
    expected={f'main_ps{ps}_G{G}_y{y}_bcut_dzf{dz}' for ps in ('0.5','0.9') for G in ('100','1000') for y in ('0.1','1') for dz in ('1e-3','5e-4')}
    require(set(meta['runs'])==expected,'all eight registered cutoff cells on both grids')
    require(meta['windows']==[[3,4],[4,5],[5,6]],'fixed windows')
    require(meta['sample_target_H0tau']==5,'fixed sample target')
    require(set(r['tag'] for r in states)==expected and set(r['tag'] for r in constraints)==expected,'raw coverage')
    cut={('0.5','100'):6.49e-6,('0.5','1000'):6.49e-9,('0.9','100'):1.39e-6,('0.9','1000'):1.39e-9}
    all_derived=[];out_runs={}
    for tag in sorted(expected):
        rm=meta['runs'][tag];raw=[r for r in states if r['tag']==tag];cs=[r for r in constraints if r['tag']==tag]
        require(len(raw)==rm['selected_count'],'selected row count')
        require([r['index'] for r in raw]==list(range(rm['selected_first_index'],rm['selected_last_index']+1)),'no omitted native records')
        require(all(3<=r['H0tau']<=6 for r in raw),'selected native time window')
        require(all(a['H0tau']<b['H0tau'] for a,b in zip(raw,raw[1:])),'state time monotonicity')
        require(rm['before_window']['H0tau']<3 and rm['after_window']['H0tau']>6,'outside window brackets')
        require(cs[0]['H0tau']<3 and cs[-1]['H0tau']>6,'measured constraint anchors bracket all windows')
        require(all(c['Hmax_near']<=.05 and c['Mmax_near']<=.05 for c in cs),'measured reliable prefix')
        require(rm['measured_prefix_end_H0tau']==cs[-1]['H0tau']>6,'strict measured prefix margin')
        if rm['first_measured_constraint_failure_H0tau'] is not None:
            require(rm['first_measured_constraint_failure_H0tau']>cs[-1]['H0tau'],'first failure after admitted anchors')
        mm=rm['matter'];ps=str(mm['phistar']);G=str(int(mm['G']))
        require(mm['b']==cut[ps,G],'registered b_cut')
        vals=[derive(r,rm['rho_b']) for r in raw];all_derived+=vals
        windows=[]
        for a,b in meta['windows']:
            vr=[v for v in vals if a<=v['H0tau']<=b]
            require(len(vr)>=2,'nonempty fixed window')
            right=min((c for c in cs if c['H0tau']>b),key=lambda c:c['H0tau'])
            checked=[c for c in cs if c['index']<=right['index']]
            window=extrema(vr);window.update(window=[a,b],right_constraint_anchor=right['H0tau'],
              measured_prefix_constraint_max=max(max(c['Hmax_near'],c['Mmax_near']) for c in checked))
            windows.append(window)
        sample=min(vals,key=lambda v:abs(v['H0tau']-5))
        out_runs[tag]=dict(cell=rm['cell'],grid=rm['grid'],selected_count=len(vals),windows=windows,sample_nearest_5=sample,
          measured_prefix_end_H0tau=cs[-1]['H0tau'],measured_prefix_constraint_max=max(max(c['Hmax_near'],c['Mmax_near']) for c in cs))
    pair=[]
    for cell in sorted({v['cell'] for v in out_runs.values()}):
        coarse=out_runs[f'main_{cell}_dzf1e-3'];fine=out_runs[f'main_{cell}_dzf5e-4']; ws=[]
        for cw,fw in zip(coarse['windows'],fine['windows']):
            ws.append(dict(window=fw['window'],relative_window_extremum_difference={k:abs(cw[k]['value']-fw[k]['value'])/abs(fw[k]['value']) for k in ('r_abs_min','Omega_max','density_multiplier','contribution_multiplier')}))
        a,b=coarse['sample_nearest_5'],fine['sample_nearest_5']
        pair.append(dict(cell=cell,windows=ws,nearest5_clock_difference=abs(a['H0tau']-b['H0tau']),
           nearest5_relative_difference={k:abs(a[k]-b[k])/abs(b[k]) for k in ('r_abs_min','Omega_max','density_multiplier')}))
    global_bounds=extrema(all_derived)
    global_bounds.update(min_H2=min(v['H2'] for v in all_derived),min_W=min(v['W'] for v in all_derived),
      min_sigma=min(v['sigma'] for v in all_derived),min_target_weyl_over_share=min(v['target_weyl_over_share'] for v in all_derived),
      max_closure_abs=max(v['closure_abs'] for v in all_derived),
      max_pair_window_relative_difference=max(x for p in pair for w in p['windows'] for x in w['relative_window_extremum_difference'].values()))
    require(global_bounds['density_multiplier']['value']>557.6,'density gap lower margin')
    require(global_bounds['contribution_multiplier']['value']>567,'contribution gap lower margin')
    require(global_bounds['Omega_max']['value']<.01854,'share envelope upper margin')
    # Exact toy checks distinguish the nonlinear contribution from density.
    sg,R,X,rb=map(F,(2,3,5,4))
    env=rb*rb*(sg*(R+X)/18+(R+X)**2/36)
    rad=rb*rb*(sg*R/18+R*R/36)
    chi=rb*rb*(sg*X/18+(X*X+2*R*X)/36)
    require(env==rad+chi,'exact quadratic envelope identity')
    require(env>rad+rb*rb*(sg*X/18+X*X/36),'omitting cross term is rejected')
    require(root_required(0,1,1)==0 and density_ratio(1,0)[1]=='infinite_requirement','zero density cases')
    return dict(status='PASS_NUMERICAL_FROZEN_STATE_ENVELOPE',input_pins=PINS,
      scope='Pinned stored states in new fixed windows, no continuum PDE error bound, no future plateau reclassification',
      state_count=len(states),run_count=len(out_runs),registered_cell_count=8,window_count_per_run=3,
      denominator='Original recorded h^2=(H/H0)^2; Weyl term also original. Envelope includes remaining X and nonlinear cross term.',
      global_bounds=global_bounds,runs=out_runs,paired_sensitivity=pair,
      conclusion='At all retained sampled states, full conversion of current chi energy still fails both same-state conservative diagnostics. This does not bound future production/scalar work or alter the registered INCONCLUSIVE verdict.')

if __name__=='__main__':print(json.dumps(compute(),indent=2,sort_keys=True,allow_nan=False))
