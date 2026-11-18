"""Read-only extraction of registered B1 cutoff states; requires numpy once.
The reproduction checker consumes the resulting CSV/JSON using stdlib only.
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np

SOURCE_PATH = 'research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source'
COLUMNS = ['T','H0tau','H_over_H0','Wy','sigma','R','X','rad','chi','vac','kin','fric']
GRIDS = [('1e-3',.001),('5e-4',.0005)]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def extract(SOURCE, OUT):
    states=[]; constraints=[]; runs={}
    for ps in ('0.5','0.9'):
      for G in ('100','1000'):
       for y in ('0.1','1'):
        for gs,gnum in GRIDS:
         cell=f'ps{ps}_G{G}_y{y}_bcut';tag=f'main_{cell}_dzf{gs}'
         sp=SOURCE/'runs/main'/f'{tag}_summary.json';npz=SOURCE/'runs/main'/f'{tag}_timeseries.npz'
         summary=json.loads(sp.read_text());packed=np.load(npz,allow_pickle=False)
         data={str(k):packed['rec'][:,i] for i,k in enumerate(packed['cols'])}
         t=data['H0tau']; finite=np.isfinite(data['Hmax_near']) & np.isfinite(data['Mmax_near'])
         bad=np.flatnonzero(finite & ((data['Hmax_near']>.05)|(data['Mmax_near']>.05)))
         before_bad=(np.arange(len(t)) < bad[0]) if len(bad) else np.ones(len(t),dtype=bool)
         measured=np.flatnonzero(finite & before_bad)
         if len(measured)==0 or not t[measured[-1]]>6: raise ValueError('no measured reliable anchor past6')
         if np.any(np.diff(t)<=0): raise ValueError('nonmonotone physical time')
         selected=np.flatnonzero((t>=3)&(t<=6)&before_bad)
         if len(selected)==0 or not np.array_equal(selected,np.arange(selected[0],selected[-1]+1)):raise ValueError('noncontiguous selected rows')
         if not selected[0]>0 or not selected[-1]+1<len(t):raise ValueError('missing strict outside brackets')
         for i in selected:
          states.append(dict(tag=tag,index=int(i),**{k:repr(float(data[k][i])) for k in COLUMNS}))
         # All measured constraint anchors in the admitted prefix, including after6.
         for i in measured:
          constraints.append(dict(tag=tag,index=int(i),H0tau=repr(float(t[i])),Hmax_near=repr(float(data['Hmax_near'][i])),Mmax_near=repr(float(data['Mmax_near'][i]))))
         m=summary['matter'];p=summary['params']
         if (m['phistar'],m['G'],m['y'],p['dzf'])!=(float(ps),float(G),float(y),gnum):raise ValueError('registered label mismatch')
         runs[tag]=dict(cell=cell,grid=gs,rho_b=summary['rho_b'],matter={k:m[k] for k in ('phistar','G','y','b')},
           stop_reason=summary['stop_reason'],run_end_H0tau=float(t[-1]),measured_prefix_end_H0tau=float(t[measured[-1]]),
           first_measured_constraint_failure_H0tau=float(t[bad[0]]) if len(bad) else None,
           selected_count=len(selected),selected_first_index=int(selected[0]),selected_last_index=int(selected[-1]),
           before_window=dict(index=int(selected[0]-1),H0tau=float(t[selected[0]-1])),after_window=dict(index=int(selected[-1]+1),H0tau=float(t[selected[-1]+1])),
           source_summary=str(sp.relative_to(SOURCE)),source_summary_sha256=sha(sp),source_timeseries=str(npz.relative_to(SOURCE)),source_timeseries_sha256=sha(npz))
    for filename,rows,fields in [('STATES.csv',states,['tag','index']+COLUMNS),('CONSTRAINTS.csv',constraints,['tag','index','H0tau','Hmax_near','Mmax_near'])]:
      with (OUT/filename).open('w',newline='') as f:
       w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
    metadata=dict(protocol='New descriptive fixed windows [3,4],[4,5],[5,6]; nearest retained state to5; fixed before result evaluation; no original plateau-rule change.',
      source_root=SOURCE_PATH,source_registration_sha256=sha(SOURCE/'REGISTRATION.md'),source_analyzer_sha256=sha(SOURCE/'b1_analyze.py'),source_solver_sha256=sha(SOURCE/'evolve_b1.py'),
      windows=[[3,4],[4,5],[5,6]],sample_target_H0tau=5,registered_cutoff_cells=8,grids=2,runs=runs)
    (OUT/'INPUTS.json').write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n')
    for n in ('STATES.csv','CONSTRAINTS.csv','INPUTS.json'):print(n,sha(OUT/n))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path,help='Directory of the original B1 simulation histories')
    parser.add_argument('output',type=Path,help='Directory for the derived tables')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    extract(args.source,args.output)
