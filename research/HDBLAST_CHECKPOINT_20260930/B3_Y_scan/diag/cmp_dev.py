import sys, numpy as np, json
sys.path.insert(0,'.')
import a1_analyze as A
def get(p):
    s,ts=A.load(p); return A.derived(s,ts)
a1='/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main/'
R={'A1_5e-4_dzc2e-3_L16':get(a1+'main_dstar_Y3_dc1e-2_dzf5e-4_summary.json'),'A1_1e-3_dzc4e-3_L16':get(a1+'main_dstar_Y3_dc1e-2_dzf1e-3_summary.json'),
   'dev_5e-4_dzc4e-3_L10':get(sys.argv[1])}
W0=8.3703289e7
out={}
for tt in [5,8,10,13.73,15.73,17.73,19.73,21.73,23.73,27.73]:
    row={}
    for k,d in R.items():
        tau=d['H0tau']; keep=np.concatenate(([True],np.diff(tau)>1e-12))
        row[k]=dict(Wa4_drift=float(np.interp(tt,tau[keep],d['Wa4'][keep])/W0-1), H=float(np.interp(tt,tau[keep],d['H_over_H0'][keep])), phi=float(np.interp(tt,tau[keep],d['phi_b'][keep])))
    out['%.2f'%tt]=row
    print(tt, {k:('%.5f %.6f'%(v['Wa4_drift'],v['H'])) for k,v in row.items()})
json.dump(out,open(sys.argv[2],'w'),indent=1)
