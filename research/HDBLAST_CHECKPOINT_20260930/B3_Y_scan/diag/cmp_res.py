import sys, glob, numpy as np
sys.path.insert(0,'.')
import a1_analyze as A
base='/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main/'
tags=sys.argv[1].split(',')
taus=np.array([6,7,8,9,10,11,12,13,14,15,17,20,23,26,30,35,40,45])
D={}
for t in tags:
    s,ts=A.load(base+t+'_summary.json'); d=A.derived(s,ts); D[t]=d
print('tau   '+'  '.join('%26s'%t[-22:] for t in tags))
for tt in taus:
    row=[]
    for t in tags:
        d=D[t]; tau=d['H0tau']
        if tt>tau[-1]: row.append('%26s'%'-'); continue
        i=np.searchsorted(tau,tt)
        row.append('Wa4=%.4e r=%.4f H=%.4f lp=%.1f'%(d['Wa4'][i],d['r'][i],d['H_over_H0'][i],d['lapse'][i]))
    print('%5.1f '%tt+'  '.join(row))
