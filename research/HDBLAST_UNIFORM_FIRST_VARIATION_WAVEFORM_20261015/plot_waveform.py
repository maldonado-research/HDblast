"""Render already-certified cells. Floating arithmetic here is display only."""
import csv
from fractions import Fraction
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/hdblast-waveform-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'computed/WAVEFORM_ENCLOSURE.json').read_text())
scale=data['cell_data']['value_denominator']
grid=data['cell_data']['time_grid_denominator']
rows=list(csv.DictReader((ROOT/'computed/WAVEFORM_CELLS.csv').open()))
xs=[];lows=[];highs=[]
for row in rows:
    xs.extend([int(row['left_grid'])/grid,int(row['right_grid'])/grid])
    lows.extend([int(row['ell_S_lo'])/scale]*2)
    highs.extend([int(row['ell_S_hi'])/scale]*2)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.hashsalt':'hdblast-waveform-20261015'})
fig,(ax,tail)=plt.subplots(2,1,figsize=(10,7.4),gridspec_kw={'height_ratios':[3,1.35]})
fig.suptitle('Conditional first-variation shell waveform',fontsize=17,x=.51,y=.97)
ax.fill_between(xs,lows,highs,color='#33998d',alpha=.28,label='Certified continuous-cell envelope')
u=np.linspace(1,2,3001);q=np.zeros_like(u);qp=np.zeros_like(u)
inside=(u>1)&(u<2);v=1/((u[inside]-1)*(2-u[inside]))
q[inside]=np.exp(4-v);qp[inside]=(3-2*u[inside])*v*v*q[inside]
ell=float(Fraction(data['state']['ell']))
g=sum(float(Fraction(x)) for x in data['physical_gain_interval'])/2
illustration=g*(qp+ell*np.exp(-ell*u)*q)
ax.plot(u,illustration,color='#176a78',lw=1.5,ls='--',label='Illustrative comparison at midpoint gain')
ax.axhline(0,color='#444444',lw=.7);ax.axvline(1,color='#777777',lw=.7,ls=':');ax.axvline(2,color='#777777',lw=.7,ls=':')
ax.set(xlim=(0,3),ylabel=r'$\ell S(\ell x)$: scalar quadrupole coefficient',xlabel=r'$x=t/\ell$ (model coordinates)')
ax.legend(loc='upper right',fontsize=9,frameon=False);ax.grid(alpha=.15)
last=rows[-1]
tail.fill_between([2,3],[int(last['ell_S_lo'])/scale]*2,[int(last['ell_S_hi'])/scale]*2,color='#33998d',alpha=.28)
tail.axhline(0,color='#176a78',ls='--',lw=1)
tail.set(xlim=(2,3),ylabel=r'$\ell S(\ell x)$',xlabel=r'Post-pulse window: $x=t/\ell$')
tail.ticklabel_format(axis='y',style='sci',scilimits=(0,0));tail.grid(alpha=.15)
tail.set_title('Certified small tail band; tail shape is unresolved',fontsize=10)
fig.text(.09,.035,r'$\delta=10^{-3}$, $\ell=10^{-10}$. The gain interval dominates the visible pulse envelope.',fontsize=9)
fig.text(.09,.013,'Linear first variation under retained analytic hypotheses; no finite-amplitude nonlinear waveform or observed time scale.',fontsize=8.5)
fig.subplots_adjust(left=.11,right=.97,top=.91,bottom=.13,hspace=.60)
out=ROOT/'figures';out.mkdir(exist_ok=True)
fig.savefig(out/'SHELL_FIRST_VARIATION_WAVEFORM.png',dpi=200,metadata={'Software':'Matplotlib; certified source cells rendered separately'})
fig.savefig(out/'SHELL_FIRST_VARIATION_WAVEFORM.svg',metadata={'Date':'2026-10-15','Description':'Conditional linear first-variation envelope; display only'})
plt.close(fig)
print(json.dumps({'status':'RENDERED_FROM_CERTIFIED_CELLS','figures':2,'display_arithmetic_only':True}))
