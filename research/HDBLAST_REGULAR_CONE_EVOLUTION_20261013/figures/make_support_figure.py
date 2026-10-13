#!/usr/bin/env python3
"""Geometric finite-propagation schematic; no numerical physical solution."""
from pathlib import Path
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
T=np.linspace(-1,3,401)
outer=1.5+np.abs(T-1)
inner=np.maximum(0,.5-np.abs(T-1))
shell=np.sqrt(4.5**2+T**2)
with (args.output/'support-geometry.csv').open('w',newline='') as f:
    out=csv.writer(f)
    out.writerow(['T_over_t0','support_inner_bound','support_outer_bound','shell_R_for_rb_over_t0_4p5'])
    out.writerows(zip(T,inner,outer,shell))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(7.2,5.4),layout='constrained')
ax.fill_betweenx(T,inner,outer,color='#bddbe8',alpha=.85,label='Possible causal support')
ax.plot(outer,T,color='#3a7392',lw=1.2)
ax.plot(shell,T,color='#804629',lw=2,label='Timelike shell (schematic)')
ax.plot(np.abs(T),T,color='#353c46',lw=1.5,ls='--',label='Characteristic cone')
ax.plot([.5,1.5],[1,1],color='#164766',lw=5,solid_capstyle='round',label='Compact initial-data annulus')
ax.scatter([1],[1],color='#962b36',s=45,zorder=5)
ax.annotate('Exact flux jets\nat $T=R=t_0$',(1,1),xytext=(2.35,.5),arrowprops={'arrowstyle':'->','color':'#962b36'},color='#76232c')
ax.text(3.4,2.05,'Fields vanish\nnear the shell',ha='center',va='center',color='#704427')
ax.scatter([0],[0],color='#353c46',s=20,zorder=5)
ax.set(xlabel=r'$R/t_0$',ylabel=r'$T/t_0$',xlim=(0,5.6),ylim=(-1,3))
ax.set_title('A local coupled solution before shell contact',loc='left',fontweight='bold')
ax.legend(loc='lower right',fontsize=8,frameon=True,framealpha=.96)
ax.grid(alpha=.15)
fig.savefig(args.output/'support-geometry.png',dpi=180)
fig.savefig(args.output/'support-geometry.pdf',metadata={'Title':'Finite propagation support schematic','Author':'Ricardo Maldonado with OpenAI Codex assistance'})
plt.close(fig)
