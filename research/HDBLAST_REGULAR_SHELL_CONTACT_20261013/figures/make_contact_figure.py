"""Geometric support schematic, not a computed perturbation solution."""
from pathlib import Path
import argparse,csv,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent);args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
# An illustrative admissible geometric scale. No numerical background bound is inferred.
tau_star=math.pi/8; a=tau_star/4; depth=tau_star/3
scaled=np.linspace(0,1,401);tau=scaled*tau_star
lower=np.maximum(0,a-tau)/tau_star;upper=(depth+tau)/tau_star
xx=np.linspace(0,1.45,401);end=(tau_star/2+np.arcsin(np.sin(tau_star/2)*np.sin(xx*tau_star)))/tau_star
with (args.output/'contact-support.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['tau_over_tau_star','minimum_possible_x_over_tau_star','maximum_possible_x_over_tau_star']);w.writerows(zip(scaled,lower,upper))
fig,ax=plt.subplots(figsize=(7.2,4.8),layout='constrained')
ax.fill_betweenx(scaled,lower,upper,color='#d7e7f2',label='Possible causal support')
ax.plot(lower,scaled,color='#226487',lw=1.6);ax.plot(upper,scaled,color='#226487',lw=1.6)
ax.plot([a/tau_star,depth/tau_star],[0,0],color='#132b3b',lw=6,clip_on=False,label='Compact initial annulus')
ax.axvline(0,color='#9d3c28',lw=2,label='Shell')
ax.axhline(a/tau_star,color='#9d3c28',ls=':',lw=1.3)
ax.plot(xx,end,color='#42484e',ls='--',lw=1.6,label='End of fully covered ambient slab')
ax.annotate('Arrival infimum',xy=(0,.25),xytext=(.55,.21),arrowprops={'arrowstyle':'->','color':'#9d3c28'},color='#793120',fontsize=10)
ax.text(.09,.7,'Nonzero shell trace occurs\narbitrarily soon after arrival.',fontsize=9,color='#793120')
ax.set(xlim=(-.08,1.45),ylim=(-.025,1.02),xlabel=r'Inward collar distance $x/\tau_*$',ylabel=r'Conformal time $\tau/\tau_*$',title='Controlled shell contact')
ax.legend(loc='upper right',fontsize=8,framealpha=.95);ax.grid(alpha=.16)
fig.savefig(args.output/'contact-support.png',dpi=170)
fig.savefig(args.output/'contact-support.pdf',metadata={'CreationDate':None,'ModDate':None})
plt.close(fig)
