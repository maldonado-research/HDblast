"""Plot the equal-energy timing comparison and a test of a fitted particle stream."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(root, output):
    source=json.loads((root/'source/SOURCE.json').read_bytes())
    counts=json.loads((root/'orbits/COARSE.json').read_bytes())['runs']
    trial=next(r for r in counts if r['history']=='dedicated_half' and r['pairs']==8192 and r['case']=='uniform_local_pressure_match')
    diff=100.*(np.array(trial['observed'])-np.array(trial['thermal']))
    fig, (a,b)=plt.subplots(1,2,figsize=(10,4.5),gridspec_kw={'width_ratios':[.85,1.]})
    vals=[source[name]['maximum_phi'] for name in ('simultaneous_time_control','delayed_matched_native')]
    bars=a.bar([0,1],vals,color=['#267783','#cb7331'],width=.55,zorder=2)
    a.set(ylim=(0,.84),yticks=np.arange(0,.81,.2),xticks=[0,1],xticklabels=['Both pulses\nearly','1/1,600 of energy\nheld for later'],ylabel='Largest scalar field value through wall time 8',title='When can particles form?')
    a.bar_label(bars,fmt='%.3f',padding=3)
    a.axhline(.53,color='#222222',ls='--',lw=1)
    a.text(.98,.54,'Decay threshold',ha='right',va='bottom',fontsize=9,color='#333333',transform=a.get_yaxis_transform())
    a.spines[['top','right']].set_visible(False); a.grid(axis='y',color='#e1e5e6'); a.set_axisbelow(True)
    limit=np.ceil(abs(diff).max())
    im=b.imshow(diff,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    b.set(xticks=range(4),xticklabels=['−1 to −½','−½ to 0','0 to ½','½ to 1'],
          yticks=range(5),yticklabels=['0–1','1–2','2–3','3–5','5 or more'],
          xlabel='Direction cosine along the stream',ylabel='Kinetic energy / thermal energy scale',title='Particle counts after fitting both pressures')
    for i,row in enumerate(diff):
        for j,val in enumerate(row):
            b.text(j,i,f'{val:+.1f}',ha='center',va='center',fontsize=9,color='white' if abs(val)>limit*.56 else '#212121')
    bar=fig.colorbar(im,ax=b,fraction=.044,pad=.025)
    bar.set_label('Simulated minus thermal count (percentage points)',fontsize=8)
    fig.subplots_adjust(bottom=.24,top=.85,left=.075,right=.96,wspace=.42)
    fig.text(.075,.055,'Left: equal total emitted energy.   Right: separate-interaction history, prescribed two-pulse pressure fit at 1% electrical work.',fontsize=9)
    fig.savefig(output,dpi=160,bbox_inches='tight',facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--output',type=Path,default=Path('polarization-results.png'))
    a=p.parse_args();render(a.study,a.output)
