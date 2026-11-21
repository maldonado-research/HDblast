"""Plot five-dimensional wall transfer and two independent pressure diagnostics."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make(study, output):
    spectrum=json.loads((study/'geometry/SPECTRUM.json').read_bytes())
    frame=json.loads((study/'orbits/FRAME.json').read_bytes())
    fig, axes=plt.subplots(1,2,figsize=(10,4.7),gridspec_kw={'width_ratios':(1.15,1)})
    a,b=axes; blue='#23657a'; gold='#af5427'
    modes=np.arange(2,41)
    energy=[abs(spectrum['dedicated_half']['modes'][str(k)]['birth']['E']) for k in modes]
    a.plot(modes,energy,color=blue,linewidth=2.3)
    a.scatter([7],[energy[5]],s=35,color=gold,zorder=4)
    a.annotate('Mode 7',(7,energy[5]),xytext=(12,.93),arrowprops=dict(arrowstyle='-',color=gold))
    a.set(xlabel='Angular mode on the shell',ylabel='Electric field at first particle birth\n(per unit launch amplitude)',ylim=(-.02,1.09),xlim=(1,40),
          title='Which part of the pulse survives?')
    names=['Mean energy','Middle pressure share','Smallest pressure share']
    keys=['mean_energy_over_boosted_prediction','middle_principal_stress_fraction_over_prediction','least_principal_stress_fraction_over_prediction']
    for j,(channel,color,title) in enumerate((('dedicated_half',blue,'Separate interactions'),('shared_half',gold,'Shared interaction'))):
        row=frame[channel]['comparisons_by_resolution']['fine']
        ys=np.arange(3)+(j-.5)*.19
        b.scatter([row[k] for k in keys],ys,s=40,color=color,label=title,zorder=4)
    b.axvline(1.,color='#555555',linestyle=':',linewidth=1.)
    b.set(xscale='log',xlim=(8e-6,4.),ylim=(2.4,-.5),yticks=np.arange(3),yticklabels=names,
          xlabel='Ratio to the moving thermal prediction',title='What reaches the particles?')
    b.set_xticks([1e-5,1e-3,.1,1.],['0.00001','0.001','0.1','1'])
    fig.legend(*b.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.11),ncol=2,frameon=False,fontsize=9)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y' if ax is a else 'x',color='#dedede',lw=.6)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=.1,right=.97,bottom=.33,top=.87,wspace=.70)
    fig.text(.10,.045,'Pressure measurements use two sampled single-polarization waves in a local patch at 1% of total absolute electrical work.',fontsize=8)
    fig.savefig(output,dpi=165,bbox_inches='tight',facecolor='white')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=Path(__file__).parent)
    p.add_argument('--output',type=Path,default=Path('prebirth-response.png'))
    a=p.parse_args();make(a.study,a.output)
