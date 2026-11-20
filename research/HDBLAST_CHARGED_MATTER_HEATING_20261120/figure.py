"""Draw the charged-fluid transfer and the activation-time control from CHARGED.json."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def draw():
    data = json.loads((ROOT/'coupled/CHARGED.json').read_bytes())
    runs = data['runs']
    primary = runs['primary']
    curve = dict(zip(data['trace_columns'], np.array(primary['trace']).T, strict=True))
    show = (curve['t'] >= 4.4) & (curve['t'] <= primary['maximum_tau'])
    times = curve['t'][show]
    colors = ('#245b73','#a24725','#4e5960')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.labelcolor':colors[2],
                         'axes.edgecolor':'#b5bfc2','xtick.color':colors[2],'ytick.color':colors[2]})
    fig,(ax,bx) = plt.subplots(1,2,figsize=(11.3,5.3),gridspec_kw={'width_ratios':[1.2,1.]})
    fig.subplots_adjust(left=.08,right=.98,bottom=.22,top=.78,wspace=.31)
    fig.text(.08,.955,'When can the wave heat the gas?',fontsize=20,weight='bold',va='top')
    fig.text(.08,.862,'A closed universe with a pulse arriving from a fifth dimension',color=colors[2])
    ax.plot(times,100*curve['sheet'][show]/curve['e'][show],label='Electromagnetic field',
            color=colors[0],lw=2.4)
    ax.plot(times,100*curve['intrinsic_random'][show]/curve['e'][show],label='Particle internal energy',
            color=colors[1],lw=2.4)
    ax.axvline(data['parameters']['primary']['activation'],color=colors[2],ls=':',lw=1.2)
    ax.text(10.64,35,'Charges become\nmobile',fontsize=8,color=colors[2])
    ax.set_title('Energy during expansion',loc='left',weight='bold')
    ax.set_xlabel('Time on the wall (model units)')
    ax.set_ylabel('Share of added wall energy (%)')
    ax.legend(frameon=False,loc='upper left',fontsize=8.6,bbox_to_anchor=(.0,.57))
    names=['activation_4p7','activation_7','activation_9p5','coarse','activation_11','activation_11p5']
    x=[data['parameters'][n]['activation'] for n in names]
    y=[100*runs[n]['maximum_intrinsic_fraction'] for n in names]
    bx.plot(x,y,'o-',color=colors[1],lw=1.4,ms=5)
    bx.set_title('Same gas, different activation times',loc='left',weight='bold')
    bx.set_xlabel('When charges become mobile (model units)')
    bx.set_ylabel('Peak internal energy share (%)')
    bx.set_xticks([5,7,9,11])
    for a in (ax,bx):
        a.set_ylim(-3,103);a.grid(axis='y',alpha=.19)
    fig.text(.08,.080,'The right panel varies only the activation time on a common numerical grid. The gas starts cold and globally neutral.',fontsize=9,color=colors[2])
    fig.text(.08,.035,'Internal energy follows an assumed collision law. The expanding wall reaches a maximum radius and then contracts.',fontsize=9,color=colors[2])
    fig.savefig(ROOT/'charged-heating.png',dpi=175)
    plt.close(fig)


if __name__=='__main__':
    draw()
