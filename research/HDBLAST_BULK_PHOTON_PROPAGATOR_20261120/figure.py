"""Draw the fixed-wall comparison from the released numerical outputs."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator

ROOT=Path(__file__).resolve().parent
WALL='#166a78'
WARP='#ad4a20'
INK='#202e37'
MUTED='#54636d'


def read(path):
    return json.loads((ROOT/path).read_bytes())


def numeric_rows(rows):
    return sorted((float(key),row) for key,row in rows.items())


def set_ticks(ax,ticks,axis='x'):
    target=ax.xaxis if axis=='x' else ax.yaxis
    target.set_major_locator(FixedLocator(ticks))
    target.set_major_formatter(FormatStrFormatter('%g'))


def main():
    electro=read('euclidean/EUCLIDEAN.json')
    retarded=read('retarded/RETARDED.json')
    forces=read('forces/FORCES.json')
    plt.rcParams.update({'font.size':11,'axes.labelcolor':INK,'axes.titlecolor':INK,
                         'xtick.color':MUTED,'ytick.color':MUTED,'text.color':INK,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#aebcc4','font.family':'DejaVu Sans',
                         'legend.handlelength':1.5})
    fig,axs=plt.subplots(1,3,figsize=(14.2,5.25))
    fig.subplots_adjust(left=.066,right=.982,top=.73,bottom=.235,wspace=.42)
    fig.text(.066,.942,'How strongly can a fixed wall feel the extra dimension?',fontsize=19.5,
             weight='bold',ha='left',va='top')
    fig.text(.066,.853,'The same wall geometry, probed with two different electromagnetic couplings',
             fontsize=11.3,color=MUTED)
    models=[('leaves','branches',WALL,'o',r'Wall weight  $G=(1-\phi^2)^2$'),
            ('selected_weight_leaves','alternate_weight',WARP,'s',r'Null-wave uplift  $G=e^{2A}$')]
    for ax in axs:
        ax.tick_params(length=3)
        ax.grid(axis='y',which='major',lw=.65,alpha=.45,color='#b1bec5')
        ax.set_axisbelow(True)
    a,b,c=axs
    for group,band,color,marker,label in models:
        rows=numeric_rows(electro[group]['0.0']['momenta'])
        rows=[(q,row) for q,row in rows if .01<=q<=32]
        a.loglog([q for q,_ in rows],[100*row['extra_response_fraction'] for _,row in rows],
                 color=color,marker=marker,lw=1.8,ms=4.5,label=label)
        spectral=[row for row in retarded[band]['0']['samples'] if .01<=row['mass']<=2.]
        b.loglog([row['mass'] for row in spectral],[row['density'] for row in spectral],
                 color=color,marker=marker,lw=1.8,ms=4.5)
    a.set_title('Continuum share of response',fontsize=13,weight='bold',loc='left',pad=11)
    a.set_xlabel(r'Spacelike momentum  $Q L_0$')
    a.set_ylabel('Share of wall response (%)')
    set_ticks(a,[.01,.1,1,10])
    set_ticks(a,[.01,.1,1,10,100],'y')
    a.set_ylim(.0018,170)
    a.legend(loc='lower right',fontsize=8.2,frameon=False)
    b.set_title('Modes felt by a wall source',fontsize=13,weight='bold',loc='left',pad=11)
    b.set_xlabel(r'Continuum mass  $m L_0$')
    b.set_ylabel(r'Spectral density  $\rho(m)$')
    set_ticks(b,[.01,.1,1])
    b.set_yticks([1e-12,1e-9,1e-6,1e-3,1])
    b.yaxis.set_minor_locator(NullLocator())
    b.text(.016,.000030,r'Uplift: $\rho\propto m$',color=WARP,fontsize=9)
    b.text(.048,1.5e-11,r'Wall: $\rho\propto m^{6.2}$',color=WALL,fontsize=9)
    b.set_ylim(8e-15,3)
    styles=[('0.5','#164954','^','-'),('0',WALL,'o','--'),('-0.5','#508697','s',':')]
    for name,color,marker,style in styles:
        rows=numeric_rows(forces['leaves'][name]['force']['samples'])
        rows=[(radius,row) for radius,row in rows if radius<=8]
        c.loglog([radius for radius,_ in rows],
                 [100*row['excess_over_zero_mode'] for _,row in rows],ls=style,lw=1.8,
                 color=color,marker=marker,ms=4.5,label=rf'$\phi_b={float(name):g}$')
    c.axhline(1,color=MUTED,lw=.9,ls=(0,(4,3)),zorder=0)
    c.set_title('Force between wall charges',fontsize=13,weight='bold',loc='left',pad=11)
    c.set_xlabel(r'Charge separation  $R/L_0$')
    c.set_ylabel('Excess over the Coulomb force (%)')
    c.set_xlim(1.9,8.5)
    set_ticks(c,[2,3,4,6,8])
    set_ticks(c,[.001,.01,.1,1,10],'y')
    c.legend(loc='lower left',frameon=False,fontsize=8.7,ncols=1)
    fig.text(.066,.102,r'First two panels use the central cut $\phi_b=0$. Last panel uses the wall weight at three cuts. The dashed horizontal line marks 1%.',
             fontsize=9.3,color=INK)
    fig.text(.066,.052,r'$L_0$ is the model length unit. All curves are calculations on fixed tuned leaves with an even photon field and no additional brane kinetic term.',
             fontsize=9.1,color=MUTED)
    fig.savefig(ROOT/'bulk-photon.png',dpi=220,bbox_inches='tight',pad_inches=.16,
                metadata={'Title':'Photon response in the extra dimension on fixed tuned HDblast leaves'})
    plt.close(fig)


if __name__=='__main__':
    main()
