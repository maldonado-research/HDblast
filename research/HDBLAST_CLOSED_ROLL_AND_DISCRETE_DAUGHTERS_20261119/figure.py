"""Render the field-crossing and static daughter-mode comparisons."""
from pathlib import Path
import argparse,json

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

ROOT=Path(__file__).resolve().parent


def make(output):
    p=json.loads((ROOT/'patches/PATCHES.json').read_text())
    s=json.loads((ROOT/'spectrum/SPECTRUM.json').read_text())
    fig,(left,right)=plt.subplots(1,2,figsize=(9.3,4.7))
    blue,orange='#165976','#a54926'
    for amp,color in [(0.001,blue),(-0.001,orange)]:
        r=next(r for r in p['runs'] if r['amplitude']==amp)
        rows=r['fractions_past_minimum']
        left.loglog([v['delta_x'] for v in rows],[100*v['fraction'] for v in rows],
                    marker='o',ms=4,color=color,label='Polar caps' if amp>0 else 'Equatorial belt')
    left.set_title('Where the field crosses first')
    left.set_xlabel('Proper time since first contact (model units)')
    left.set_ylabel('Volume past the minimum (%)')
    left.legend(loc='upper left',frameon=False)
    for role,color in [('unit',blue),('representative',orange)]:
        series=[r for r in s['cases'] if r['branch']=='phi_m=0.5' and r['radius_role']==role]
        series.sort(key=lambda r:r['cycles'])
        right.plot([r['cycles'] for r in series],[r['ratio_energy_exact_to_flat_continuum'] for r in series],
                   marker='o',ms=5,color=color,label=f"Radius {series[0]['radius_Htop_inverse']:g}")
    right.axhline(1,ls='--',lw=1,color='#717b7f')
    right.set_xticks([4,8,16])
    right.set_ylim(bottom=0)
    right.set_title('A separate static-sphere test')
    right.set_xlabel('Cycles of the imposed oscillation')
    right.set_ylabel('Energy relative to the flat continuum')
    right.legend(loc='upper left',frameon=False)
    for ax in [left,right]:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(color='#e0e6e8',lw=.7)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=.10,bottom=.15,right=.98,top=.86,wspace=.4)
    fig.savefig(output,dpi=165,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    make(p.parse_args().output)
