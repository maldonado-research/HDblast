"""Plot the closed-wall calculations and the equal-energy pressure comparison."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
INK = '#1e303d'
GRAY = '#516674'
BLUE = '#196575'
RED = '#a64a29'


def main():
    data = json.loads((ROOT/'junction/JUNCTION.json').read_bytes())
    cases = [row for key, row in data['cases'].items() if key.startswith('near_turnaround/')]
    cases.sort(key=lambda row: row['fraction_absorbed_of_initial_mass_parameter'])
    samples = [(100*row['fraction_absorbed_of_initial_mass_parameter'], row['end']['H']) for row in cases]
    plt.rcParams.update({'font.size': 11, 'text.color': INK, 'axes.labelcolor': INK,
                         'axes.titlecolor': INK, 'xtick.color': GRAY, 'ytick.color': GRAY,
                         'axes.edgecolor': '#b3bec5', 'axes.spines.right': False, 'axes.spines.top': False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5,5.5))
    fig.subplots_adjust(left=.093,right=.97,bottom=.22,top=.71,wspace=.37)
    fig.text(.093,.955,'How a pulse changes a model universe',size=20,weight='bold',va='top')
    fig.text(.093,.856,'Two gravitational tests in units of the anti-de Sitter length and time',color=GRAY)
    ax.axhline(0,color='#758a96',lw=1)
    ax.plot(*zip(*samples),color=BLUE,marker='o',lw=2,ms=7)
    ax.set_title('1   Add an incoming radiation pulse',loc='left',weight='bold',size=12.8)
    ax.set_xlabel('Initial bulk mass carried away on each side (%)')
    ax.set_ylabel('Expansion rate at the end of the pulse')
    ax.set_xticks([0,15,50,100]); ax.set_yticks([-.10,-.05,0,.025])
    ax.text(.97,.94,'Expanding',ha='right',va='top',transform=ax.transAxes,color=GRAY)
    ax.text(.05,.10,'Contracting',transform=ax.transAxes,color=GRAY)
    ax.grid(axis='y',alpha=.15)
    late = data['radiation_exit']['runs']['baseline']
    early = data['scalar_capture']['cases']['primary']
    before = [dict(t=0,a=2.),early['pulse_end'],early['plateau_entry'],early['bounce'],early['terminal']]
    trace = late['samples']
    begin, maximum = late['events']['rad90'][0],late['events']['maximum'][0]
    bx.axvspan(begin['t'],maximum['t'],color='#d8eae6',alpha=.85)
    bx.plot([s['t'] for s in before],[s['a'] for s in before],color=RED,ls=':',lw=1.6)
    bx.plot([s['t'] for s in trace],[s['a'] for s in trace],color=RED,lw=2)
    for row in [early['bounce'], begin, maximum]:
        bx.scatter(row['t'],row['a'],color=RED,s=25,zorder=3)
    bx.annotate('Bounce',(early['bounce']['t'],early['bounce']['a']),xytext=(10,-1),
                textcoords='offset points',color=GRAY,size=9)
    bx.annotate('Maximum radius',(maximum['t'],maximum['a']),xytext=(5,-10),
                textcoords='offset points',color=GRAY,size=8.7)
    bx.text(.47,.43,'Expanding with\n>90% radiation',ha='center',va='center',
            transform=bx.transAxes,color=BLUE,size=9.4)
    bx.set_title('2   Store the pulse energy in a field',loc='left',weight='bold',size=12.8)
    bx.set_xlabel('Time on the wall')
    bx.set_ylabel('Radius')
    bx.set_ylim(1.55,5.26)
    bx.grid(axis='y',alpha=.15)
    fig.text(.093,.066,'Left panel uses a slowly expanding wall. Right panel starts at rest. The dotted line joins early event samples.',size=9.1,color=GRAY)
    fig.text(.093,.033,'Shading marks expansion with over 90% of added wall energy in the radiation fluid. Both calculations use model units.',size=9.1,color=GRAY)
    fig.savefig(ROOT/'incident-wall.png',dpi=170,facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    main()
