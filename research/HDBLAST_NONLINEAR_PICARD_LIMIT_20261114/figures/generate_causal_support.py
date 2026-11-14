"""Draw the proved support envelope; this is not a wave-amplitude simulation.

Requires Matplotlib only for this optional static figure. The exact algebra
checks do not import this file. Outputs must be fresh and outside the capsule.
"""
import argparse
import io
from fractions import Fraction as F
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def output_path(value, capsule):
    p=value.absolute()
    if '..' in p.parts or p.is_relative_to(capsule) or p.exists() or p.is_symlink():
        raise ValueError('Fresh figure output outside capsule required')
    if not p.parent.is_dir() or any(q.is_symlink() for q in p.parents):
        raise ValueError('Regular existing output directory required')
    return p

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--preview-png',type=Path)
    args=parser.parse_args()
    capsule=Path(__file__).resolve().parent.parent
    out=output_path(args.output,capsule)
    preview=output_path(args.preview_png,capsule) if args.preview_png else None
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','svg.hashsalt':'hdblast-causal-support-20261114'})
    speed=F(20,19);final=F(13,8);arrival=F(19,40)
    times=sorted(set([final*F(i,240) for i in range(241)]+[arrival]))
    low=[max(F(0),F(1,2)-speed*s) for s in times]
    high=[F(5,2)+speed*s for s in times]
    fig,ax=plt.subplots(figsize=(9.2,5.5),layout='constrained')
    fig.patch.set_facecolor('white');ax.set_facecolor('#f2f4f7')
    x=[float(s) for s in times]
    ax.fill_between(x,[float(v) for v in low],[float(v) for v in high],color='#8db4d8',alpha=.7,label='Permitted support envelope')
    ax.plot(x,[float(v) for v in high],color='#164b78',lw=2)
    ax.plot(x,[float(v) for v in low],color='#164b78',lw=2)
    ax.axhline(0,color='#182838',lw=1.5)
    ax.axhline(6,color='#182838',lw=1.5,ls='--')
    ax.axvline(float(final),color='#687786',lw=1,ls=':')
    ax.plot([0,0],[.5,2.5],color='#164b78',lw=5,clip_on=False,label='Exact initial support annulus')
    ax.scatter([float(arrival)],[0],color='#164b78',s=32,zorder=4)
    ax.text(.10,4.6,'Perturbation is zero\noutside the envelope',color='#435564')
    ax.text(.83,1.9,'Support upper bound\nfor the reduced solution',ha='center',color='#123957')
    ax.text(.80,6.10,'Artificial outer wall',color='#182838')
    ax.text(.03,-.31,'Shell',color='#182838')
    arrow_x=1.75
    ax.annotate('',xy=(arrow_x,6),xytext=(arrow_x,float(F(80,19))),arrowprops={'arrowstyle':'<->','color':'#435564','lw':1.3},annotation_clip=False)
    ax.text(arrow_x+.035,5.12,'34/19\nmargin',va='center',fontsize=10,color='#435564')
    ax.set(xlim=(0,float(final)),ylim=(-.38,6.42),xlabel='Elapsed coordinate ratio s / ell',ylabel='Radial coordinate ratio z / ell')
    ax.set_xticks([0,float(arrival),float(final)],['0','19/40','13/8'])
    ax.set_yticks([0,.5,2.5,float(F(80,19)),6],['0','1/2','5/2','80/19','6'])
    ax.grid(alpha=.18)
    ax.text(.035,1.5,'Initial support bound',rotation=90,va='center',fontsize=9,color='#164b78')
    ax.set_title('Conditional causal support bound: radial speed at most 20/19',pad=12)
    fig.text(.5,-.025,'Shading bounds possible support. It gives no amplitude or guaranteed arrival.',ha='center',fontsize=10,color='#435564')
    buffer=io.StringIO()
    fig.savefig(buffer,format='svg',bbox_inches='tight',metadata={'Date':None,'Title':'Conditional reduced-system causal support envelope'})
    svg='\n'.join(line.rstrip() for line in buffer.getvalue().splitlines())+'\n'
    with out.open('x',encoding='utf-8') as stream:stream.write(svg)
    if preview:
        with preview.open('xb') as stream:fig.savefig(stream,format='png',dpi=170,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main()
