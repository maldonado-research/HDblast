#!/usr/bin/env python3
"""Visualize certified whole-cell enclosures; floats are display-only."""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import tempfile


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    if a.output_dir.exists():
        raise ValueError('new figure output directory required')
    root=Path(__file__).resolve().parent.parent
    pin=json.loads((root/'response/PINS.json').read_text())['initial']
    raw=(root/'response'/pin['path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=pin['sha256']:
        raise ValueError('initial profile input hash mismatch')
    rows=json.loads(raw)['profile_cells']
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'hdblast-numerical-mpl'))
    os.environ.setdefault('XDG_CACHE_HOME',str(Path(tempfile.gettempdir())/'hdblast-numerical-fontcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    x=[float(Fraction(v['r'][0])) for v in rows]+[2.0]
    fig,axs=plt.subplots(2,1,figsize=(7.2,6.2),sharex=True,constrained_layout=True)
    for ax,key,scale,label in [
        (axs[0],'metric_angular_mean',1e22,'Angular metric mean coefficient / $10^{22}$'),
        (axs[1],'scalar_mean',1e31,'Scalar mean coefficient / $10^{31}$')]:
        lower=[float(Fraction(v[key][0]))/scale for v in rows]
        upper=[float(Fraction(v[key][1]))/scale for v in rows]
        ax.fill_between(x,lower+[lower[-1]],upper+[upper[-1]],step='post',
                        color='#245882',alpha=.55,label='Whole-cell enclosure')
        ax.axhline(0,color='black',lw=.6)
        ax.set_ylabel(label)
        ax.grid(alpha=.2)
        ax.legend(loc='upper right',frameon=False)
    axs[1].set_xlabel(r'Initial radial coordinate $r=z/\ell$; $t=0$, $\ell=10^{-6}$')
    axs[1].set_xlim(1,2)
    fig.suptitle('Selected formal second initial coefficients',fontsize=13)
    a.output_dir.mkdir(parents=True)
    fig.savefig(a.output_dir/'INITIAL_PROFILE_ENCLOSURES.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(a.output_dir/'INITIAL_PROFILE_ENCLOSURES.png',dpi=180,
                metadata={'Description':'Actual whole-cell enclosures; second parameter derivatives, not force or evolved waveform.'})
    plt.close(fig)
    record={'input_sha256':pin['sha256'],'whole_cells':len(rows),
            'scope':'Display of certified cell enclosures. Float conversion and plotting are not acceptance arithmetic. Neither curve is a physical force or an evolved second response.',
            'files':[]}
    for path in sorted(a.output_dir.iterdir()):
        b=path.read_bytes();record['files'].append({'name':path.name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
    (a.output_dir/'FIGURE_RECORD.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    main()
