"""Generate continuous-cell waveform and norm enclosures using exact arithmetic."""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('waveform_intervals', HERE/'waveform_intervals.py')
iv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(iv)


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def rational(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:rational(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [rational(v) for v in value]
    return value


def calculate(output, cells=4096):
    if not isinstance(cells,int) or cells < 16 or cells > 16384 or cells % 2:
        raise ValueError('Even cell count between 16 and 16384 required')
    output=Path(output)
    require(not output.exists(),'Fresh output directory required')
    output.mkdir(parents=True)
    gain_path=HERE/'archive/PHYSICAL_GAIN_ENCLOSURE.json'
    gain=json.loads(gain_path.read_text())
    require(gain['status'].startswith('PASS'),'Actual-background gain audit must pass')
    box=gain['correlated_retained_eta_box']
    G=tuple(F(x) for x in box['G']['outward_rational_at_12_decimal_places'])
    K=tuple(F(x) for x in box['K']['outward_rational_at_12_decimal_places'])
    require(G[0]<G[1]<0 and 0<K[0]<K[1],'Correct signed physical gain and positive K required')
    constants=iv.sharp_waveform_error()
    ell,eta,tail,zeta=(constants[k] for k in ['ell','published_eta','published_tail','norm_error_upper'])
    require(constants['eta_upper']<eta,'Global uniform error constant is inadequate')
    require(constants['eta_tail_upper']<tail,'Tail error constant is inadequate')
    require(constants['eta_upper']**2+constants['eta_tail_upper']**2<zeta**2,'Observation norm error is inadequate')
    error,tail_iv=iv.interval(-eta,eta),iv.interval(-tail,tail)
    gain_iv=iv.interval(*G)
    rows=[]
    lower_sum=upper_sum=0
    ranges=[(0,cells)]+[(cells+k,cells+k+1) for k in range(cells)]+[(2*cells,3*cells)]
    for lnum,rnum in ranges:
        left,right=F(lnum,cells),F(rnum,cells)
        comparison=iv.comparison_cell(left,right,ell)
        if right<=1:
            actual=(0,0)
        elif left>=2:
            actual=tail_iv
        else:
            actual=iv.add(comparison,error)
            squared=iv.square(comparison)
            lower_sum+=squared[0];upper_sum+=squared[1]
        # cosh(t)^(-3/2) in [1-t^2,1] for 0<=t<=1/2.
        factor=iv.interval(1-(ell*right)**2,1)
        physical=iv.mul(gain_iv,iv.mul(factor,actual))
        require(comparison[0]<=comparison[1] and actual[0]<=actual[1] and physical[0]<=physical[1],
                'Reversed computed enclosure')
        rows.append([lnum,rnum,*comparison,*actual,*physical])
    integral=(F(lower_sum,cells*iv.SCALE),F(upper_sum,cells*iv.SCALE))
    low_norm=iv.sqrt_bounds(integral[0])[0]
    high_norm=iv.sqrt_bounds(integral[1])[1]
    true_norm=(max(F(0),low_norm-zeta),high_norm+zeta)
    true_integral=tuple(x*x for x in true_norm)
    C=(K[0]*true_integral[0]/ell,K[1]*true_integral[1]/ell)
    witnesses={}
    for label,x in [('rising',F(5,4)),('falling',F(7,4))]:
        h=iv.comparison_cell(x,x,ell)
        actual=iv.add(h,error)
        physical=iv.mul(gain_iv,iv.mul(iv.interval(1-(ell*x)**2,1),actual))
        witnesses[label]={'x':x,'H_interval':[F(v,iv.SCALE) for v in actual],
                          'ell_S_interval':[F(v,iv.SCALE) for v in physical]}
    require(witnesses['rising']['ell_S_interval'][1]<0,'Rising physical negative witness failed')
    require(witnesses['falling']['ell_S_interval'][0]>0,'Falling physical positive witness failed')
    path=output/'WAVEFORM_CELLS.csv'
    with path.open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n')
        writer.writerow(['left_grid','right_grid','comparison_lo','comparison_hi','H_lo','H_hi','ell_S_lo','ell_S_hi'])
        writer.writerows(rows)
    report={'schema':'hdblast-uniform-first-variation-waveform-v1','status':'PASS',
      'state':{'delta':'1/1000','ell':str(ell),'time_coordinate':'x=t/ell','domain':['0','3']},
      'arithmetic':{'kind':'Outward integer dyadic intervals and exact rational sums','bits':iv.BITS,
                    'transcendentals':'Reduced Taylor exponential with remainder and a separate large-argument tail'},
      'error_constants':constants,'physical_gain_interval':G,'physical_norm_prefactor_interval':K,
      'quadrature':{'method':'Continuous real cell ranges; all support-edge cells included',
                    'support_cells':cells,'comparison_integral_interval':integral,
                    'actual_H_integral_interval':true_integral,'C_interval':C},
      'sign_witnesses':witnesses,
      'cell_data':{'path':path.name,'rows':len(rows),'time_grid_denominator':cells,
                   'value_denominator':iv.SCALE,'bytes':path.stat().st_size,'sha256':sha256(path.read_bytes()).hexdigest()},
      'gain_input':{'path':'archive/PHYSICAL_GAIN_ENCLOSURE.json','sha256':sha256(gain_path.read_bytes()).hexdigest()},
      'source_pins':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['waveform_intervals.py','compute_waveform.py']},
      'scope':'Continuous pointwise enclosure of the linear first variation for the specified narrow state, conditional on retained analytic inputs. No numerical exact nonlinear waveform, effective epsilon interval, observational scale, external peer review or formal PDE certification.'}
    (output/'WAVEFORM_ENCLOSURE.json').write_text(json.dumps(rational(report),indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--cells',type=int,default=4096)
    args=parser.parse_args()
    result=calculate(args.output_dir,args.cells)
    print(json.dumps({'status':result['status'],'continuous_cells':result['cell_data']['rows'],
                      'quadrature':rational(result['quadrature']),'cells_sha256':result['cell_data']['sha256']}))
