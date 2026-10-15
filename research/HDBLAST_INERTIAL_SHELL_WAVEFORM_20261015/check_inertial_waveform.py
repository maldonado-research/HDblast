"""Audit the wider-state continuous certificate and its mathematical conventions.

The separate referee/check_independent_cells.py reconstructs every CSV cell.
This fast checker authenticates those data, verifies the analytic error gates,
signed witnesses, physical normalization, and complete-cell squared norm.
"""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def run(data_dir):
    checks = []

    def check(label, condition):
        if not condition:
            raise RuntimeError(label)
        checks.append(label)

    data_dir = Path(data_dir)
    raw = (data_dir/'INERTIAL_WAVEFORM_ENCLOSURE.json').read_bytes()
    d = json.loads(raw)
    check('approved continuous certificate schema', d['schema']=='hdblast-inertial-shell-waveform-v1' and d['status']=='PASS')
    spec = importlib.util.spec_from_file_location('audited_error', ROOT/'inertial_error_bounds.py')
    eb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(eb)
    b = eb.constants()
    check('all exact residual constants reproduced', all(F(d['error_certificate'][k])==v for k,v in b.items()))
    ell = F(d['state']['ell'])
    check('new wider state and observation window', ell==F(1,10**6) and F(d['state']['T'])==3*ell and F(d['state']['Z'])==8*ell)
    check('specified nonzero incoming velocity retained', d['state']['initial_velocity']=='exp(-z) Q(z/ell)')
    check('positive majorant denominator', b['r']<1 and b['s']<1 and b['determinant']>0)
    check('global actual comparison error gate', b['global_error']<F('1.1e-7'))
    check('postpulse actual comparison error gate', b['tail_error']<F('1.7e-10'))
    check('full support plus tail norm error gate', b['published_global']**2+b['published_tail']**2<b['published_norm_error']**2)
    check('causal localization strict geometry', ell/2<ell and 3*ell+F(5,4)*F(7,2)*ell<8*ell and 3*ell<F(7,2)*ell)
    check('source pins authenticated', all(sha256((ROOT/name).read_bytes()).hexdigest()==pin for name,pin in d['source_pins'].items()))
    check('all mathematical inputs authenticated', all(sha256((ROOT/p['path']).read_bytes()).hexdigest()==p['sha256'] and (ROOT/p['path']).stat().st_size==p['bytes'] for p in d['mathematical_inputs'].values()))
    N, S = d['cell_data']['time_grid_denominator'], d['arithmetic']['value_denominator']
    check('approved outward integer scale', d['arithmetic']['bits']==112 and S==2**112 and d['cell_data']['value_denominator']==S and N==4096)
    csv_raw = (data_dir/d['cell_data']['path']).read_bytes()
    check('complete CSV authenticated', len(csv_raw)==d['cell_data']['bytes'] and sha256(csv_raw).hexdigest()==d['cell_data']['sha256'])
    rows = list(csv.DictReader(csv_raw.decode().splitlines()))
    check('complete support and tail cell count', len(rows)==4353 and d['cell_data']['rows']==4353 and d['cell_data']['support_cells']==4096 and d['cell_data']['tail_cells']==256)
    ranges = [(0,N)]+[(N+j,N+j+1) for j in range(N)]+[(2*N+16*j,2*N+16*(j+1)) for j in range(256)]
    check('continuous grid covers full zero-to-three domain', all((int(r['left_grid']),int(r['right_grid']))==p for r,p in zip(rows,ranges)))
    check('prefront has exact causal zero', all(int(rows[0][k])==0 for k in ['H0_lo','H0_hi','error_upper','H_lo','H_hi','ell_S_lo','ell_S_hi']))
    check('all interval endpoints ordered', all(int(r[a])<=int(r[z]) for r in rows for a,z in [('H0_lo','H0_hi'),('H_lo','H_hi'),('ell_S_lo','ell_S_hi')]))
    check('error included in every actual interval', all(int(r['error_upper'])>=0 and int(r['H_lo'])==int(r['H0_lo'])-int(r['error_upper']) and int(r['H_hi'])==int(r['H0_hi'])+int(r['error_upper']) for r in rows))
    check('all postpulse cells have positive H and negative physical coefficient', all(int(r['H_lo'])>0 and int(r['ell_S_hi'])<0 for r in rows if int(r['left_grid'])>=2*N))
    sums = [0,0]
    for r in rows:
        lo, hi = int(r['H_lo']), int(r['H_hi'])
        width = int(r['right_grid'])-int(r['left_grid'])
        lower = 0 if lo<=0<=hi else min(lo*lo,hi*hi)//S
        upper = -((-max(lo*lo,hi*hi))//S)
        sums[0] += width*lower
        sums[1] += width*upper
    norm = tuple(F(v,N*S) for v in sums)
    check('whole-cell square integral includes entire tail', norm==tuple(map(F,d['quadrature']['actual_H_integral_interval'])) and d['quadrature']['domain']==[0,3])
    G = tuple(map(F,d['physical_gain_interval']))
    K = tuple(map(F,d['physical_norm_prefactor_interval']))
    gain = json.loads((ROOT/d['mathematical_inputs']['gain']['path']).read_text())['correlated_retained_eta_box']
    check('correlated physical gain authenticated', G==tuple(map(F,gain['G']['outward_rational_at_12_decimal_places'])) and G[0]<G[1]<0)
    check('correlated norm prefactor authenticated', K==tuple(map(F,gain['K']['outward_rational_at_12_decimal_places'])) and 0<K[0]<K[1])
    C = tuple(map(F,d['quadrature']['C_interval']))
    check('correct inverse-width norm Jacobian', C==(K[0]*norm[0]/ell,K[1]*norm[1]/ell))
    check('strict rounded wider-state response', F('6.06e11')<C[0]<C[1]<F('7.44e11'))
    h = tuple(map(F,d['tail_theorem']['actual_H_interval']))
    h0 = tuple(map(F,d['tail_theorem']['free_H0_interval']))
    tail = b['published_tail']
    check('tail theorem encloses exact comparator residual', h[0]<=h0[0]-tail and h[1]>=h0[1]+tail)
    check('strict nonzero signed tail theorem', d['tail_theorem']['domain']==[2,3] and F('3.2e-8')<h[0]<h[1]<F('5.3e-8'))
    check('negative gain reverses physical tail sign', d['tail_theorem']['physical_coefficient_sign']=='negative' and G[1]*h[0]<0)
    witnesses = d['sign_witnesses']
    check('all six witness intervals signed', all(F(witnesses[k]['H'][0])>0 for k in ['rising','central_left','late_right']) and all(F(witnesses[k]['H'][1])<0 for k in ['central_right','falling','late_left']))
    check('central zero has opposite signed endpoints', [F(witnesses[k]['x']) for k in ['central_left','central_right']]==[F('1.49998'),F('1.49999')])
    check('late zero has opposite signed endpoints', [F(witnesses[k]['x']) for k in ['late_left','late_right']]==[F('1.96'),F('1.97')])
    check('only existence claimed for two zero intervals', d['zero_existence_brackets']==[['1.49998','1.49999'],['1.96','1.97']] and d['zero_uniqueness_claimed'] is False)
    return {'status':'PASS_INERTIAL_SHELL_WAVEFORM','check_count':len(checks),'checks':checks,
            'csv_rows':len(rows),'C_interval':list(map(str,C)),
            'report_sha256':sha256(raw).hexdigest(),'csv_sha256':sha256(csv_raw).hexdigest(),
            'scope':'Finite continuous certificate checks; analytic hypotheses and nonlinear existence remain separately stated.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,default=ROOT/'computed')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        result=run(args.data_dir)
        encoded=json.dumps(result,indent=2)+'\n'
        if args.output:
            args.output.write_text(encoded)
        print(encoded,end='')
    except (RuntimeError,ValueError,OSError,KeyError) as exc:
        print('FAIL: '+str(exc),file=sys.stderr)
        raise SystemExit(1)
