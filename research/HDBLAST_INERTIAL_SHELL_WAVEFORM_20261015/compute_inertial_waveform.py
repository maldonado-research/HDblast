"""Enclose the actual wider-state first variation using an inertial free comparator."""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


iv = module('inertial_comparator', HERE/'sol/inertial_comparator.py')
eb = module('inertial_error_bounds', HERE/'inertial_error_bounds.py')


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)):
        return [encode(x) for x in value]
    return value


def authenticated_json(relative, expected):
    path = HERE/relative
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == expected, 'mathematical input authentication failed: '+relative)
    return json.loads(raw), {'path':relative,'bytes':len(raw),'sha256':expected}


def calculate(output, cells=4096, tail_cells=256):
    require(type(cells) is int and 16 <= cells <= 16384 and cells%2 == 0, 'even support grid required')
    require(type(tail_cells) is int and 1 <= tail_cells <= cells and cells%tail_cells == 0, 'tail grid must divide support grid')
    output = Path(output)
    require(not output.exists(), 'fresh output directory required')
    output.mkdir(parents=True)
    bounds = eb.constants()
    ell = bounds['ell']
    require(ell == iv.ELL, 'same wider scale required')
    background, bg_pin = authenticated_json('archive/BACKGROUND_DERIVATIVE_BOUNDS.json',
        'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c')
    require(background['candidate_state_exact']['ell'] == str(ell), 'background collar must match wider state')
    gain, gain_pin = authenticated_json('inputs/previous/archive/PHYSICAL_GAIN_ENCLOSURE.json',
        '2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4')
    box = gain['correlated_retained_eta_box']
    G = tuple(map(F, box['G']['outward_rational_at_12_decimal_places']))
    K = tuple(map(F, box['K']['outward_rational_at_12_decimal_places']))
    require(G[0] < G[1] < 0 < K[0] < K[1], 'signed gain and positive prefactor required')
    mu, gamma, shell_pin = iv.load_shell(HERE)
    mass = iv.BumpMass(cells)
    gain_iv = iv.iv(*G)
    rows = []
    comparison_sum = [0,0]
    actual_sum = [0,0]
    tail_step = cells//tail_cells
    ranges = [(0,cells)] + [(cells+j,cells+j+1) for j in range(cells)] + [
        (2*cells+j*tail_step,2*cells+(j+1)*tail_step) for j in range(tail_cells)]

    def enclosure(left,right):
        comparison = iv.comparator_cell(left,right,mu,gamma,mass)
        qhi = F(iv.q_cell(left,right)[0][1], iv.S)
        error = eb.cell_error(left,right,qhi,bounds)
        actual = iv.add(comparison,iv.iv(-error,error))
        # Positive clock factor: cosh(t)^(-3/2) belongs to [1-t^2,1].
        factor = iv.iv(1-(ell*right)**2,1)
        physical = iv.mul(gain_iv,iv.mul(factor,actual))
        require(comparison[0]<=comparison[1] and actual[0]<=actual[1] and physical[0]<=physical[1],
                'ordered continuous enclosures required')
        return comparison,actual,physical,error

    for ln,rn in ranges:
        left,right = F(ln,cells),F(rn,cells)
        comp,actual,physical,error = enclosure(left,right)
        rows.append([ln,rn,*comp,iv.iv(error)[1],*actual,*physical])
        for total,value in [(comparison_sum,comp),(actual_sum,actual)]:
            squared = iv.square(value)
            total[0] += (rn-ln)*squared[0]
            total[1] += (rn-ln)*squared[1]
    denominator = cells*iv.S
    comp_integral = tuple(F(x,denominator) for x in comparison_sum)
    actual_integral = tuple(F(x,denominator) for x in actual_sum)
    C = K[0]*actual_integral[0]/ell,K[1]*actual_integral[1]/ell
    full_tail = iv.comparator_cell(F(2),F(3),mu,gamma,mass)
    actual_tail = iv.add(full_tail,iv.iv(-bounds['published_tail'],bounds['published_tail']))
    require(F(actual_tail[0],iv.S)>F('3.2e-8') and F(actual_tail[1],iv.S)<F('5.3e-8'),
            'strict signed actual tail gate')
    witnesses = {}
    for name,x in [('rising',F(5,4)),('central_left',F('1.49998')),('central_right',F('1.49999')),
                   ('falling',F(7,4)),('late_left',F('1.96')),('late_right',F('1.97'))]:
        comp,actual,physical,error = enclosure(x,x)
        witnesses[name] = dict(x=x,H0=[F(v,iv.S) for v in comp],H=[F(v,iv.S) for v in actual],
                               ell_S=[F(v,iv.S) for v in physical],error_upper=error)
    for name in ['rising','central_left','late_right']:
        require(witnesses[name]['H'][0]>0,'positive witness failed: '+name)
    for name in ['central_right','falling','late_left']:
        require(witnesses[name]['H'][1]<0,'negative witness failed: '+name)
    csv_path = output/'INERTIAL_WAVEFORM_CELLS.csv'
    with csv_path.open('w',newline='') as stream:
        writer = csv.writer(stream,lineterminator='\n')
        writer.writerow(['left_grid','right_grid','H0_lo','H0_hi','error_upper','H_lo','H_hi','ell_S_lo','ell_S_hi'])
        writer.writerows(rows)
    report = {'schema':'hdblast-inertial-shell-waveform-v1','status':'PASS',
        'state':{'delta':'1/1000','ell':ell,'T':3*ell,'Z':8*ell,'domain':[0,3],
                 'time_coordinate':'x=t/ell=tau/(rho_b ell)','initial_velocity':'exp(-z) Q(z/ell)'},
        'arithmetic':{'bits':iv.BITS,'value_denominator':iv.S,'decision_method':'Outward integer intervals and exact rational sums'},
        'error_certificate':bounds,'shell_parameter_intervals':{'mu':iv.record(mu),'gamma':iv.record(gamma)},
        'physical_gain_interval':G,'physical_norm_prefactor_interval':K,
        'bump_integral_interval':[F(v,iv.S) for v in mass.prefix[-1]],
        'tail_theorem':{'domain':[2,3],'free_H0_interval':[F(v,iv.S) for v in full_tail],
                        'actual_H_interval':[F(v,iv.S) for v in actual_tail],
                        'published_strict_H_interval':['3.2e-8','5.3e-8'],
                        'physical_coefficient_sign':'negative'},
        'zero_existence_brackets':[['1.49998','1.49999'],['1.96','1.97']],
        'zero_uniqueness_claimed':False,'sign_witnesses':witnesses,
        'quadrature':{'method':'Outward squares of whole-cell intervals, integrated over the full domain including the tail',
                       'domain':[0,3],'comparison_H0_integral_interval':comp_integral,
                       'actual_H_integral_interval':actual_integral,'C_interval':C},
        'cell_data':{'path':csv_path.name,'support_cells':cells,'tail_cells':tail_cells,'rows':len(rows),
                     'time_grid_denominator':cells,'value_denominator':iv.S,'bytes':csv_path.stat().st_size,
                     'sha256':sha256(csv_path.read_bytes()).hexdigest()},
        'mathematical_inputs':{'background':bg_pin,'gain':gain_pin,'shell':shell_pin},
        'source_pins':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in [
                      'compute_inertial_waveform.py','inertial_error_bounds.py','sol/inertial_comparator.py']},
        'scope':'Conditional actual coupled linear first-variation shell trace enclosure at the newly specified wider state. The comparator has exact constant-boundary inertia; potential and time-varying boundary residuals are bounded separately. No finite-amplitude nonlinear waveform, effective epsilon, observed time scale or cosmological-origin inference.'}
    (output/'INERTIAL_WAVEFORM_ENCLOSURE.json').write_text(json.dumps(encode(report),indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--cells',default=4096,type=int)
    parser.add_argument('--tail-cells',default=256,type=int)
    args=parser.parse_args()
    result=calculate(args.output_dir,args.cells,args.tail_cells)
    print(json.dumps({'status':result['status'],'rows':result['cell_data']['rows'],
                      'C_interval':encode(result['quadrature']['C_interval']),
                      'tail_H_interval':encode(result['tail_theorem']['actual_H_interval'])}))
