#!/usr/bin/env python3
"""Pure exact replay for selected weighted initial Sobolev output bounds.

The closed manifest authenticates helper code. This evaluator authenticates
the seven retained input sources, then derives all new finite source jets,
polynomial coefficient envelopes and numerical joins. It writes no files.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from types import ModuleType


PINS = {
    'EXACT_INITIAL_THEOREM.md':'95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b',
    'REGISTERED_EXACT_FAMILY_BINDING.md':'69bfe846cc4fda8ec6a981969fbdd12fff7862b5ca7cf31a19d696b0a839165a',
    'CANONICAL_CONSTRAINTS.md':'4e25e993f9653b530e30e2a2d3fefcb0a006f04050534190b99a98126ba04f63',
    'EXPLICIT_COMPACT_SECOND_INITIAL_JET.md':'2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed',
    'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
    'PRIOR_SCALAR_OUTPUT.md':'67c3fe77b0fc6c2d0ce07e7466a46148e7983f93c4f14dae0344f00a2c6dc925',
    'PRIOR_ACTUAL_SOURCE.md':'e64d591585c0958dde28ad9428b8ca44fb61a9c441b33f28e7fc44693c007a4b',
}
UPSTREAM_REF='8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c'


def need(condition,message):
    if not condition:
        raise ValueError(message)


def load_json(raw):
    def unique(items):
        out={}
        for key,value in items:
            need(key not in out,'Duplicate JSON key: '+key)
            out[key]=value
        return out
    def nonfinite(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(raw.decode('utf-8'),object_pairs_hook=unique,
                      parse_float=str,parse_constant=nonfinite)


def helper(root,name):
    path=root/(name+'.py')
    need(path.is_file() and not path.is_symlink(),'Regular helper source required: '+name)
    # Compile directly so evaluate remains free of bytecode-cache writes even
    # when called by an API consumer that omitted the recommended -B flag.
    module=ModuleType('weighted_regularity_'+name)
    exec(compile(path.read_bytes(),str(path),'exec'),module.__dict__)
    return module


def evaluate(root):
    root=Path(root)
    names=[]
    def check(condition,name):
        need(condition,name)
        need(name not in names,'Duplicate mathematical check name: '+name)
        names.append(name)
    pins=load_json((root/'SOURCE_PINS.json').read_bytes())
    check(pins.get('schema')=='hdblast-weighted-regularity-pins-v1'
          and isinstance(pins.get('files'),list),'retained input pin inventory schema')
    seen=set()
    for row in pins['files']:
        local=row.get('local')
        check(isinstance(local,str) and local.startswith('inputs/')
              and local[7:] in PINS and local not in seen,
              'retained input local member '+str(local))
        seen.add(local)
        raw=(root/local).read_bytes()
        check(hashlib.sha256(raw).hexdigest()==PINS[local[7:]]==row.get('sha256')
              and len(raw)==row.get('bytes'),'retained input bytes '+local[7:])
        check(row.get('ref')==UPSTREAM_REF,'retained immutable upstream ref '+local[7:])
    check(seen=={'inputs/'+name for name in PINS},'exact seven retained upstream inputs')
    actual=load_json((root/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_bytes())
    check(actual['status']=='PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ',
          'retained conditional actual background status')
    check(actual['domain']['ell']=='1/1000000','same actual packet scale')
    check(actual['actual_model']['U_coefficients']==['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27'],
          'actual polynomial potential coefficients')
    for key,expected in {'rho':('129','131'),'a':('111/1000','14/125'),
                         'p':('-7/50000','-3/25000'),'eta':('-9/100000','-1/12500')}.items():
        check(tuple(actual['background_intervals'][key]['exact'])==expected,
              'actual all point background box '+key)
    source_module=helper(root,'derive_source_jets')
    source=source_module.derive()
    check(source_module.N>=16 and source_module.SIG==14,
          'outward rational source Taylor derivative margin')
    check(source_module.ELL==F(1,10**6) and source_module.KAPPA==F(1,10**8),
          'same explicitly distinguished defining function scale')
    check(source['bump_normalized_x_derivative_upper'][:2]==['1','16'],
          'global bump Bell recurrence first coefficients')
    check(source['support_ds_upper']=='20960000000000/9'
          and F(source['sqrt_support_ds_upper'])**2>=F(source['support_ds_upper']),
          'transformed source support length and outward square root')
    scalar=helper(root,'derive_scalar_envelopes').derive()
    constants=helper(root,'derive_regularity_constants')
    low=constants.scalar_coercivity()
    joined=constants.join(source)
    cauchy=constants.physical_cauchy(source)
    for record in (scalar,low,joined,cauchy):
        check(record['status']=='PASS' and len(record['named_checks'])==record['named_check_count'],
              'component roster '+str(record['named_checks'][0]))
        for name in record['named_checks']:
            check(True,name)
    return {'schema':'hdblast-weighted-regularity-v1','status':'PASS',
            'named_checks':names,'named_check_count':len(names),
            'exact':{'source_jets':source,'scalar_polynomial_envelopes':scalar['exact'],
                     'scalar_low_coercivity':low['exact'],'actual_high_order_joins':joined['exact'],
                     'physical_second_Cauchy_jet':cauchy['exact']},
            'scope':'Finite rational majorants, sparse polynomial envelopes and exact numerical joins only. '
                    'Retained branch, exact-family existence, continuous coercivity, dotted closure, '
                    'whole-line Fourier regularity and covariant extension arguments are premises. '
                    'Initial bounds apply only to the explicitly distinguished weight member. '
                    'No evolved second-field bound, amplitude radius, common cylinder or physical calibration follows.'}
