"""Ordered evolution finite replay, authenticated original source snapshots."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF='73b60c4f4e7013c3d4c682553274825ebf14888b'
REPOSITORY='https://github.com/maldonado-research/HDblast'


def need(condition,message):
    if not condition:
        raise ValueError(message)


def public_sources(root,check):
    pins=json.loads((root/'SOURCE_PINS.json').read_text(encoding='utf-8'))
    need(pins.get('schema')=='hdblast-ordered-evolution-pins-v1','Public pin schema mismatch')
    need(pins.get('repository')==REPOSITORY and pins.get('ref')==REF,
         'Public source repository/ref mismatch')
    files=pins.get('files')
    need(isinstance(files,list) and len(files)==8,'Eight original public inputs required')
    names=set()
    for row in files:
        local=PurePosixPath(row['local'])
        source=PurePosixPath(row['path'])
        need(not local.is_absolute() and len(local.parts)==2 and local.parts[0]=='inputs'
             and local.name not in names and '\\' not in row['local']
             and local.as_posix()==row['local'] and '..' not in local.parts,
             'Unsafe or repeated local source pin')
        need(not source.is_absolute() and source.parts[0]=='research'
             and '\\' not in row['path'] and '..' not in source.parts
             and source.as_posix()==row['path'],'Unsafe original source path')
        need(row['ref']==REF and row['url']==REPOSITORY+'/blob/'+REF+'/'+source.as_posix(),
             'Immutable public URL mismatch')
        raw=(root/local).read_bytes()
        check('source_size:'+local.name,type(row['bytes']) is int and len(raw)==row['bytes'])
        check('source_sha256:'+local.name,isinstance(row['sha256'],str)
              and re.fullmatch('[0-9a-f]{64}',row['sha256']) is not None
              and hashlib.sha256(raw).hexdigest()==row['sha256'])
        blob=b'blob '+str(len(raw)).encode('ascii')+b'\0'+raw
        check('source_git_blob:'+local.name,hashlib.sha1(blob).hexdigest()==row['git_blob_sha1'])
        names.add(local.name)
    need({p.name for p in (root/'inputs').iterdir()}==names,'Public input roster mismatch')
    return len(names)


def load_module(root,name):
    spec=importlib.util.spec_from_file_location('ordered_'+name,root/'checks'/(name+'.py'))
    need(spec is not None and spec.loader is not None,'Evaluator loader unavailable')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(root):
    data=json.loads((root/'DATA.json').read_text(encoding='utf-8'))
    need(set(data)=={'schema','symbolic','bounds','forcing','scope'} and
         data['schema']=='HDBLAST_ORDERED_EVOLUTION_DATA_V1','Data schema mismatch')
    for section in ('bounds','forcing'):
        for key,value in data[section].items():
            if key.endswith('_decade') or key.endswith('_order') or key.endswith('_factor'):
                need(type(value) is int and value>=0,'Nonnegative integer budget required: '+key)
    expected={
        'coefficient_decade':150000,'state_dimension':16,
        'bundle_connection_bound':'2','bundle_connection_derivative_bound':'4',
        'base_connection_entrysum_bound':'4','inverse_entrysum_bound':'1',
        'fiber_norm_conversion':'2','moving_frame_conversion':'4',
        'inverse_aggregate_decade':17000,'initial_first_decade':5200,
    }
    for key,value in expected.items():
        need(data['bounds'].get(key)==value,'Retained input premise mismatch: '+key)
    need(data['forcing']['coefficient_decade']==data['bounds']['coefficient_H12_decade'],
         'Forcing coefficient interface differs')
    need(data['forcing']['initial_N14_decade']==data['bounds']['initial_first_decade'],
         'Forcing initial norm interface differs')
    suites={name:[] for name in ('sources','symbolic','bounds','forcing')}
    def checker(suite):
        def check(name,condition):
            need(condition,'Arithmetic failed: '+suite+':'+name)
            need(isinstance(name,str) and name and name not in suites[suite],'Invalid or repeated check')
            suites[suite].append(name)
        return check
    count=public_sources(root,checker('sources'))
    coefficients=(root/'inputs/COEFFICIENTS_06G.md').read_text(encoding='utf-8')
    initial=(root/'inputs/INITIAL_TRANSPORT_06H.md').read_text(encoding='utf-8')
    checker('bounds')('retained_coefficient_literal_premise','C=10^150000, |α|≤16.' in coefficients)
    checker('bounds')('retained_initial_first_literal_premise','< 10^5200' in initial)
    exact={
        'symbolic':load_module(root,'formal_ordered').evaluate(data['symbolic'],checker('symbolic')),
        'bounds':load_module(root,'ordered_bounds').evaluate(root,data['bounds'],checker('bounds')),
        'forcing':load_module(root,'forcing_bounds').evaluate(data['forcing'],checker('forcing')),
    }
    names=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-ordered-evolution-v1','status':'PASS',
        'named_check_count':len(names),'named_checks':names,
        'suite_counts':{suite:len(entries) for suite,entries in suites.items()},
        'source_count':count,'exact':exact,'scope':data['scope']}
