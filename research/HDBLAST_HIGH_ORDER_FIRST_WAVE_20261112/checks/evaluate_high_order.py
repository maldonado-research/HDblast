"""High-order first wave finite replay, authenticated original source snapshots."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF='b842e66402e24cf1ebe3cc278364f6062d395cce'
REPOSITORY='https://github.com/maldonado-research/HDblast'


def need(condition,message):
    if not condition:
        raise ValueError(message)


def public_sources(root,check):
    pins=json.loads((root/'SOURCE_PINS.json').read_text(encoding='utf-8'))
    need(pins.get('schema')=='hdblast-high-order-first-wave-pins-v1','Public pin schema mismatch')
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
    need(set(data)=={'schema','symbolic','bounds','scope'} and data['schema']=='HDBLAST_HIGH_ORDER_FIRST_WAVE_DATA_V1','Data schema mismatch')
    for key,value in data['bounds'].items():
        if key.endswith('_decade') or key.endswith('_order'):
            need(type(value) is int and value>=0,'Nonnegative integer budget required: '+key)
    for key,value in {'raw_coefficient_order':16,'coefficient_decade':152000,'initial_input_decade':5200,'base_K_decade':150006,'prior_forcing_decade':304240}.items():
        need(data['bounds'].get(key)==value,'Retained analytic input mismatch: '+key)
    suites={s:[] for s in ('sources','rotations','bounds')}
    def checker(suite):
        def check(name,condition):
            need(condition,'Arithmetic failed: '+suite+':'+name)
            need(name not in suites[suite],'Duplicate check')
            suites[suite].append(name)
        return check
    count=public_sources(root,checker('sources'))
    exact={'rotations':load_module(root,'formal_rotations').evaluate(data['symbolic'],checker('rotations')),
           'bounds':load_module(root,'high_order_bounds').evaluate(data['bounds'],checker('bounds'))}
    names=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-high-order-first-wave-v1','status':'PASS','named_check_count':len(names),'named_checks':names,
      'suite_counts':{s:len(entries) for s,entries in suites.items()},'source_count':count,'exact':exact,'scope':data['scope']}
