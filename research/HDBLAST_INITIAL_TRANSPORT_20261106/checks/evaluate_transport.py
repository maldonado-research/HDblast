"""Initial transport finite replay, authenticated original source snapshots."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF='66ce4c07cd9b16cb67f5e9f147402d16781c676e'
REPOSITORY='https://github.com/maldonado-research/HDblast'


def need(condition,message):
    if not condition:
        raise ValueError(message)


def public_sources(root,check):
    pins=json.loads((root/'SOURCE_PINS.json').read_text(encoding='utf-8'))
    need(pins.get('schema')=='hdblast-initial-transport-pins-v1','Public pin schema mismatch')
    need(pins.get('repository')==REPOSITORY and pins.get('ref')==REF,
         'Public source repository/ref mismatch')
    files=pins.get('files')
    need(isinstance(files,list) and len(files)==9,'Nine original public inputs required')
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
    spec=importlib.util.spec_from_file_location('transport_'+name,root/'checks'/(name+'.py'))
    need(spec is not None and spec.loader is not None,'Evaluator loader unavailable')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(root):
    data=json.loads((root/'DATA.json').read_text(encoding='utf-8'))
    need(set(data)=={'schema','symbolic','bounds','scope'} and
         data['schema']=='HDBLAST_INITIAL_TRANSPORT_DATA_V1','Data schema mismatch')
    for key,value in data['bounds'].items():
        if key.endswith('_decade') or key.endswith('_order') or key.endswith('_factor'):
            need(type(value) is int and value>=0,'Nonnegative integer budget required: '+key)
    need(data['bounds']['flow_total_order']==18 and data['bounds']['flow_decade']==126,
         'Retained flow input mismatch')
    need(data['bounds']['raw_second_decade']==4000,'Retained raw W2 input mismatch')
    need(data['bounds']['conversion_branch_factor']>=10000,'Conversion branch factor too small')
    need(data['bounds']['normalization_decade']>=210,'Normalization factor budget too small')
    need(data['bounds']['warped_background_y_decade']>=80 and
         data['bounds']['background_z_decade']>=100 and data['bounds']['angular_decade']>=200,
         'Retained geometric premise budget too small')
    suites={name:[] for name in ('sources','symbolic','bounds')}
    def checker(suite):
        def check(name,condition):
            need(condition,'Arithmetic failed: '+suite+':'+name)
            need(isinstance(name,str) and name and name not in suites[suite],'Invalid or repeated check')
            suites[suite].append(name)
        return check
    count=public_sources(root,checker('sources'))
    gauge=(root/'inputs/GAUGE_CAUCHY_COMPLETION_06G.md').read_text(encoding='utf-8')
    flow=(root/'inputs/BOUNDARY_FLOW_PROOF_06F.md').read_text(encoding='utf-8')
    checker('bounds')('retained_raw_W2_literal_premise','W2 <10^3641 <10^4000 per copy.' in gauge)
    checker('bounds')('retained_flow_literal_premise','<10^126,     m+k<=18.' in flow)
    exact={
        'symbolic':load_module(root,'formal_transport').evaluate(data['symbolic'],checker('symbolic')),
        'bounds':load_module(root,'transport_bounds').evaluate(root,data['bounds'],checker('bounds')),
    }
    names=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-initial-transport-v1','status':'PASS',
        'named_check_count':len(names),'named_checks':names,
        'suite_counts':{suite:len(entries) for suite,entries in suites.items()},
        'source_count':count,'exact':exact,'scope':data['scope']}
