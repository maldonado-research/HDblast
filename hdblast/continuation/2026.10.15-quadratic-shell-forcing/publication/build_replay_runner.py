"""Embed an explicitly approved build-time configuration; never execute scientific code."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path, PurePosixPath
import pprint
import re

CAPSULE='HDBLAST_QUADRATIC_SHELL_FORCING_20261015'

def require(ok,message):
    if not ok:raise RuntimeError(message)

def canonical(name):
    require(isinstance(name,str) and name and '\\' not in name,'Relative string path required')
    p=PurePosixPath(name)
    require(not p.is_absolute() and '..' not in p.parts and p.as_posix()==name,'Canonical relative path required')
    return name

def json_path(path):
    require(isinstance(path,list) and path and all((type(k) is int and k>=0) or (isinstance(k,str) and k) for k in path),'Nonempty JSON path required')
    return path

def load_template(path):
    raw=path.read_bytes();spec=importlib.util.spec_from_file_location('replay_template',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return raw,module

def build(configuration,source,template):
    require(configuration['schema']=='hdblast-build-time-replay-config-v1'
            and configuration['status']=='APPROVED_FOR_FROZEN_REPLAY'
            and configuration['capsule']==CAPSULE,'Explicit current approved build-time configuration required')
    require(re.fullmatch(r'3\.\d+\.\d+',configuration['python_version']) is not None,'Exact tested Python version required')
    dependencies=configuration['dependency_versions']
    require(isinstance(dependencies,dict) and all(re.fullmatch(r'[A-Za-z0-9_.-]+',k) and isinstance(v,str) and v for k,v in dependencies.items()),'Exact dependency versions required')
    source=source.resolve();pins={}
    def read(name):
        canonical(name);p=source/name
        require(p.is_file() and not any(q.is_symlink() for q in [p,*p.parents]) and p.resolve().is_relative_to(source),'Real source dependency required: '+name)
        raw=p.read_bytes();pin={'path':name,'bytes':len(raw),'sha256':sha256(raw).hexdigest()}
        require(name not in pins or pins[name]==pin,'Source changed during build-time inspection')
        pins[name]=pin;return raw
    checkers=[]
    for rawchecker in configuration['checkers']:
        c=dict(rawchecker);name=c['name']
        require(re.fullmatch('[a-z][a-z0-9_]*',name) is not None,'Safe unique checker name required')
        c['script']=canonical(c['script']);c.setdefault('dependencies',[])
        require(isinstance(c['dependencies'],list) and len(set(c['dependencies']))==len(c['dependencies']) and c['script'] not in c['dependencies'],'Unique dependency list excluding checker required')
        for dependency in [c['script'],*c['dependencies']]:read(dependency)
        c.setdefault('arguments',[]);require(all(isinstance(a,str) for a in c['arguments']),'String arguments required')
        c.setdefault('output',None)
        if c['output'] is not None:canonical(c['output'])
        c.setdefault('additional_outputs',[])
        outputs=[c['output']] if c['output'] else []
        for pair in c['additional_outputs']:
            require(isinstance(pair,list) and len(pair)==2,'Additional [output,golden] pair required');canonical(pair[0]);read(pair[1]);outputs.append(pair[0])
        require(len(outputs)==len(set(outputs)) and not set(outputs)&set([c['script'],*c['dependencies']]) and not set(outputs)&{'fresh.json','execution.log'},'Generated outputs must be fresh, distinct, nonreserved paths')
        c.setdefault('byte_equal',False);require(type(c['byte_equal']) is bool,'Boolean byte_equal required')
        c.setdefault('status_path',['status']);json_path(c['status_path']);require(isinstance(c['status'],str) and c['status'],'Explicit status required')
        c.setdefault('mode_fields',[])
        fieldpaths=[]
        for field in c['mode_fields']:
            json_path(field['path']);fieldpaths.append(field['path'])
            require(type(field['normal']) in (int,bool,str) and type(field['normal']) is type(field['optimized']) and field['normal']!=field['optimized'],'Explicit typed distinct execution-mode values required')
        require(len({tuple(x) for x in fieldpaths})==len(fieldpaths),'Unique execution-mode paths required')
        require(not any(a!=b and a==b[:len(a)] for a in fieldpaths for b in fieldpaths),'Overlapping execution-mode fields')
        c.setdefault('checks',None);c.setdefault('embedded_checks',None)
        for key in ('checks','embedded_checks'):
            if c[key] is not None:
                json_path(c[key]['path']);require(c[key]['kind'] in ('names','true_map'),'Supported named inventory kind required')
        c.setdefault('required_fields',[])
        for row in c['required_fields']:json_path(row['path']);require('equals' in row,'Required-field comparison value missing')
        for key in ('check_count_path','embedded_count_path'):
            c.setdefault(key,None)
            if c[key] is not None:json_path(c[key])
        normal=json.loads(read(c['golden']));optimized=json.loads(read(c.get('optimized_golden',c['golden'])))
        c['check_names']=template.inventory_names(normal,c['checks']);c['embedded_control_names']=template.inventory_names(normal,c['embedded_checks'])
        c['groups']=len(c['check_names']);c['embedded']=len(c['embedded_control_names'])
        for mode,data in [('normal',normal),('optimized',optimized)]:template.positive_inventory(c,data,mode)
        for field in c['mode_fields']:
            for data in (normal,optimized):
                parent=data
                for key in field['path'][:-1]:parent=parent[key]
                del parent[field['path'][-1]]
        require(template.json_equal(normal,optimized),'Complete golden scientific outputs differ by mode: '+name)
        checkers.append(c)
    named={c['name']:c for c in checkers}
    require(named and len(named)==len(checkers) and sum(c['groups'] for c in checkers)>0,'Unique positive checkers with actual enumerated groups required')
    controls=[];ids=set()
    for original in configuration['controls']:
        c=dict(original);key=(c['checker'],c['name'])
        require(c['checker'] in named and key not in ids and re.fullmatch('[A-Za-z][A-Za-z0-9_-]*',c['name']) is not None,'Unique control for known checker required');ids.add(key)
        require(isinstance(c['reason'],str) and c['reason'] and isinstance(c['mechanism'],str) and c['mechanism'],'Explicit control reason/mechanism required')
        require(c['failure_format'] in ('fail_line','traceback','json'),'Exact failure diagnostic format required')
        if c['failure_format']=='traceback':
            require(re.fullmatch('[A-Za-z][A-Za-z0-9_.]*',c['exception']) is not None,'Expected exception class required')
            prefix=c.get('diagnostic_prefix','')
            require(isinstance(prefix,str) and (not prefix or prefix.endswith('\n')), 'Exact optional diagnostic prefix required')
        if c['failure_format']=='json':json_path(c['failure_path']);json_path(c.get('status_path',['status']))
        if 'arguments' in c:require(isinstance(c['arguments'],list) and all(isinstance(x,str) for x in c['arguments']),'Control string arguments required')
        checker=named[c['checker']]
        if 'target' in c:
            require(c['target'] in checker['dependencies'] and c['target']!=checker['script'],'Corruption must target copied dependency, retaining independent checker')
            recipes=int('append' in c)+int('json_path' in c)+int('before' in c)
            require(recipes==1,'Exactly one explicit input-mutation recipe required')
            raw=read(c['target'])
            if 'append' in c:require(isinstance(c['append'],str) and c['append'],'Nonempty append required')
            elif 'json_path' in c:
                json_path(c['json_path']);template.at_path(json.loads(raw),c['json_path'])
                require(('value' in c)^('factor' in c),'Exactly one JSON replacement value or rational factor required')
                if 'factor' in c:require(template.F(c['factor'])!=1,'Nontrivial rational factor required')
            else:require(isinstance(c['before'],str) and raw.decode().count(c['before'])==1 and c['after']!=c['before'],'Unique changed formula replacement required')
        else:require(c.get('arguments') and c['arguments']!=checker['arguments'],'Distinct production CLI injection arguments required')
        controls.append(c)
    require(controls,'Meaningful external controls required')
    # Recheck all source bytes after reading both goldens and data dependencies.
    for name,pin in list(pins.items()):require(sha256(read(name)).hexdigest()==pin['sha256'],'Source input changed while building')
    counts={'successful_checker_runs':2*len(checkers),'external_corruption_processes_rejected':2*len(controls),'embedded_wrong_formula_controls_rejected':2*sum(c['embedded'] for c in checkers),'finite_positive_check_groups_each_mode':sum(c['groups'] for c in checkers)}
    roster={'schema':'hdblast-replay-roster-v1','status':'APPROVED_FOR_FROZEN_REPLAY','capsule':CAPSULE,'modes':['normal','optimized'],'counts':counts,
      'positive_processes':[{'name':c['name'],'script':c['script'],'finite_positive_groups_each_mode':c['groups'],'embedded_controls_each_mode':c['embedded'],'finite_positive_check_names':c['check_names'],'embedded_control_names':c['embedded_control_names']} for c in checkers],
      'external_controls':[{'checker':c['checker'],'name':c['name'],'expected_reason':c['reason'],'mechanism':c['mechanism']} for c in controls]}
    return checkers,controls,roster,sorted(pins.values(),key=lambda r:r['path'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True);parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--replace',action='store_true')
    args=parser.parse_args();here=Path(__file__).resolve().parent
    raw=args.config.read_bytes();configuration=json.loads(raw);template_raw,template=load_template(here/'runner_template.py')
    checkers,controls,roster,pins=build(configuration,args.source_root,template)
    constants={'CAPSULE':CAPSULE,'PYTHON_VERSION':configuration['python_version'],'DEPENDENCY_VERSIONS':configuration['dependency_versions'],'CONFIG_SHA256':sha256(raw).hexdigest(),'CHECKERS':checkers,'CONTROLS':controls}
    source=template_raw.decode();start=source.index('# BEGIN EMBEDDED APPROVED CONFIGURATION');end=source.index('# END EMBEDDED APPROVED CONFIGURATION')+len('# END EMBEDDED APPROVED CONFIGURATION')
    source=source[:start]+'# BEGIN EMBEDDED APPROVED CONFIGURATION\n'+'\n'.join(k+' = '+pprint.pformat(v,sort_dicts=False,width=100) for k,v in constants.items())+'\n# END EMBEDDED APPROVED CONFIGURATION'+source[end:]
    compile(source,'reproduce.py','exec')
    output=args.output_dir.absolute();require(not any(p.is_symlink() for p in [output,*output.parents]),'Real output directory required')
    require(not output.resolve().is_relative_to(args.source_root.resolve()) or output.resolve().is_relative_to(here.resolve()),'Build only into publication or separate scratch; do not overwrite science')
    files={'reproduce.py':source.encode(),'REPLAY_ROSTER.json':(json.dumps(roster,indent=2)+'\n').encode()}
    for name in files:require(args.replace or not (output/name).exists(),'Refuse existing output without --replace')
    output.mkdir(parents=True,exist_ok=True)
    for name,data in files.items():(output/name).write_bytes(data)
    receipt={'status':'PASS_BUILD_TIME_CONFIGURATION_EMBEDDED','config':{'bytes':len(raw),'sha256':sha256(raw).hexdigest()},'template_sha256':sha256(template_raw).hexdigest(),'builder_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'source_input_pins':pins,'counts':roster['counts'],'outputs':[{'path':k,'bytes':len(v),'sha256':sha256(v).hexdigest()} for k,v in files.items()],'scientific_processes_executed':0,'scope':'Source generation and complete existing-golden inspection only; no fresh scientific replay or publication certified.'}
    (output/'RUNNER_BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'counts':roster['counts'],'outputs':receipt['outputs']},indent=2))

if __name__=='__main__':main()
