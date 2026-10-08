"""Source-free audit of actual exported certificates and authorization/custody.

Only standard-library code is imported. No retained archive is opened, no
production module is imported or executed, and no physical source is rebuilt.
This verifies exported algebra and bound receipts, not original input truth or
the Taylor coefficients' physical-source provenance by a second construction.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial, prod
import os
from pathlib import Path
import re
import stat

BASE = Path('/workspace/hdblast-research-work/continuation-trajectory-20261008')
HERE = BASE/'actual-independent-review/receipt-audit'
FROZEN = BASE/'implementation-frozen-005'
PACKAGE = Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW')
CORE = PACKAGE/'core'
REGISTRATION = PACKAGE/'FULL_REGISTRATION.json'
GO = BASE/'root/public-readback005/PUBLIC_GO.json'
REG_SHA = '33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796'
GO_SHA = 'edded4f67840bc54e95673c5f43b71f6804e373d13211c90ba48ca99f72a8389'
COMMIT = '3fc078553d1c9b10f11099870195c7a317951a71'
SCOPE = 'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE'
SOURCES = ('positive_B','signed_uB')
CAPSULES = tuple(s+'/'+g for s in SOURCES for g in ('coarse','fine'))
SELECTED = ('k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy','u_2.npy','w_2.npy','u_3.npy','w_3.npy')
DECODE_ORDER = ('observation_eta.npy','k.npy','momentum_weights.npy','u_1.npy','w_1.npy','u_2.npy','w_2.npy','u_3.npy','w_3.npy')
STABLE = ('SOURCE_CERTIFICATE.json','TARGET_BUDGETS.json','NODE_CERTIFICATES.jsonl.gz',
          'DATA.json','SCIENCE_SUMMARY.json','SOURCE_ATTEMPTS.jsonl','DECODE_ATTEMPTS.jsonl',
          'INPUT_AUTHENTICATION.json','OUTPUT_READBACK.json')
LIMITS = {'address_space_bytes':536870912,'aggregate_output_bytes':134217728,
          'core_bytes':0,'cpu_seconds':900,'custodian_receipt_reserve_bytes':65536,
          'file_size_bytes':134217728,'process_limit':0,'wall_seconds':900}
DENIED = {'socket','socketpair','connect','bind','listen','accept','accept4','fork',
          'vfork','clone','clone3','execve','execveat','io_uring_setup','io_uring_enter','io_uring_register'}
CONTRACT = {'anchor':'-9/2','arithmetic_bits':1024,'endpoints':['-4','-7/2'],
            'export_bits':96,'export_gate':'1/1000000000000000000','momentum_degree':2048,
            'output_bytes':134217728,'rss_kib':524288,'selected_decode_members':list(SELECTED),
            'source_cells':64,'source_coefficient_bits':512,'source_degree':24,'wall_seconds':900}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pairs(items):
    out = {}
    for key, value in items:
        require(key not in out, 'duplicate JSON key: '+key)
        out[key] = value
    return out


def parse(raw):
    return json.loads(raw, object_pairs_hook=pairs,
                      parse_constant=lambda s: (_ for _ in ()).throw(ValueError('nonfinite JSON '+s)))


def canonical(value):
    if isinstance(value, F):
        return str(value.numerator)+'/'+str(value.denominator)
    if isinstance(value, dict):
        return {k:canonical(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)):
        return [canonical(v) for v in value]
    return value


def safe_relative(text):
    require(type(text) is str and text and '\\' not in text, 'relative path string')
    p = Path(text)
    require(not p.is_absolute() and '..' not in p.parts and p.as_posix()==text and text!='.', 'unsafe relative path')
    return p


def signature(s):
    return (s.st_dev,s.st_ino,s.st_mode,s.st_nlink,s.st_size,s.st_mtime_ns,s.st_ctime_ns)


def capture(path, limit=134217728):
    path = Path(path).absolute()
    for p in (path,*path.parents):
        require(not stat.S_ISLNK(p.lstat().st_mode), 'symlink ancestry: '+str(p))
    fd = os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK|os.O_CLOEXEC)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink==1 and before.st_size<=limit,'bounded unaliased file: '+str(path))
        parts=[];total=0
        while True:
            raw=os.read(fd,1<<20)
            if not raw:break
            total+=len(raw);require(total<=limit,'read limit');parts.append(raw)
        require(signature(before)==signature(os.fstat(fd))==signature(path.lstat()) and total==before.st_size,'file changed: '+str(path))
        return b''.join(parts)
    finally:os.close(fd)


def pin(raw):
    return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}


def load(path):
    return parse(capture(path,8<<20))


def inventory(root):
    root=Path(root); files={};dirs=set();identities=set()
    require(root.is_dir() and not root.is_symlink(),'real inventory root')
    def error(exc):raise exc
    for parent,names,leaves in os.walk(root,followlinks=False,onerror=error):
        for name in names+leaves:
            path=Path(parent)/name; info=path.lstat(); rel=path.relative_to(root).as_posix()
            safe_relative(rel);identity=(info.st_dev,info.st_ino)
            require(identity not in identities,'aliased inventory inode');identities.add(identity)
            require(not stat.S_ISLNK(info.st_mode),'symlink inventory entry')
            if stat.S_ISDIR(info.st_mode):dirs.add(rel)
            else:files[rel]=pin(capture(path))
    return files,dirs


def q(value):
    require(type(value) is str and re.fullmatch(r'-?(0|[1-9][0-9]*)/[1-9][0-9]*',value),'canonical exact rational syntax')
    result=F(value);require(canonical(result)==value,'unreduced/noncanonical ratio')
    return result


def round_up(value):
    scale=1<<512
    return F(-(-value.numerator*scale//value.denominator),scale)


def control_audit():
    regraw=capture(REGISTRATION); goraw=capture(GO)
    require(pin(regraw)['sha256']==REG_SHA and pin(goraw)['sha256']==GO_SHA,'external control hashes')
    reg=parse(regraw);go=parse(goraw)
    require(reg['scope']==SCOPE and reg['schema_version']==1 and reg['contract']==CONTRACT,'registration scope/contract')
    require(go['status']=='PASS_READBACK_PUBLIC_BYTE_GO' and go['scope']==SCOPE and go['freeze_commit']==COMMIT,'public GO scope/status/commit')
    for key in ('decode_go','source_go','independent_review_pass','public_bytes_verified'):
        require(go[key] is True,'GO boolean: '+key)
    require(go['registration_sha256']==REG_SHA and go['registered_file_pins']==reg['files'],'GO source-byte binding')
    require(go['public_visibility']=={'anonymous_immutable_README_verified':True,'repository_private':False},'public visibility evidence')
    require(set(go['real_operations_before_readback'])=={'physical_source_callbacks','retained_array_decodes','stored_endpoint_comparisons','observational_likelihood_evaluations'} and
            all(type(v) is int and v==0 for v in go['real_operations_before_readback'].values()),'pre-readback zero operations')
    require(go['metric_calibration']=='FAIL' and go['full_continuous_certificate']=='UNRESOLVED' and go['higher_dimensional_origin']=='NOT_ESTABLISHED' and
            go['external_peer_review'] is False and go['external_novelty']=='NOT_ASSESSED','unsupported claims preserved')
    remote=go['remote_file_receipts'];require(len(remote)==go['remote_files_verified']==57,'remote count')
    require(remote['FULL_REGISTRATION.json']['sha256']==REG_SHA and remote['FULL_REGISTRATION.json']['bytes']==len(regraw),'remote registration pin')
    for name,p in reg['files'].items():
        require({k:remote['core/'+name][k] for k in ('bytes','sha256')}==p,'remote source pin '+name)
    expected_dirs={parent.as_posix() for name in reg['files'] for parent in safe_relative(name).parents if parent.as_posix()!='.'}
    source_trees={}
    for name,root in [('execution_core',CORE),('frozen_005',FROZEN)]:
        files,dirs=inventory(root)
        require(files==reg['files'] and dirs==expected_dirs,'complete source tree '+name)
        source_trees[name]={'root':str(root),'file_count':len(files),'directory_count':len(dirs)}
    return reg,go,{'registration':pin(regraw),'GO':pin(goraw),'freeze_commit':COMMIT,
                  'source_trees':source_trees,'remote_receipts_checked_offline':57,
                  'new_public_network_readback':False,'source_files':reg['files']}


def source_accounting(directory):
    cert=load(directory/'SOURCE_CERTIFICATE.json'); budgets=load(directory/'TARGET_BUDGETS.json')
    require({k:v for k,v in cert.items() if k!='rows'}=={'schema_version':1,'manufactured':False,'source_count':2,'cell_count_per_source':64,'degree':24},'actual source certificate schema')
    rows=cert['rows'];require(len(rows)==128,'128 source rows')
    half=F(1,128);radius=F(1,8);majorant=F(64);unit=F(1,1<<512)
    tail=majorant*(half/radius)**25/(1-half/radius)
    require(tail==F(1,15*(1<<90)),'analytic tail identity')
    expected_attempts=[];models={s:[] for s in SOURCES}
    for i,row in enumerate(rows):
        source=SOURCES[i//64];j=i%64;center=F(-9,2)+F(2*j+1,128)
        require(row['source']==source and type(row['cell']) is int and row['cell']==j and q(row['center'])==center,'source order/center')
        require(q(row['half'])==half and q(row['radius'])==radius and q(row['cauchy_majorant'])==majorant,'source geometry')
        require(type(row['exp_degree']) is int and row['exp_degree']==200 and type(row['coefficient_bits']) is int and row['coefficient_bits']==512,'coefficient scheme')
        co=[q(c) for c in row['coefficients']]
        require(len(co)==25 and all((c*(1<<512)).denominator==1 for c in co),'25 dyadic512 coefficients')
        require(hashlib.sha256(json.dumps(row['coefficients'],separators=(',',':')).encode()).hexdigest()==row['coefficient_sha256'],'coefficient vector SHA')
        ce=q(row['coefficient_error']);rounding=q(row['source_error_rounding']);error=q(row['source_error'])
        require(q(row['analytic_tail'])==tail and ce>=0 and 0<=rounding<unit and error==tail+ce+rounding and error==round_up(tail+ce),'complete source error rounded upward exactly once')
        models[source].append((center,co,error))
        expected_attempts.append({'attempt':i,'source':source,'center':canonical(center),'scope':'LATER_ENDPOINT_SOURCE64CELLS_DEGREE24'})
    require([parse(line) for line in capture(directory/'SOURCE_ATTEMPTS.jsonl').splitlines()]==expected_attempts,'128 exact ordered source attempts')
    require(set(budgets)==set(SOURCES),'budget source coverage')
    recomputed={}
    for source,models_for_source in models.items():
        require(set(budgets[source])=={'-4','-7/2'},'two budget endpoints')
        for endpoint in (F(-4),F(-7,2)):
            active=[r for r in models_for_source if r[0]+half<=endpoint]
            duration=endpoint+F(9,2)
            mass=sum((2*half*sum((abs(c)*half**j for j,c in enumerate(co)),F(0)) for center,co,e in active),F(0))
            werror=sum((2*half*e for center,co,e in active),F(0))
            uerror=sum((2*half*(endpoint-center)*e for center,co,e in active),F(0))
            x=2*256*duration
            kernel=x**2049/F(factorial(2049))/(1-x/F(2050))
            expected={'endpoint':str(endpoint),'source_polynomial_L1_upper':str(mass),
                      'source_error_U':str(uerror),'source_error_W':str(werror),
                      'kernel_tail_U':str(round_up(mass*duration*kernel)),'kernel_tail_W':str(round_up(mass*kernel)),
                      'arithmetic':'Arb1024bits contained in exported endpoint rectangles',
                      'integration':'finite polynomial moments integrated exactly; degree2048 exponential tail rounded outward to dyadic512'}
            require(budgets[source][str(endpoint)]==expected,'recomputed complete target budget '+source+' '+str(endpoint))
            recomputed[source+'/'+str(endpoint)]={'cells':len(active),**expected,
                'source_rounding_W_addition_strict_upper':canonical(duration*unit),
                'source_rounding_U_addition_strict_upper':canonical(duration**2*unit/2)}
    return {'source_cells_verified':128,'coefficient_vectors_hashed':128,'dyadic_coefficients_checked':3200,
            'source_attempts_verified':128,'analytic_tail':canonical(tail),'budget_vectors_verified':4,
            'coefficient_construction_recomputed':False,'coefficient_error_physical_premise_recomputed':False,
            'target_budgets':recomputed}


def authentication_accounting(directory):
    spec=load(CORE/'provenance_decoder/INPUT_SPEC.json');auth=load(directory/'INPUT_AUTHENTICATION.json')
    require(auth['manufactured'] is False and len(auth['capsules'])==4,'actual input auth mode')
    expected=[];attempts=[];slots=0;member_count=0;node_count=0
    for ci,capsule in enumerate(spec['capsules']):
        ident=capsule['source']+'/'+capsule['grid'];require(ident==CAPSULES[ci],'capsule roster')
        reader={k:v for k,v in capsule['reader_inputspec'].items() if k!='members'}
        reader['members']=[{k:v for k,v in row.items() if k!='payload_bytes'} for row in capsule['reader_inputspec']['members']]
        require(len(reader['members'])==19 and len({m['name'] for m in reader['members']})==19,'19 unique members')
        selected=[m for m in reader['members'] if m['decode'] is True]
        require({m['name'] for m in selected}==set(SELECTED) and capsule['selected_decode_members']==list(SELECTED),'nine selected members')
        n={'coarse':8192,'fine':16384}[capsule['grid']];require(capsule['full_node_count']==n,'full node plan')
        for member in selected:
            times=member['name']=='observation_eta.npy';complex_=member['name'].startswith(('u_','w_'))
            require(member['shape']==[6 if times else n] and member['descr']==('<c32' if complex_ else '<f16'),'selected member schema')
            slots+=prod(member['shape'])*(2 if complex_ else 1)
        digest=hashlib.sha256(json.dumps(reader,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode('ascii')).hexdigest()
        expected.append({'capsule':ident,'file_sha256':reader['file_sha256'],'input_spec_sha256':digest,
                         'authenticated_members':19,'selected_decode_members':list(SELECTED)})
        for name in DECODE_ORDER:attempts.append({'attempt':len(attempts),'capsule':ident,'member':name})
        member_count+=19;node_count+=n
    require(auth['capsules']==expected,'input authentication canonical metadata bindings')
    require([parse(line) for line in capture(directory/'DECODE_ATTEMPTS.jsonl').splitlines()]==attempts,'36 exact ordered decode attempts')
    require(slots==688152 and member_count==76 and node_count==49152,'derived global input coverage')
    return {'authenticated_member_receipts':member_count,'selected_array_attempts_verified':len(attempts),
            'planned_real_slots_derived_from_shapes':slots,'planned_nodes':node_count,
            'canonical_inputspec_digests':[r['input_spec_sha256'] for r in expected],
            'all_member_authentication_before_decode':'Static worker/codec sequencing plus authenticated successful execution; no second archive read.',
            'original_represented_values_independently_checked':False}


def verify_entry_bindings(entry):
    require(entry['status']=='PASS_REGISTERED_ENDPOINT_ENTRY' and entry['manufactured'] is False and entry['scope']==SCOPE,'actual entry scope')
    require(entry['registration_sha256']==REG_SHA and entry['go_sha256']==GO_SHA,'entry external controls')


def run_audit(directory,optimized,reg):
    directory=Path(directory).absolute();files,dirs=inventory(directory)
    require(not dirs and set(files)==set(STABLE)|{'RESOURCE.json','ENTRY_RECEIPT.json','EXECUTION.json','child.log'},'complete flat actual output roster')
    entry=load(directory/'ENTRY_RECEIPT.json'); execution=load(directory/'EXECUTION.json')
    readback=load(directory/'OUTPUT_READBACK.json');summary=load(directory/'SCIENCE_SUMMARY.json');resource=load(directory/'RESOURCE.json')
    verify_entry_bindings(entry)
    require(entry['outputs']=={n:files[n] for n in (*STABLE,'RESOURCE.json')},'ten exact entry output pins')
    require(execution['status']=='PASS_BOUNDED_REGISTERED_EXECUTION' and execution['manufactured_only'] is False and execution['optimized'] is optimized,'actual bounded mode')
    require(execution['limits']==LIMITS and all(type(v) is int for v in execution['limits'].values()),'exact hard resource policy')
    for field in ('exit_code','wait_status','wrapper_exit_code'):
        require(type(execution[field]) is int and execution[field]==0,'zero '+field)
    require(type(execution['uid']) is int and execution['uid']>0,'nonroot custody')
    for field in ('stop_reason','monitor_error','preserved_child_log'):require(execution[field] is None,'clean custody '+field)
    require(execution['entry_receipt']==files['ENTRY_RECEIPT.json'],'custodian entry SHA')
    require(execution['captured_child_log_sha256']==execution['observed_child_log_sha256']==files['child.log']['sha256'],'held child log SHA')
    require(execution['aggregate_limit_enforcement']=='POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE','aggregate cap semantics')
    for field,cap in [('wall_seconds',900),('cpu_seconds_wait4',900),('peak_rss_kib_wait4',524288)]:
        require(type(execution[field]) in (int,float) and 0<=execution[field]<=cap,'authoritative resource '+field)
    command=execution['command'];require(type(command) is list and all(type(v) is str for v in command),'recorded worker command')
    prefix=['-I','-B']+(['-O'] if optimized else []);script_index=1+len(prefix)
    require(command[1:script_index]==prefix and Path(command[0]).is_absolute(),'isolated actual Python command')
    args=command[script_index+1:]
    require(args[::2]==['--root','--output','--registration','--registration-sha256','--repository-root','--go','--go-sha256'] and len(args)==14,'exact actual worker argument roster')
    arguments=dict(zip(args[::2],args[1::2]));recorded_root=Path(arguments['--root'])
    require(command[script_index]==str(recorded_root/'execution/run_endpoint.py'),'recorded script/source-root binding')
    for option in ('--root','--output','--registration','--repository-root','--go'):
        require(Path(arguments[option]).is_absolute(),'absolute historical command path')
    require(execution['custodian_output_directory']==arguments['--output'],'custodian held output/command binding')
    require(arguments['--registration-sha256']==REG_SHA and arguments['--go-sha256']==GO_SHA,'command external control pins')
    # Historical pathnames remain evidence, not new read locations. Copied saved
    # outputs can therefore be audited at a new path using the same byte pins.
    denied=entry['kernel_denied_syscalls'];require(type(denied) is list and len(denied)==len(DENIED) and set(denied)==DENIED,'full kernel-denial evidence')
    tree=hashlib.sha256();size=0;count=0
    for name,p in sorted(files.items()):
        if name=='EXECUTION.json':continue
        tree.update(name.encode('utf-8',errors='surrogateescape')+b'\0'+str(p['bytes']).encode()+b'\0'+bytes.fromhex(p['sha256']));size+=p['bytes'];count+=1
    require(execution['worker_output']=={'bytes':size,'directory_count':0,'file_count':count,'tree_sha256':tree.hexdigest()},'complete custodian output tree')
    require(size<=134217728-65536 and sum(p['bytes'] for p in files.values())<=134217728,'worker and total output caps')
    require(resource['output_bytes']==sum(p['bytes'] for n,p in files.items() if n not in ('RESOURCE.json','ENTRY_RECEIPT.json','EXECUTION.json')),'worker resource byte snapshot')
    require(0<=resource['wall_seconds']<=900 and 0<=resource['peak_rss_kib']<=524288,'worker resource bounds')
    for receipt in (entry,readback,summary):
        require(type(receipt['nodes']) is int and receipt['nodes']==49152 and type(receipt['endpoint_comparisons']) is int and receipt['endpoint_comparisons']==98304,'node/comparison receipt coverage')
    require(readback['status']=='PASS_STANDALONE_ENDPOINT_READBACK' and readback['manufactured'] is False and readback['prefix_cases']==24 and readback['source_coefficient_vectors']==128,'readback scope/coverage')
    for key in ('source_evaluations','retained_arrays_opened'):
        require(type(readback[key]) is int and readback[key]==0,'readback zero extra operations')
    require(summary['status']=='FINITE_RETAINED_ENDPOINT_ERRORS_ENCLOSED' and summary['selected_arrays']==36 and summary['source_callbacks']==128 and summary['decoded_real_slots']==688152,'science receipt coverage')
    require(summary['incoming_error_counted'] is False and summary['prior_metric_calibration']=='FAIL_UNCHANGED' and summary['original_unsaved_continuous_solver_trajectory']=='UNRECOVERABLE_FROM_RETAINED_MODE_SNAPSHOTS' and summary['continuous_momentum_pressure_contact_time_UV_certificate']=='UNRESOLVED','science scope restrictions')
    require(readback['maximum_complete_endpoint_export_L1_radius']==summary['maximum_complete_endpoint_export_L1_radius'],'summary/readback radius agreement')
    require(all(0<=F(v)<=F(1,10**18) for v in summary['maximum_complete_endpoint_export_L1_radius'].values()),'complete export radius bound')
    return {'directory':str(directory),'mode':'optimized' if optimized else 'normal','recorded_custodian_directory':execution['custodian_output_directory'],
            'recorded_command':command,
            'source_accounting':source_accounting(directory),'input_authentication':authentication_accounting(directory),
            'wall_seconds':execution['wall_seconds'],'cpu_seconds_wait4':execution['cpu_seconds_wait4'],
            'peak_rss_kib_wait4':execution['peak_rss_kib_wait4'],'worker_output':execution['worker_output'],
            'file_pins':files,'stable_scientific_files':{name:files[name] for name in STABLE}}


def main():
    global CORE,FROZEN,REGISTRATION,GO,HERE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core',type=Path,default=CORE)
    parser.add_argument('--frozen-source',type=Path,default=FROZEN)
    parser.add_argument('--registration',type=Path,default=REGISTRATION)
    parser.add_argument('--go',type=Path,default=GO)
    parser.add_argument('--normal',type=Path,default=BASE/'root/actual-normal005-001')
    parser.add_argument('--optimized',type=Path)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    args=parser.parse_args()
    CORE=args.core.absolute();FROZEN=args.frozen_source.absolute()
    REGISTRATION=args.registration.absolute();GO=args.go.absolute();HERE=args.output_dir.absolute()
    require(HERE.is_dir() and not HERE.is_symlink(),'existing real review output directory')
    reg,go,controls=control_audit()
    normal=run_audit(args.normal,False,reg);runs={'normal':normal}
    if args.optimized:
        optimized=run_audit(args.optimized,True,reg)
        require(normal['stable_scientific_files']==optimized['stable_scientific_files'],'nine normal/optimized scientific artifacts byte-identical')
        runs['optimized']=optimized
    result={'status':'PASS_INDEPENDENT_ACTUAL_RECEIPTS_AND_EXPORTED_SOURCE_ACCOUNTING',
            'paired_scientific_identity_checked':bool(args.optimized),'stable_files_compared':9 if args.optimized else 0,
            'controls':controls,'runs':runs,'audit_script':pin(capture(Path(__file__))),
            'audit_operations':{'production_modules_imported':0,'retained_arrays_opened':0,'retained_values_decoded':0,
                                'physical_source_callbacks':0,'producer_executions':0,'network_requests':0,'GO_issued':False},
            'limits_of_review':[
                'Original represented values cannot be independently checked without the prohibited duplicate retained-array decode.',
                'Exported coefficient/error algebra is independently checked; physical coefficients and the coefficient-error premise are not reconstructed.',
                'Reported complete input authentication relies on frozen worker/codec sequencing and successful pinned custody; archives are not reopened.',
                'This audit does not independently evolve the target flow or verify all-node arithmetic; that is addressed separately.',
                'Normal/optimized byte identity is a control diagnostic, not an independent mathematical proof.',
                'GO public-readback evidence is checked offline against the pinned receipt; this audit performs no second public network readback.']}
    name='ACTUAL_PAIR_RECEIPT_AUDIT.json' if args.optimized else 'ACTUAL_NORMAL_RECEIPT_AUDIT.json'
    (HERE/name).write_text(json.dumps(canonical(result),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'receipt':str(HERE/name),'sha256':pin(capture(HERE/name))['sha256'],
                      'runs':list(runs),'source_rows_per_run':128,'arrays_per_run':36,'slots_per_run':688152},sort_keys=True))


if __name__=='__main__':main()
