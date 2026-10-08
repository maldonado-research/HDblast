"""Registered all-node endpoint worker; import only after outer byte guard."""
from contextlib import ExitStack
from dataclasses import asdict
from fractions import Fraction as Q
from itertools import zip_longest
from pathlib import Path
import gzip
import hashlib
import io
import json
import math
import time
import zipfile
from registration_guard import require,safe_file,load,sha,reauthenticate,resource_guard
import exact_binary80 as codec
import fabricated_fixtures as fixture
from input_spec_adapter import hydrate_reader_spec,SELECTED_MEMBERS,metadata_receipt
from endpoint_engine import EndpointFlow,ENDPOINTS,SCALE
from endpoint_aggregate import EPSILON,PI,Accumulator,node_terms
from later_source import build_models
from kernel_guard import real_directory,require_worker_limits

CAPSULES=tuple(s+'/'+g for s in ('positive_B','signed_uB') for g in ('coarse','fine'))
OBSERVATIONS=(Q(-11,2),Q(-9,2),Q(-4),Q(-7,2),Q(-5,2),Q(-3,2))
PREFIXES={'coarse':{64:2048,128:4096,256:8192},'fine':{64:4096,128:8192,256:16384}}

def canonical(value):
    if type(value) is Q:return str(value.numerator)+'/'+str(value.denominator)
    if type(value) is dict:return {k:canonical(v) for k,v in value.items()}
    if type(value) in (list,tuple):return [canonical(v) for v in value]
    return value

def json_bytes(value):return (json.dumps(canonical(value),sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def write_json(path,value):
    with Path(path).open('xb') as out:out.write(json_bytes(value))

def whole_hash(path,bytes_expected,digest):
    require(path.stat().st_size==bytes_expected,'opaque input size differs')
    h=hashlib.sha256();count=0
    with path.open('rb') as source:
        while block:=source.read(1<<20):
            count+=len(block);require(count<=bytes_expected,'opaque input grew');h.update(block)
    require(count==bytes_expected and h.hexdigest()==digest,'opaque input complete hash differs')

def inventory(path,pin):
    with zipfile.ZipFile(path) as z:
        require(z.comment==b'','archive comment differs')
        found=[{'bytes':i.file_size,'compressed_bytes':i.compress_size,'compression':i.compress_type,
            'crc32_hex':f'{i.CRC:08x}','create_system':i.create_system,'external_attr':i.external_attr,
            'extra_hex':i.extra.hex(),'flag_bits':i.flag_bits,'header_offset':i.header_offset,
            'name':i.filename,'timestamp':list(i.date_time)} for i in z.infolist()]
        require(found==pin,'exact ordered ZIP metadata differs')

def authenticate_inputs(v,repo,document,stack,output):
    require(v['manufactured'] is False and v['go'] is not None,'actual decode GO required')
    repo=Path(repo).absolute();require(repo.is_dir() and not repo.is_symlink(),'real retained repository required')
    for label in ('producer','manifest'):
        path=safe_file(repo,document[label+'_repository_path'])
        whole_hash(path,path.stat().st_size,document[label+'_sha256'])
    for label,name in (('layout_evidence','LAYOUT_PROVENANCE_BASIS.json'),
                       ('original_zip_inventories','ORIGINAL_ZIP_INVENTORIES.json')):
        pin=document[label]
        whole_hash(safe_file(v['root'],'provenance_decoder/'+name),pin['bytes'],pin['sha256'])
    static=document['prior_static_audit']
    whole_hash(safe_file(repo,static['repository_path']),static['bytes'],static['sha256'])
    limits=codec.Limits(**document['explicit_limits']);snapshots={};proof=[]
    for c in document['capsules']:
        ident=c['source']+'/'+c['grid'];spec=hydrate_reader_spec(c['reader_inputspec'])
        original=c['original'];whole_hash(safe_file(repo,original['repository_path']),original['bytes'],original['sha256'])
        path=safe_file(repo,c['repository_path']);whole_hash(path,spec.file_bytes,spec.file_sha256)
        inventory(path,c['strict_zip_inventory'])
        snapshots[ident]=stack.enter_context(codec.verify_inputs(path,spec,limits,temporary_directory=output))
        proof.append({'capsule':ident,'file_sha256':spec.file_sha256,
                      'input_spec_sha256':snapshots[ident].input_spec_sha256,'authenticated_members':19,
                      'selected_decode_members':list(SELECTED_MEMBERS)})
    require(tuple(snapshots)==CAPSULES,'all four opaque snapshots required before decode')
    reauthenticate(v)
    return snapshots,proof

def slot(q,padding):
    require(type(q) is Q and q.denominator&(q.denominator-1)==0,'exact manufactured dyadic required')
    if q==0:return fixture.slot(0,0,False,padding)
    a=abs(q);e=a.numerator.bit_length()-1-(a.denominator.bit_length()-1)
    significand=a*Q(2)**(63-e)
    require(significand.denominator==1 and 1<<63<=significand<1<<64,'manufactured binary80 precision overflow')
    return fixture.slot(significand.numerator,e+16383,q<0,padding)

def manufacture_inputs(document,stack,output):
    directory=output/'manufactured-inputs';directory.mkdir();snapshots={};proof=[]
    limits=codec.Limits(**document['explicit_limits'])
    for c in document['capsules']:
        sid=0 if c['source']=='positive_B' else 1;n=c['full_node_count'];den=64 if c['grid']=='coarse' else 128
        padding=bytes([sid+1])*6;entries=[]
        for old in c['reader_inputspec']['members']:
            name,descr,shape=old['name'],old['descr'],tuple(old['shape'])
            if name=='k.npy':payload=b''.join(slot(Q(2*i+1,den),padding) for i in range(n))
            elif name=='momentum_weights.npy':payload=slot(Q(2,den),padding)*n
            elif name=='observation_eta.npy':payload=b''.join(slot(q,padding) for q in OBSERVATIONS)
            elif name.startswith(('u_','w_')):
                suffix=int(name[2]);divisor=1<<30 if name.startswith('u_') else 1<<31
                payload=b''.join(slot(Q((sid+1)*suffix*(i%17+1),divisor),padding)+
                                  slot(Q(-(sid+1)*suffix*(i%13+1),divisor),padding) for i in range(n))
            else:payload=b'\0'*(math.prod(shape)*{'<f16':16,'<c32':32,'<i8':8,'<U8':32}[descr])
            data=fixture.npy(payload,descr,shape);entries.append((name,data,descr,shape,old['decode']))
        buffer=io.BytesIO();members=[]
        with zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_STORED) as z:
            for name,data,descr,shape,decode in entries:
                info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_STORED
                z.writestr(info,data);members.append(fixture.member(name,data,descr,shape,decode,0))
        raw=buffer.getvalue();spec=fixture.input_spec(raw,members,'npz')
        path=directory/(c['source']+'_'+c['grid']+'.npz')
        with path.open('xb') as out:out.write(raw)
        ident=c['source']+'/'+c['grid']
        snapshots[ident]=stack.enter_context(codec.verify_inputs(path,spec,limits,temporary_directory=directory))
        proof.append({'capsule':ident,'manufactured':True,'file_sha256':spec.file_sha256,
                      'input_spec_sha256':snapshots[ident].input_spec_sha256,'reader_inputspec':asdict(spec)})
    require(tuple(snapshots)==CAPSULES,'all four manufactured capsules required')
    return snapshots,proof

def execute(v,output,repository_root=None):
    started=time.monotonic();root=v['root'];manufactured=v['manufactured'];output=Path(output).absolute()
    require_worker_limits()
    require(root not in output.parents and root!=output,'output outside frozen source required')
    output=real_directory(output)
    require(sorted(p.name for p in output.iterdir())==['child.log'] and
            (output/'child.log').is_file() and not (output/'child.log').is_symlink(),
            'only custodian-created child.log allowed in fresh output')
    document=load(safe_file(root,'provenance_decoder/INPUT_SPEC.json').read_bytes());metadata_receipt(document)
    require(document['epsilon_exact_ratio']==str(EPSILON) and document['pi_exact_ratio']==str(PI),'represented constants differ')
    check=lambda:resource_guard(started,output)
    source_count=decode_count=0
    with (output/'SOURCE_ATTEMPTS.jsonl').open('xb') as source_journal,(output/'DECODE_ATTEMPTS.jsonl').open('xb') as decode_journal,ExitStack() as stack:
        if manufactured:
            require(repository_root is None,'manufactured branch cannot accept repository root')
            snapshots,proof=manufacture_inputs(document,stack,output)
        else:snapshots,proof=authenticate_inputs(v,repository_root,document,stack,output)
        def authorize_source(kind,payload):
            nonlocal source_count
            require(not manufactured and v['go']['source_go'] is True,'actual source GO absent')
            expected_source=('positive_B','signed_uB')[source_count//64]
            expected_center=Q(-9,2)+Q(2*(source_count%64)+1,128)
            require(source_count<128 and kind=='later_source_taylor' and payload=={'source':expected_source,
                'center':canonical(expected_center),'scope':'LATER_ENDPOINT_SOURCE64CELLS_DEGREE24'},'source callback plan differs')
            source_journal.write(json_bytes({'attempt':source_count,**payload}));source_journal.flush();source_count+=1
        models,certificate=build_models(None if manufactured else authorize_source,manufactured,check)
        require(manufactured or source_count==128,'all128 registered physical source callbacks required')
        write_json(output/'SOURCE_CERTIFICATE.json',certificate)
        targets={source:{endpoint:EndpointFlow(rows,endpoint,check) for endpoint in ENDPOINTS} for source,rows in models.items()}
        write_json(output/'TARGET_BUDGETS.json',{s:{str(t):flow.budget() for t,flow in d.items()} for s,d in targets.items()})
        nodes=0;comparisons=0;scalar_count=0;negative_zero_count=0;rows=[];maximum_radius={'U':Q(0),'W':Q(0)}
        with (output/'NODE_CERTIFICATES.jsonl.gz').open('xb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as compressed:
            for capsule in CAPSULES:
                check();source,grid=capsule.split('/');n=PREFIXES[grid][256];verified=snapshots[capsule]
                def selected(name):
                    nonlocal decode_count
                    require(name in SELECTED_MEMBERS,'unregistered selected member')
                    if not manufactured:
                        require(v['go']['decode_go'] is True,'actual decode GO absent')
                        decode_journal.write(json_bytes({'attempt':decode_count,'capsule':capsule,'member':name}));decode_journal.flush()
                    decode_count+=1
                    return verified.iter_values(name)
                times=tuple(selected('observation_eta.npy'))
                require(len(times)==6 and all(type(x) is codec.ExactReal for x in times) and
                    tuple(x.value for x in times)==OBSERVATIONS,'six exact observation times differ')
                scalar_count+=6;negative_zero_count+=sum(x.negative_zero for x in times)
                streams=[selected(name) for name in ('k.npy','momentum_weights.npy','u_1.npy','w_1.npy','u_2.npy','w_2.npy','u_3.npy','w_3.npy')]
                sentinel=object();acc={t:Accumulator() for t in ENDPOINTS};bycount={c:k for k,c in PREFIXES[grid].items()}
                for index,values in enumerate(zip_longest(*streams,fillvalue=sentinel)):
                    require(index<n and all(x is not sentinel for x in values),'full vectors misaligned')
                    kval,weight,*modes=values
                    require(type(kval) is codec.ExactReal and type(weight) is codec.ExactReal and all(type(x) is codec.ExactComplex for x in modes),'decoded scalar types differ')
                    k,w=kval.value,weight.value
                    require(0<k<256 and w>0,'nonpositive/out-of-band momentum or weight')
                    for K,count in PREFIXES[grid].items():
                        require((index<count and k<K) or (index>=count and k>K),'exact cutoff mask/prefix differs')
                    normalized=[(x.real.value/EPSILON,x.imag.value/EPSILON) for x in modes]
                    flags=[kval.negative_zero,weight.negative_zero]+[v.negative_zero for x in modes for v in (x.real,x.imag)]
                    scalar_count+=14;negative_zero_count+=sum(flags)
                    entry={'capsule':capsule,'index':index,'k':k,'weight':w,
                           'incoming_U':normalized[0],'incoming_W':normalized[1],
                           'negative_zero_flags':flags,'endpoints':{}}
                    for j,endpoint in enumerate(ENDPOINTS):
                        target=targets[source][endpoint].evaluate(k,normalized[0],normalized[1])
                        saved_u,saved_w=normalized[2+2*j],normalized[3+2*j]
                        terms=node_terms(k,w,endpoint,saved_u,saved_w,target,normalized[0],normalized[1]);acc[endpoint].add(k,w,terms)
                        for part in maximum_radius:
                            radius=Q(sum(target[part][a][1]-target[part][a][0] for a in ('real','imag')),2*SCALE)
                            maximum_radius[part]=max(maximum_radius[part],radius)
                        entry['endpoints'][str(endpoint)]={'saved_U':saved_u,'saved_W':saved_w,'target_dyadic96':target,
                                                         'delta_U':terms['delta_U'],'delta_W':terms['delta_W'],
                                                         'canonical_drift':terms['canonical_drift']}
                        comparisons+=1
                    compressed.write(json_bytes(entry));nodes+=1
                    if index+1 in bycount:
                        rows.extend(acc[t].snapshot(capsule,bycount[index+1],t) for t in ENDPOINTS)
                    if index%256==0:check()
                require(all(a.count==n for a in acc.values()),'complete capsule node count differs')
        require(nodes==49152 and comparisons==98304 and decode_count==36 and scalar_count==688152 and len(rows)==24,
                'fixed all-node counters differ')
        write_json(output/'DATA.json',{'scope':'SAVED_ENDPOINT_MINUS_EXACT_FLOW_FROM_EXACT_REPRESENTED_INCOMING_STATE',
                    'manufactured':manufactured,'nodes':nodes,'endpoint_comparisons':comparisons,'rows':rows})
        write_json(output/'SCIENCE_SUMMARY.json',{'status':'MANUFACTURED_PIPELINE_COMPLETE' if manufactured else 'FINITE_RETAINED_ENDPOINT_ERRORS_ENCLOSED',
            'maximum_complete_endpoint_export_L1_radius':maximum_radius,'nodes':nodes,'endpoint_comparisons':comparisons,
            'selected_arrays':decode_count,'decoded_real_slots':scalar_count,'negative_zero_flags':negative_zero_count,
            'source_callbacks':source_count,'incoming_error_counted':False,
            'prior_metric_calibration':'FAIL_UNCHANGED','original_unsaved_continuous_solver_trajectory':'UNRECOVERABLE_FROM_RETAINED_MODE_SNAPSHOTS',
            'continuous_momentum_pressure_contact_time_UV_certificate':'UNRESOLVED',
            'interpretation':'Endpoint discrepancy encloses total subsequent integration/arithmetic error; no separation of unsaved historical per-step components.'})
        write_json(output/'INPUT_AUTHENTICATION.json',{'manufactured':manufactured,'capsules':proof})
        from validate_outputs import validate
        write_json(output/'OUTPUT_READBACK.json',validate(output))
        reauthenticate(v);write_json(output/'RESOURCE.json',check())
        files=('SOURCE_CERTIFICATE.json','TARGET_BUDGETS.json','NODE_CERTIFICATES.jsonl.gz','DATA.json',
               'SCIENCE_SUMMARY.json','RESOURCE.json','SOURCE_ATTEMPTS.jsonl','DECODE_ATTEMPTS.jsonl',
               'INPUT_AUTHENTICATION.json','OUTPUT_READBACK.json')
        output_pins={}
        for name in files:
            raw=(output/name).read_bytes();output_pins[name]={'bytes':len(raw),'sha256':sha(raw)}
        write_json(output/'ENTRY_RECEIPT.json',{'status':'PASS_MANUFACTURED_ENDPOINT_ENTRY' if manufactured else 'PASS_REGISTERED_ENDPOINT_ENTRY',
            'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','manufactured':manufactured,
            'registration_sha256':v['registration_sha256'],'go_sha256':v['go_sha256'],
            'kernel_denied_syscalls':v['kernel_denied_syscalls'],'nodes':nodes,'endpoint_comparisons':comparisons,
            'outputs':output_pins})
