"""All-node worker, imported only after outer registration authentication.

Actual decode and source callbacks have separate authorization journals. The
fabricated path manufactures its own bytes and never opens a repository input.
"""
from contextlib import ExitStack
from dataclasses import asdict
from fractions import Fraction as Q
from pathlib import Path
from itertools import zip_longest
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import resource
import struct
import sys
import time
import zipfile
from registration_guard import (CAPSULES,SOURCES,SourceAuthorization,DecodeAuthorization,
                               require,safe_file,real_directory,load,sha,reauthenticate)
from source_certificate import import_prior,build_models,canonical

SELECTED = ('k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy')
OBSERVATIONS = (Q(-11,2),Q(-9,2),Q(-4),Q(-7,2),Q(-5,2),Q(-3,2))
LIMIT_VALUES = (8<<20,1<<20,8<<20,19,4096,3,65536,64,65536,65536)

def exact_json(value):
    return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()

def write_json(path,value):
    with Path(path).open('xb') as out:
        out.write(exact_json(value))

def import_registered(root,name,path):
    spec = importlib.util.spec_from_file_location(name,safe_file(root,path))
    require(spec is not None and spec.loader is not None,'Registered module import failed')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

def import_numeric(root):
    # Authenticated source directories only; all candidate Python files were
    # inventoried before reaching this function.
    sys.path.insert(0,str(root/'decoder'))
    codec=import_registered(root,'exact_binary80','decoder/exact_binary80.py')
    fixture=import_registered(root,'registered_fabricated_fixtures','decoder/fabricated_fixtures.py')
    engine=import_registered(root,'registered_entire_target','theory/entire_target_engine.py')
    transport=import_registered(root,'registered_discrete_transport','transport/discrete_transport.py')
    prior=import_prior(root)
    require(engine.BITS == 512 and engine.MOMENTUM_DEGREE == 832 and engine.EXPORT_BITS == 96 and
            engine.EXPORTED_L1_GATE == Q(1,10**18),'Fixed target engine constants differ')
    return codec,fixture,engine,transport,prior

def _stream_hash(path,expected_bytes,expected_sha):
    require(path.stat().st_size == expected_bytes,'Opaque input size differs: '+str(path))
    h=hashlib.sha256(); size=0
    with path.open('rb') as stream:
        while block:=stream.read(1<<20):
            h.update(block);size+=len(block)
            require(size<=expected_bytes,'Opaque input grew')
    require(size==expected_bytes and h.hexdigest()==expected_sha,'Opaque input SHA256 differs: '+str(path))

def _strict_inventory(path,pinned):
    with zipfile.ZipFile(path,'r') as z:
        require(z.comment==b'','Unpinned ZIP archive comment')
        found=[{'bytes':i.file_size,'compressed_bytes':i.compress_size,'compression':i.compress_type,
                'crc32_hex':f'{i.CRC:08x}','create_system':i.create_system,'external_attr':i.external_attr,
                'extra_hex':i.extra.hex(),'flag_bits':i.flag_bits,'header_offset':i.header_offset,
                'name':i.filename,'timestamp':list(i.date_time)} for i in z.infolist()]
        require(found==pinned,'Exact pinned ZIP metadata inventory differs')

def _reader_spec(codec,data):
    members=[]
    for m in data['members']:
        clean={k:v for k,v in m.items() if k!='payload_bytes'}
        clean['shape']=tuple(clean['shape']);clean['npy_version']=tuple(clean['npy_version'])
        members.append(codec.MemberSpec(**clean))
    return codec.InputSpec(**{k:v for k,v in data.items() if k!='members'},members=tuple(members))

def authenticate_actual_inputs(v,repository_root,codec,stack,temporary):
    require(not v['fabricated_only'],'Fabricated entry cannot open retained inputs')
    root=v['root'];repo=real_directory(repository_root)
    require(root!=repo and root not in repo.parents,'Actual repository root required')
    document=load(safe_file(root,'provenance/INPUT_SPEC.json').read_bytes())
    require(document['epsilon_exact_ratio']==v['contract']['represented_epsilon'] and
            document['pi_exact_ratio']==v['contract']['represented_Pi'] and
            document['full_distinct_capsule_node_count']==49152 and
            document['all_case_prefix_node_appearances']==86016,'InputSpec constants/universe differ')
    for field,sizefield in (('producer',None),('manifest',None)):
        path=safe_file(repo,document[field+'_repository_path'])
        _stream_hash(path,path.stat().st_size,document[field+'_sha256'])
    for field in ('layout_evidence','original_zip_inventories'):
        p=document[field]
        _stream_hash(safe_file(root,p['path']),p['bytes'],p['sha256'])
    static=document['prior_static_audit']
    _stream_hash(safe_file(repo,static['repository_path']),static['bytes'],static['sha256'])
    snapshots={}; provenance=[]
    limits=codec.Limits(*LIMIT_VALUES)
    require(len(document['capsules'])==4,'Four explicit capsules required')
    for ident,capsule in zip(CAPSULES,document['capsules']):
        require(ident==capsule['source']+'/'+capsule['grid'] and capsule['selected_decode_members']==list(SELECTED),
                'Pinned capsule identity/selection differs')
        n=8192 if capsule['grid']=='coarse' else 16384
        require(capsule['full_node_count']==n,'Full capsule shape changed')
        spec=_reader_spec(codec,capsule['reader_inputspec'])
        require(len(spec.members)==19 and {m.name for m in spec.members if m.decode}==set(SELECTED),
                'Exactly nineteen authenticated members and five selected arrays required')
        for m in spec.members:
            if m.decode:
                require(m.shape==((6,) if m.name=='observation_eta.npy' else (n,)) and
                        m.descr==('<c32' if m.name in ('u_1.npy','w_1.npy') else '<f16'),
                        'Selected member shape/alignment differs')
        original=capsule['original']
        _stream_hash(safe_file(repo,original['repository_path']),original['bytes'],original['sha256'])
        path=safe_file(repo,capsule['repository_path'])
        _stream_hash(path,spec.file_bytes,spec.file_sha256)
        _strict_inventory(path,capsule['strict_zip_inventory'])
        snapshots[ident]=stack.enter_context(codec.verify_inputs(path,spec,limits,temporary_directory=temporary))
        provenance.append({'capsule':ident,'repository_path':capsule['repository_path'],
                           'file_sha256':spec.file_sha256,'input_spec_sha256':snapshots[ident].input_spec_sha256,
                           'full_member_count':19,'authenticated_members':19,'selected_member_names':list(SELECTED),
                           'original_archive_stream_sha256':original['sha256']})
    # No selected-member iterator has been requested until ALL 76 members have
    # passed complete opaque snapshot/hash/header authentication.
    require(tuple(snapshots)==CAPSULES,'Complete four-capsule opaque snapshot missing')
    reauthenticate(v)
    return document,snapshots,provenance

def _slot_exact(q,fixture,padding=b'\x00'*6,negative_zero=False):
    require(type(q) is Q and q.denominator&(q.denominator-1)==0,'Fabricated binary80 must be exact dyadic')
    if q==0:
        return fixture.slot(0,0,negative_zero,padding)
    a=abs(q);e=a.numerator.bit_length()-1-(a.denominator.bit_length()-1)
    significand=a*(Q(2)**(63-e))
    require(significand.denominator==1 and 1<<63<=significand.numerator<1<<64,'Fabricated value exceeds binary80 precision')
    return fixture.slot(significand.numerator,e+16383,q<0,padding)

def manufacture_inputs(root,output,codec,fixture,stack):
    # Only the registered JSON shape metadata is read. No repository root or
    # historical NPY/NPZ pathname participates in this branch.
    metadata=load(safe_file(root,'provenance/INPUT_SPEC.json').read_bytes())
    manufactured=output/'manufactured-inputs';manufactured.mkdir()
    snapshots={};proof=[]
    for sid,source in enumerate(SOURCES):
        for grid,n,den in (('coarse',8192,64),('fine',16384,128)):
            capsule=source+'/'+grid;padding=bytes([sid+1])*6
            p=next(c for c in metadata['capsules'] if c['source']==source and c['grid']==grid)
            entries=[]
            for old in p['reader_inputspec']['members']:
                name,descr,shape=old['name'],old['descr'],tuple(old['shape'])
                if name=='k.npy':
                    payload=b''.join(_slot_exact(Q(2*i+1,den),fixture,padding) for i in range(n))
                elif name=='momentum_weights.npy':
                    payload=_slot_exact(Q(2,den),fixture,padding)*n
                elif name=='observation_eta.npy':
                    payload=b''.join(_slot_exact(q,fixture,padding) for q in OBSERVATIONS)
                elif name in ('u_1.npy','w_1.npy'):
                    divisor=1<<30 if name=='u_1.npy' else 1<<31
                    payload=b''.join(_slot_exact(Q((sid+1)*(i%17+1),divisor),fixture,padding)+
                                     _slot_exact(Q(0) if i==0 else Q(-(sid+1)*(i%13+1),divisor),fixture,padding,
                                                 negative_zero=(i==0)) for i in range(n))
                else:
                    width={'<f16':16,'<c32':32,'<i8':8,'<U8':32}[descr]
                    payload=b'\x00'*(math.prod(shape)*width)
                data=fixture.npy(payload,descr,shape)
                entries.append((name,data,descr,shape,name in SELECTED))
            # Explicit timestamps make the whole manufactured capsule hash
            # reproducible across normal/O runs launched at different times.
            buffer=io.BytesIO();members=[]
            with zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_STORED) as z:
                for name,data,descr,shape,decode in entries:
                    info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
                    info.compress_type=zipfile.ZIP_STORED
                    z.writestr(info,data)
                    members.append(fixture.member(name,data,descr,shape,decode,zipfile.ZIP_STORED))
            data=buffer.getvalue();spec=fixture.input_spec(data,members,'npz')
            path=manufactured/(source+'_'+grid+'.npz')
            with path.open('xb') as out:out.write(data)
            snapshots[capsule]=stack.enter_context(codec.verify_inputs(path,spec,codec.Limits(*LIMIT_VALUES),temporary_directory=manufactured))
            proof.append({'capsule':capsule,'fabricated':True,'file_sha256':spec.file_sha256,
                          'input_spec_sha256':snapshots[capsule].input_spec_sha256,'authenticated_members':19,
                          'selected_member_names':list(SELECTED),'reader_inputspec':asdict(spec)})
    require(tuple(snapshots)==CAPSULES,'Four manufactured snapshots required')
    return metadata,snapshots,proof

def target_radius(target):
    radii={}
    for part in ('U','W'):
        radii[part]=Q(sum(target[part][axis][1]-target[part][axis][0] for axis in ('real','imag')),2*(1<<96))
        require(radii[part]<=Q(1,10**18),'Actual complete target export gate failed')
    return radii

def diagnostics(root,targets,prior_route,fabricated):
    data={method:load(safe_file(root,'prior/'+method+'_DATA.json').read_bytes()) for method in ('primary','independent')}
    report=[];pairs=0;components=0
    for source in SOURCES:
        for k in prior_route.MOMENTA:
            t=targets[source].evaluate(k,None if fabricated else source)
            target_radius(t)
            entry={'source':source,'momentum':canonical(k),'target':t,'method_intersections':{}}
            if not fabricated:
                for method,old in data.items():
                    matches=[r for r in old['whole_rows'] if r['source']==source and Q(r['momentum'])==k]
                    require(len(matches)==1,'Prior diagnostic row missing/duplicate')
                    r=matches[0]; checks={}
                    for part in ('U','W'):
                        for axis in ('real','imag'):
                            lo,hi=(Q(x,1<<96) for x in t[part][axis])
                            priorlo,priorhi=Q(r[part][axis]['lo']),Q(r[part][axis]['hi'])
                            require(max(lo,priorlo)<=min(hi,priorhi),'Prior method target rectangle disjoint')
                            checks[part+'/'+axis]=True;components+=1
                        pairs+=1
                    entry['method_intersections'][method]=checks
            report.append(entry)
    require(len(report)==18,'Exactly eighteen fixed diagnostic target evaluations required')
    # 18 target rows x 2 methods x 2 complex rectangles = 72 intersections;
    # each method separately has 36 rectangle pairs and 72 real components.
    require(fabricated or (pairs==72 and components==144),'Complete two-method diagnostic comparison required')
    return {'rows':report,'target_evaluations':18,'methods':2,'rectangle_pairs_per_method':36,
            'real_component_intersections_per_method':72,'total_rectangle_intersections':pairs,
            'total_real_component_intersections':components,'fabricated':fabricated,
            'certifies_all_nodes_by_interpolation':False}

def execute(v,output,modules,repository_root=None):
    codec,fixture,engine,transport,prior=modules
    started=time.monotonic();root=v['root'];fabricated=v['fabricated_only']
    require(transport.EPSILON==Q(v['contract']['represented_epsilon']) and
            transport.PI==Q(v['contract']['represented_Pi']),'Exact represented normalization constants differ')
    output=real_directory(output)
    require(root!=output and root not in output.parents,'Fresh external output required')
    with ExitStack() as stack:
        if fabricated:
            require(repository_root is None,'Fabricated entry must not accept a repository root')
            document,snapshots,provenance=manufacture_inputs(root,output,codec,fixture,stack)
            for name in ('SOURCE_ATTEMPTS.jsonl','DECODE_ATTEMPTS.jsonl'):(output/name).open('xb').close()
            sourceauth=decodeauth=None
        else:
            require(repository_root is not None,'Actual repository root required')
            document,snapshots,provenance=authenticate_actual_inputs(v,repository_root,codec,stack,output)
            sourceauth=SourceAuthorization(v,output/'SOURCE_ATTEMPTS.jsonl',[r['center'] for r in prior.panels()])
            decodeauth=DecodeAuthorization(v,output/'DECODE_ATTEMPTS.jsonl',snapshots)
        models,certificate=build_models(root,prior,sourceauth,fabricated)
        if not fabricated:require(sourceauth.complete()==22,'Actual physical callback count differs')
        write_json(output/'SOURCE_CERTIFICATE.json',certificate)
        targets={s:engine.EntireTarget(models[s]) for s in SOURCES}
        diagnostic=diagnostics(root,targets,prior,fabricated)
        write_json(output/'DIAGNOSTIC_CROSSCHECK.json',diagnostic)
        coordinates={};summaries=[];maxrad={'U':Q(0),'W':Q(0)}
        arrays=slots=nodes=0;zero_flags=0
        with (output/'NODE_TARGETS.jsonl.gz').open('xb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as compressed:
                for capsule in CAPSULES:
                    source,grid=capsule.split('/');n=8192 if grid=='coarse' else 16384
                    verified=snapshots[capsule]
                    def selected(name):
                        nonlocal arrays
                        if decodeauth:decodeauth.before(capsule,name)
                        arrays+=1
                        return verified.iter_values(name)
                    times=tuple(selected('observation_eta.npy'))
                    require(len(times)==6 and tuple(q.value for q in times)==OBSERVATIONS,'Exact six-entry observation map differs')
                    require(all(type(q) is codec.ExactReal for q in times),'Observation array scalar type differs')
                    slots+=6;acc=transport.PrefixAccumulator();coordinates[capsule]=[]
                    sentinel=object(); plan=transport.PREFIXES[grid]
                    bycount={count:cutoff for cutoff,count in plan.items()}
                    streams=[selected(name) for name in ('k.npy','momentum_weights.npy','u_1.npy','w_1.npy')]
                    for index,values in enumerate(zip_longest(*streams,fillvalue=sentinel)):
                        require(index<n and all(x is not sentinel for x in values),'Full selected vectors misaligned')
                        k,w,u,vv=values
                        require(type(k) is codec.ExactReal and type(w) is codec.ExactReal and
                                type(u) is codec.ExactComplex and type(vv) is codec.ExactComplex,'Selected component types differ')
                        require(0<k.value<256 and w.value>0,'Nonpositive/out-of-band exact k or weight')
                        for cutoff,count in plan.items():
                            require((index<count and k.value<cutoff) or (index>=count and k.value>cutoff),
                                    'Exact cutoff mask differs from fixed positional prefix')
                        coordinates[capsule].append((k.value,w.value))
                        target=targets[source].evaluate(k.value,None if fabricated else source)
                        radii=target_radius(target)
                        for part in maxrad:maxrad[part]=max(maxrad[part],radii[part])
                        uraw=(u.real.value,u.imag.value);wraw=(vv.real.value,vv.imag.value)
                        saved=engine.stored_state_coordinates(k.value,uraw,wraw,transport.EPSILON)
                        error=transport.from_saved_target(k.value,w.value,saved['U_saved'],saved['W_saved'],target)
                        acc.add(error)
                        flags=[q.negative_zero for q in (k,w,u.real,u.imag,vv.real,vv.imag)]
                        zero_flags+=sum(flags)
                        row={'capsule':capsule,'index':index,'k':canonical(k.value),'weight':canonical(w.value),
                             'stored_raw_u':[canonical(q) for q in uraw],'stored_raw_w':[canonical(q) for q in wraw],
                             'negative_zero_flags':flags,'target_dyadic96':target}
                        compressed.write(exact_json(row));nodes+=1;slots+=6
                        if index+1 in bycount:summaries.append(acc.snapshot(capsule,Q(bycount[index+1])))
                    require(acc.count==n,'Incomplete full retained node vector')
        require(nodes==49152 and arrays==20 and slots==294936 and len(summaries)==12,'Fixed all-node counters differ')
        if not fabricated:require(decodeauth.complete()==20,'Actual selected decode attempt count differs')
        equality={}
        for grid in ('coarse','fine'):
            a,b=coordinates['positive_B/'+grid],coordinates['signed_uB/'+grid]
            require(len(a)==len(b),'Cross-source coordinate vectors misaligned')
            equality[grid]={'k_values_equal':all(x[0]==y[0] for x,y in zip(a,b)),
                            'weight_values_equal':all(x[1]==y[1] for x,y in zip(a,b)),
                            'compared_node_pairs':len(a),'basis':'EXACT_DECODED_VALUE_PAIRS_AFTER_SELECTED_DECODE'}
        counts={'physical_source_callbacks':0 if fabricated else 22,'retained_arrays_decoded':0 if fabricated else 20,
                'fabricated_arrays_decoded':20 if fabricated else 0,'selected_scalar_binary80_slots':slots,
                'actual_nodes':0 if fabricated else nodes,'fabricated_nodes':nodes if fabricated else 0,
                'target_node_evaluations':nodes,'diagnostic_target_evaluations':18,
                'total_target_evaluations':nodes+18,'physical_trajectories':0,'likelihood_evaluations':0,
                'unselected_member_decodes':0,'authenticated_capsules':4,'authenticated_members':76}
        scientific={'schema_version':1,'scope':'FABRICATED_ONLY' if fabricated else v['contract']['output_scope'],
                    'counter':counts,'represented_epsilon':canonical(transport.EPSILON),'represented_Pi':canonical(transport.PI),
                    'export_endpoint_exponent':-96,'max_complete_target_export_L1_radius':{k:canonical(q) for k,q in maxrad.items()},
                    'source_coordinate_value_equality':equality,'negative_zero_component_flags':zero_flags,
                    'row_count':12,'distinct_capsule_nodes':49152,'prefix_node_appearances':86016,
                    'rows':transport.encode_exact(summaries),'input_authentication':provenance,
                    'source_budgets':[{'source':s,'source_L1_majorant':canonical(targets[s].source_l1),
                                       'source_error_L1':canonical(targets[s].source_error_l1),
                                       'kernel_tails':{k:canonical(q) for k,q in targets[s].tails.items()}} for s in SOURCES],
                    'state_projection':False,'diagnostic_interpolation_certifies_all_nodes':False,
                    'full_twelve_case_continuous_pressure_contacts':'UNRESOLVED','metric_calibration':'FAIL',
                    'higher_dimensional_Big_Bang_origin':'NOT_ESTABLISHED','external_novelty':'NOT_ASSESSED'}
        write_json(output/'DATA.json',scientific)
        write_json(output/'SCIENCE_SUMMARY.json',scientific)
        write_json(output/'RESOURCE.json',{'wall_seconds':time.monotonic()-started,
                                          'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
        reauthenticate(v)
        return scientific
