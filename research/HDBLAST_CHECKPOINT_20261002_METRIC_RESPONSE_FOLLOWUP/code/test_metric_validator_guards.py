"""Pre-freeze synthetic algebra/format checks; no registered source/mode evaluation."""
import ast
import copy
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
import validate_metric as v
import raw_metric_audit as r

checks=[]
def check(name,condition):
    v.require(condition,name);checks.append(name)
def rejects(name,function):
    try:function()
    except RuntimeError:checks.append(name)
    else:raise RuntimeError('Mutation survived: '+name)

rejects('explicit guards remain enabled',lambda:v.require(False,'synthetic false'))
rejects('scalar NaN',lambda:v.close(float('nan'),0,1,'NaN'))
rejects('scalar infinity',lambda:v.close(float('inf'),0,1,'infinity'))
rejects('negative tolerance',lambda:v.close(0,0,-1,'tolerance'))
rejects('scalar alteration',lambda:v.close(1,2,.1,'changed'))
rejects('raw shape',lambda:r.near_array(np.zeros(2),np.zeros(3),'shape'))
rejects('raw nonfinite',lambda:r.near_array(np.array([float('nan')]),np.zeros(1),'NaN'))
rejects('raw changed integrand',lambda:r.near_array(np.ones(1),np.zeros(1),'changed'))
large=np.array([r.LD('123456789012345678901.2345')],dtype=r.LD)
roundtrip=np.array([float(large[0])],dtype=r.LD)
check('large synthetic binary64 JSON roundtrip',r.near_json(roundtrip,large,'JSON roundtrip')>=0)
rejects('JSON alteration beyond conversion allowance',lambda:r.near_json(roundtrip+1e8,large,'JSON mutation'))
rejects('JSON nonfinite',lambda:r.near_json([float('nan')],[0],'JSON NaN'))
rejects('real dtype downgrade',lambda:r.require_real(np.zeros(2,dtype=np.float64),'dtype'))
rejects('complex dtype downgrade',lambda:r.require_complex(np.zeros(2,dtype=np.complex128),'dtype'))
rejects('odd Simpson endpoint',lambda:r.simpson(np.zeros(4),3,.1))
rejects('out-of-range Simpson endpoint',lambda:r.simpson(np.zeros(4),4,.1))
rejects('zero Simpson endpoint',lambda:r.simpson(np.zeros(4),0,.1))
x=np.arange(9,dtype=r.LD)/8
check('synthetic cubic Simpson exact',abs(r.simpson(x**3,8,r.LD(1)/8)-r.LD(1)/4)<1e-18)
rejects('short public commit',lambda:v.hexadecimal('a'*39,40,'commit'))
rejects('nonhex public commit',lambda:v.hexadecimal('z'*40,40,'commit'))
rejects('nonfinite nested provenance',lambda:v.finite_tree({'raw':[0,{'v':float('nan')}]},'synthetic'))
rejects('nonfinite tail value',lambda:v.finite_tree({'tail':float('inf')},'synthetic'))

# Arbitrary synthetic vectors; no B/uB source, source sample, free-mode
# construction, physical time or mode evolution is evaluated.
vec=np.array([1+2j,-2+1j],dtype=r.CD);vp=np.array([2-1j,3+2j],dtype=r.CD)
dv=np.array([.2+.1j,-.3+.2j],dtype=r.CD);dvp=np.array([-.1+.3j,.4-.2j],dtype=r.CD)
k=np.array([.7,1.3],dtype=r.LD);L,M,h0,h1,epsilon=map(r.LD,('.4','.9','.6','-.3','.2'))
terms=r.operator_terms(vec,vp,dv,dvp,k,L,M,h0,h1,epsilon)
manualD=vp-L*vec;manualDeltaD=dvp-L*dv-epsilon*h1*vec
check('synthetic metric D operator identity',np.max(np.abs(terms['dD']-manualDeltaD))<1e-18)
check('synthetic density operator direct bilinear',np.max(np.abs(terms['deltaR']-(2*np.real(np.conj(manualD)*manualDeltaD)/epsilon+(k*k+M)*2*np.real(np.conj(vec)*dv)/epsilon+2*M*h0*np.real(np.conj(vec)*vec))/2))<2e-18)
check('synthetic opposite mass contacts',np.max(np.abs(terms['massR']+terms['massP']))==0)
rejects('omitted density metric D operator',lambda:r.near_array(terms['deltaR']-terms['operatorR'],terms['deltaR'],'D operator'))
rejects('omitted pressure metric D operator',lambda:r.near_array(terms['deltaP']-terms['operatorP'],terms['deltaP'],'D operator'))
rejects('omitted density metric mass operator',lambda:r.near_array(terms['deltaR']-terms['massR'],terms['deltaR'],'mass'))
rejects('wrong pressure mass sign',lambda:r.near_array(terms['deltaP']-2*terms['massP'],terms['deltaP'],'mass sign'))
baseline=np.array([.1,-.2],dtype=r.LD);sub=np.array([.04,.07],dtype=r.LD)
response=terms['deltaR']-sub-4*h0*baseline
rejects('omitted metric stress prefactor',lambda:r.near_array(terms['deltaR']-sub,response,'prefactor'))
rejects('changed full subtraction',lambda:r.near_array(response+sub,response,'subtraction'))
variance=np.array([.3,.5],dtype=r.LD)-2*h0*baseline
rejects('omitted metric variance prefactor',lambda:r.near_array(variance+2*h0*baseline,variance,'variance prefactor'))
rejects('altered literal baseline',lambda:r.literal_baseline_check(np.ones(2,dtype=r.LD),np.zeros(2,dtype=r.LD),np.ones(2,dtype=r.LD),np.ones(2,dtype=r.LD),'literal'))
check('matching literal baseline allowed',r.literal_baseline_check(np.zeros(2,dtype=r.LD),np.zeros(2,dtype=r.LD),np.ones(2,dtype=r.LD),np.ones(2,dtype=r.LD),'literal')==0)

experiment=json.loads((v.ROOT/'EXPERIMENT.json').read_text())
v.check_experiment(experiment);checks.append('prospective experiment contract')
bad=copy.deepcopy(experiment);bad['gates']['cross_rho']*=10
rejects('changed prospective gate',lambda:v.check_experiment(bad))
bad=copy.deepcopy(experiment);bad['geometry']['geometry_response']=False
rejects('wrong metric channel',lambda:v.check_experiment(bad))
bad=copy.deepcopy(experiment);bad['cutoffs']=[64,128,255]
rejects('changed cutoff',lambda:v.check_experiment(bad))
rejects('duplicate rows',lambda:v.row_index({'rows':[{'source':'positive_B','eta':-5.5},{'source':'positive_B','eta':-5.5}]},experiment,'primary'))
rejects('changed source rows',lambda:v.row_index({'rows':[{'source':'alternative','eta':-5.5}]},experiment,'primary'))
rejects('missing observation rows',lambda:v.row_index({'rows':[]},experiment,'primary'))
records={};v.mutation(records,'tiny sensitivity',1e-14,2e-7,'synthetic operator')
check('small mutation not claimed above operational gate',not records['tiny sensitivity']['operational_gate_exceeded_somewhere'])
v.guard_control(records,'failed state',lambda:v.require(False,'state'))
check('synthetic state mutation honestly classified',records['failed state']['rejected'] and 'no altered-state physical run' in records['failed state']['kind'])

with tempfile.TemporaryDirectory(prefix='metric-validator-guards-',dir='/tmp') as temp:
    root=Path(temp);rejects('existing output directory',lambda:v.validate_output_path(root))
    file=root/'existing';file.write_text('{}');rejects('existing output file',lambda:v.validate_output_path(file))
    rejects('in-checkpoint output',lambda:v.validate_output_path(v.ROOT/'not-created-guard'))
    check('fresh external output',v.validate_output_path(root/'fresh')==root/'fresh')
    rejects('archive path traversal',lambda:r.archive_path(root,{'path':'../outside','sha256':'0'*64}))
    rejects('archive absolute path',lambda:r.archive_path(root,{'path':str(file),'sha256':'0'*64}))
    rejects('archive content hash',lambda:r.archive_path(root,{'path':'existing','sha256':'0'*64}))
    check('archive actual hash accepted',r.archive_path(root,{'path':'existing','sha256':r.sha(file)})==file)
    # Synthetic source files and manifests exercise actual provenance guards
    # without importing/evaluating a producer or defining a quantum state.
    previous=v.ROOT;v.ROOT=root
    try:
        (root/'code').mkdir();(root/'independent').mkdir()
        for name in ['metric_primary.py','source_jet.py','metric_contact_coefficients.py']:(root/'code'/name).write_text('synthetic source artifact\n')
        (root/'EXPERIMENT.json').write_text(json.dumps(experiment))
        entries=[]
        for name in ['forced_metric.py','metric_wkb.py','stable_baselines.py']:
            f=root/'independent'/name;f.write_text('synthetic independent source artifact\n');entries.append({'path':name,'sha256':v.sha(f)})
        manifest=root/'independent/MANIFEST.json';manifest.write_text(json.dumps({'files':entries}))
        pin='a'*64;freeze='b'*40
        pp={'status':'completed','schema_version':1,'route':'scalar_canonical_memory_and_independently_integrated_rational_full_metric_contacts','normalization':experiment['normalization'],'epsilon':1e-4,'mass_law_b':1,'elapsed_seconds':1,'provenance':{'registration_sha256':pin,'public_freeze_commit':freeze,'producer_sha256':v.sha(root/'code/metric_primary.py'),'source_jet_sha256':v.sha(root/'code/source_jet.py'),'contact_coefficients_sha256':v.sha(root/'code/metric_contact_coefficients.py'),'experiment_sha256':v.sha(root/'EXPERIMENT.json'),'dependencies':{'numpy':'2.2.6','scipy':'1.15.3','mpmath':'1.3.0','sympy':'1.14.0'},'python':'3.12.99'}}
        v.primary_provenance(pp,pin,freeze,experiment);checks.append('synthetic primary provenance accepted')
        bad=copy.deepcopy(pp);bad['provenance']['registration_sha256']='c'*64
        rejects('changed registration provenance',lambda:v.primary_provenance(bad,pin,freeze,experiment))
        bad=copy.deepcopy(pp);bad['provenance']['dependencies']['numpy']='2.3.5'
        rejects('wrong physical NumPy provenance',lambda:v.primary_provenance(bad,pin,freeze,experiment))
        bad=copy.deepcopy(pp);bad['elapsed_seconds']=901
        rejects('exceeded primary resources',lambda:v.primary_provenance(bad,pin,freeze,experiment))
        normal={'q':'a0^2 deltaQ/epsilon','q_prime':"(a0^2 deltaQ)'/epsilon",'q_second':"(a0^2 deltaQ)''/epsilon",'rho':'a0^4 delta_rho/epsilon','p':'a0^4 delta_p/epsilon','Q0':'physical Q0,K','rho0':'physical rho0,K','p0':'physical p0,K','current':'a0^2 delta_j/epsilon=q, fixed phi, b=1, x_phi=2'}
        reg={'frozen_configuration':experiment,'independent_manifest_sha256':v.sha(manifest)}
        modes={'status':'passed_internal_gates','schema_version':1,'route':'independent_metric_forced_modes_direct_minimal_stress','gate_failures':[],'freeze_commit':freeze,'epsilon':1e-4,'mass_law_b':1,'configuration':copy.deepcopy(experiment['independent_configuration']),'provenance':{'registration_sha256':pin,'manifest_sha256':v.sha(manifest),'files':entries},'runtime':{'python':'3.12.99','numpy':'2.2.6','python_optimization':0,'real_dtype':np.dtype(r.LD).name,'complex_dtype':np.dtype(r.CD).name,'longdouble_nmant':np.finfo(r.LD).nmant,'longdouble_eps':str(np.finfo(r.LD).eps)},'normalizations':normal,'resources':{'elapsed_seconds':1,'peak_rss_kib':1024}}
        v.independent_provenance(modes,pin,freeze,experiment,reg);checks.append('synthetic independent provenance accepted')
        for name,change in [('changed state',lambda d:d['configuration'].__setitem__('state','alternative')),('Wronskian projection',lambda d:d['configuration'].__setitem__('wronskian_projection',True)),('wrong independent NumPy',lambda d:d['runtime'].__setitem__('numpy','2.3.5')),('optimized physical producer',lambda d:d['runtime'].__setitem__('python_optimization',1)),('exceeded independent wall budget',lambda d:d['resources'].__setitem__('elapsed_seconds',901)),('exceeded memory budget',lambda d:d['resources'].__setitem__('peak_rss_kib',262145)),('wrong manifest provenance',lambda d:d['provenance'].__setitem__('manifest_sha256','d'*64))]:
            bad=copy.deepcopy(modes);change(bad);rejects(name,lambda bad=bad:v.independent_provenance(bad,pin,freeze,experiment,reg))
    finally:v.ROOT=previous
for name in ['validate_metric.py','raw_metric_audit.py']:
    tree=ast.parse((Path(__file__).parent/name).read_text())
    check('no assert statements '+name,not any(isinstance(n,ast.Assert) for n in ast.walk(tree)))
    check('no physical producer imports '+name,not any(isinstance(n,ast.ImportFrom) and n.module in ['metric_primary','forced_metric','source_jet','stable_baselines'] for n in ast.walk(tree)))
print(json.dumps({'status':'PASS','checks':len(checks),'check_names':checks,'python_optimization':sys.flags.optimize,'physical_evaluations':0,'scope':'Synthetic vectors, cubic ledger, malformed metadata/arrays, source pins and static checks only.'},indent=2))
