#!/usr/bin/env python3
"""Prospective homogeneous metric response: scalar memory plus exact contacts.

Importing evaluates no sources, modes, contacts or physical integrals.
All physical evaluation requires the pinned prospective registration first.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings
import mpmath
import numpy
import scipy
import sympy
from scipy.integrate import IntegrationWarning, quad
from source_jet import source_jet_over_epsilon
from metric_contact_coefficients import local_coefficients

ROOT=Path(__file__).resolve().parents[1]
PREF=1/(8*math.pi**2)
DEADLINE=None
AUTHORIZED=False
QUANTITIES=('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')


def require(condition,message):
    if not condition: raise RuntimeError(message)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')


def physical_guard():
    require(AUTHORIZED,'Physical evaluation requires verified prospective registration')
    if DEADLINE is not None:
        require(time.monotonic()<=DEADLINE,'Frozen primary wall-time budget exceeded')


def verify_registration(pin):
    path=ROOT/'FULL_REGISTRATION.json'
    require(sha(path)==pin,'Prospective registration SHA256 mismatch')
    registration=json.loads(path.read_text())
    for name,digest in registration['files'].items():
        target=(ROOT/name).resolve()
        require(target.is_relative_to(ROOT),'Registered path escapes checkpoint')
        require(sha(target)==digest,'Frozen file changed: '+name)
    required={'code/metric_primary.py','code/source_jet.py','code/metric_contact_coefficients.py','code/derive_metric_contacts.py','EXPERIMENT.json'}
    require(required<=set(registration['files']),'Required primary input missing from manifest')
    experiment=json.loads((ROOT/'EXPERIMENT.json').read_text())
    require(registration['frozen_configuration']==experiment,'Frozen configuration differs from experiment')
    require(registration['frozen_gates']==experiment['gates'],'Frozen gates differ from experiment')
    require(registration['independent_manifest_sha256']==sha(ROOT/'independent/MANIFEST.json'),'Independent manifest pin differs')
    return registration,experiment


def metric_jet(source,eta):
    physical_guard()
    return source_jet_over_epsilon(source,eta)


def canonical_forcing_jet(jet,eta):
    """Pure arithmetic g^(n), n=0..3, with exact source h0..h5."""
    L=-1/eta
    return [4*sum(math.comb(n,j)*math.factorial(n-j+1)*L**(n-j+2)*jet[j] for j in range(n+1))
            -2*sum(math.comb(n,j)*math.factorial(n-j)*L**(n-j+1)*jet[j+1] for j in range(n+1))-jet[n+2]
            for n in range(4)]


def forcing_jet(source,eta):
    return canonical_forcing_jet(metric_jet(source,eta),eta)


def integrate(function,left,right,config,**kwargs):
    physical_guard()
    if right<=left: return 0.0,0.0
    with warnings.catch_warnings():
        warnings.simplefilter('error',IntegrationWarning)
        value,error=quad(function,left,right,epsabs=config['epsabs'],epsrel=config['epsrel'],limit=config['limit'],**kwargs)
    physical_guard()
    require(math.isfinite(value) and math.isfinite(error),'Nonfinite primary quadrature')
    return float(value),float(error)


def log_history(source,eta,order,jet,config):
    physical_guard()
    require(order in (1,2,3),'Unregistered derivative order')
    if eta<=-5: return 0.0,0.0
    right=min(eta,-3.0)
    derivative=lambda t:forcing_jet(source,t)[order]
    if right==eta:
        value,error=integrate(derivative,-5.0,right,config,weight='alg-logb',wvar=(0.0,0.0))
    else:
        value,error=integrate(lambda t:derivative(t)*math.log(eta-t),-5.0,right,config)
    constant=math.log(math.sqrt(2)/(-eta))+float(numpy.euler_gamma)+1
    return value+constant*jet[order-1],error


def continuum_q_jet(source,eta,jet,config):
    physical_guard()
    memories=[log_history(source,eta,n,jet,config) for n in (1,2,3)]
    L=-1/eta
    values=[-PREF*memories[0][0],-PREF*(memories[1][0]+L*jet[0]),-PREF*(memories[2][0]+2*L*jet[1]+L*L*jet[0])]
    return values,[PREF*m[1] for m in memories],memories


def finite_history(source,eta,cutoff,order,config):
    physical_guard()
    if eta<=-5:return 0.0,0.0,0
    left,right=max(0.0,eta+3.0),eta+5.0
    count=max(1,math.ceil(2*cutoff*(right-left)/math.pi))
    panel_config=dict(config,epsabs=config['epsabs']/count)
    def function(tau):
        if tau==0:return 0.0
        return forcing_jet(source,eta-tau)[order]*2*math.sin(cutoff*tau)**2/tau
    parts=[integrate(function,left+(right-left)*j/count,left+(right-left)*(j+1)/count,panel_config) for j in range(count)]
    return math.fsum(p[0] for p in parts),math.fsum(p[1] for p in parts),count


def finite_q_jet(source,eta,cutoff,jet,config):
    physical_guard()
    L=-1/eta
    v=cutoff/math.hypot(cutoff,math.sqrt(2)*L)
    A=math.asinh(cutoff/(math.sqrt(2)*L))-v
    Ap=-L*v**3
    App=L*L*v**3*(2-3*v*v)
    memories=[finite_history(source,eta,cutoff,n,config) for n in range(3)]
    values=[PREF*(jet[0]*A-memories[0][0]),PREF*(jet[1]*A+jet[0]*Ap-memories[1][0]),PREF*(jet[2]*A+2*jet[1]*Ap+jet[0]*App-memories[2][0])]
    return values,[PREF*m[1] for m in memories],memories


def direct_response(eta,hjet,gjet,qjet,qerrors,cutoff=None):
    """Two direct subtraction reductions; neither stress is defined by Ward.

    Full-metric bridge and baselines use separately derived rational primitives.
    """
    physical_guard()
    L=a=-1/eta
    v=1.0 if cutoff is None else cutoff/math.hypot(cutoff,math.sqrt(2)*L)
    coefficients=local_coefficients(L,v,hjet)
    Q0=coefficients['Q0']
    Q0p=-L*v*(1-v*v)*(v*(2+v)/(2*(1+v)**2)-v*v/16-5*v**4/16)/(2*math.pi**2)
    J5=v**3/6
    J7=(v**3/3-v**5/5)/4
    J9=(v**3/3-2*v**5/5+v**7/7)/8
    g,gp,gpp=gjet[:3]
    q,qp,qpp=qjet
    density_contact=(3*L*L*g-L*gp)*v**3/(96*math.pi**2)
    pressure_contact=(70*L*L*g*J9-(30*L*L*g+10*L*gp)*J7+(gpp-L*gp-9*L*L*g)*J5)/(48*math.pi**2)
    scalar_rho=(3*L*L*q-L*qp)/2+a*a*Q0*g/2+density_contact
    scalar_p=(qpp-3*L*qp-3*L*L*q)/6-a*a*Q0*g/6+pressure_contact
    final_q=q+coefficients['q_local']
    values={'q':final_q,'q_prime':qp+coefficients['q_local_prime'],'q_second':qpp+coefficients['q_local_second'],'rho':scalar_rho+coefficients['rho_local'],'p':scalar_p+coefficients['p_local'],'Q0':Q0,'rho0':coefficients['rho0'],'p0':coefficients['p0'],'current':final_q}
    e0,e1,e2=qerrors
    estimates={'q':e0,'q_prime':e1,'q_second':e2,'rho':(3*L*L*e0+L*e1)/2,'p':(e2+3*L*e1+3*L*L*e0)/6,'Q0':0.0,'rho0':0.0,'p0':0.0,'current':e0}
    density_contact_prime=((6*L**3*g+2*L*L*gp-L*gpp)*v**3-3*L*(3*L*L*g-L*gp)*v**3*(1-v*v))/(96*math.pi**2)
    scalar_rho_prime=3*L**3*q+L*L*qp-L*qpp/2+a*a*(L*Q0*g+Q0p*g/2+Q0*gp/2)+density_contact_prime
    density_prime=scalar_rho_prime+coefficients['rho_local_prime']
    require(all(math.isfinite(x) for x in values.values()),'Nonfinite primary metric response')
    return {'values':values,'quadrature_error_estimates':estimates,'local_coefficients':dict(coefficients,v=v),'components':{'scalar_memory_q':q,'scalar_memory_q_prime':qp,'scalar_memory_q_second':qpp,'scalar_mass_density':scalar_rho,'scalar_mass_pressure':scalar_p,'full_metric_local_variance':coefficients['q_local'],'full_metric_local_density':coefficients['rho_local'],'full_metric_local_pressure':coefficients['p_local']},'density_derivative':{'scaled_density_prime_over_epsilon':density_prime,'quadrature_error_estimate':3*L**3*e0+L*L*e1+L*e2/2},'trace_from_direct_stresses':-values['rho']+3*values['p']}


def primary_settings(experiment):
    geometry=experiment['geometry']
    require(geometry['H']==1 and geometry['r']==2 and geometry['xi']==0,'Producer requires H1,r2,xi0')
    require(geometry['x']==2 and geometry['incoming_state']=='fixed BD','Fixed physical mass/incoming state required')
    require(all(geometry[key] is True for key in ('physical_mass_squared_fixed','finite_reference_fixed','physical_phi_fixed')),'Metric/reference/scalar variation contract differs')
    require(geometry['geometry_response'] is True,'Metric response must be enabled')
    require(experiment['initial_eta']==-6 and experiment['final_eta']==-1.5,'Unregistered domain')
    require(experiment['source']['ids']==['positive_B','signed_uB'] and experiment['source']['support']==[-5,-3],'Unregistered metric source')
    require(experiment['source']['epsilon']==1e-4 and experiment['source']['mass_law_translation_b']==1,'Unregistered epsilon or current convention')
    require(experiment['source']['physical_scalar_perturbation'] is False,'Physical scalar must remain fixed')
    require(experiment['observation_eta']==[-5.5,-4.5,-4,-3.5,-2.5,-1.5] and experiment['cutoffs']==[64,128,256],'Unregistered observations')
    require(experiment['primary']['producer_wall_seconds']==900,'Unregistered wall budget')
    require(experiment['primary']['epsabs']==1e-12 and experiment['primary']['epsrel']==1e-12 and experiment['primary']['limit']==300,'Unregistered quadrature settings')
    require(experiment['primary']['finite_integral_panel_width_max']=='pi/(2K)' and experiment['primary']['warning_is_failure'] is True,'Unregistered panel or warning policy')
    require(set(experiment['normalization'])==set(QUANTITIES),'Quantity normalization contract differs')
    return experiment


def run(experiment,progress):
    physical_guard()
    rows=[]
    config=experiment['primary']
    for source in experiment['source']['ids']:
        for eta in experiment['observation_eta']:
            physical_guard()
            progress({'rows':rows,'active':{'source':source,'eta':eta}})
            hjet=metric_jet(source,eta)
            gjet=canonical_forcing_jet(hjet,eta)
            q,errors,memories=continuum_q_jet(source,eta,gjet,config)
            continuum=direct_response(eta,hjet,gjet,q,errors)
            continuum['log_history_integrals_over_epsilon']=[{'forcing_derivative_order':n+1,'value':m[0],'quadrature_error_estimate':m[1]} for n,m in enumerate(memories)]
            row={'source':source,'eta':eta,'a':-1/eta,'source_jet_over_epsilon':{f'h{n}':x for n,x in enumerate(hjet)},'canonical_forcing_jet_over_epsilon':{f'g{n}':x for n,x in enumerate(gjet)},'continuum':continuum,'finite_k':[]}
            for K in experiment['cutoffs']:
                q,errors,memories=finite_q_jet(source,eta,K,gjet,config)
                finite=direct_response(eta,hjet,gjet,q,errors,K)
                finite.update(K=K,finite_history_integrals_over_epsilon=[{'forcing_derivative_order':n,'value':m[0],'quadrature_error_estimate':m[1],'time_panels':m[2]} for n,m in enumerate(memories)])
                row['finite_k'].append(finite)
                progress({'rows':rows,'active':row})
            rows.append(row)
            progress({'rows':rows,'active':None})
    return {'schema_version':1,'route':'scalar_canonical_memory_and_independently_integrated_rational_full_metric_contacts','normalization':experiment['normalization'],'epsilon':experiment['source']['epsilon'],'mass_law_b':experiment['source']['mass_law_translation_b'],'rows':rows,'stress_definition':'Density and pressure separately use direct full metric W0/W2/W4 subtraction reductions. Neither is defined by conservation or trace.','uncertainty_scope':'Finite memory quadratures report empirical error estimates. Rational contacts and independent baseline primitives have floating-point arithmetic error. No finite-integral interval certification is claimed.','scope':'Homogeneous linear metric response at fixed physical x=r=2 and fixed incoming BD state. No coupled evolution, heating, stability or shifted-root/state calibration claim.'}


def main():
    global AUTHORIZED,DEADLINE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--public-freeze-commit',required=True)
    args=parser.parse_args()
    registration,experiment=verify_registration(args.registration_sha256)
    require(len(args.public_freeze_commit)==40 and all(c in '0123456789abcdef' for c in args.public_freeze_commit),'Full explicit public freeze commit required')
    primary_settings(experiment)
    versions={'numpy':'2.2.6','scipy':'1.15.3','mpmath':'1.3.0','sympy':'1.14.0'}
    require(platform.python_version_tuple()[:2]==('3','12'),'Frozen Python 3.12 required')
    require(experiment['runtime']['python_major_minor']==[3,12] and experiment['runtime']['physical_producers_optimized'] is False and sys.flags.optimize==0,'Registered physical runtime must be unoptimized Python3.12')
    for name,module in (('numpy',numpy),('scipy',scipy),('mpmath',mpmath),('sympy',sympy)):
        require(module.__version__==versions[name] and experiment['runtime'][name]==versions[name],'Frozen dependency mismatch: '+name)
    output=args.output.resolve()
    require(not output.exists(),'Refusing to overwrite a previous run')
    require(not output.is_relative_to(ROOT),'Fresh output must be outside frozen checkpoint')
    output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    DEADLINE=start+experiment['primary']['producer_wall_seconds']
    AUTHORIZED=True
    provenance={'started_at_utc':datetime.now(timezone.utc).isoformat(),'command':sys.argv,'registration_sha256':args.registration_sha256,'public_freeze_commit':args.public_freeze_commit,'producer_sha256':sha(__file__),'source_jet_sha256':sha(ROOT/'code/source_jet.py'),'contact_coefficients_sha256':sha(ROOT/'code/metric_contact_coefficients.py'),'experiment_sha256':sha(ROOT/'EXPERIMENT.json'),'python':platform.python_version(),'dependencies':{name:module.__version__ for name,module in (('numpy',numpy),('scipy',scipy),('mpmath',mpmath),('sympy',sympy))}}
    write_json(output/'started.json',provenance)
    try:
        def progress(data):write_json(output/'partial_results.json',dict(data,provenance=provenance,elapsed_seconds=time.monotonic()-start,status='running'))
        result=run(experiment,progress)
        physical_guard()
        result.update(provenance=provenance,elapsed_seconds=time.monotonic()-start,status='completed')
        write_json(output/'results.json',result)
        write_json(output/'EXECUTION.json',{'status':'completed','elapsed_seconds':result['elapsed_seconds'],'results_sha256':sha(output/'results.json')})
        print('Primary metric response written to',output,flush=True)
    except BaseException as exc:
        write_json(output/'failure.json',{'exception':repr(exc),'traceback':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start})
        raise


if __name__=='__main__':main()
