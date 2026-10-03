"""Saved-mode reconstruction and separate direct-history metric Ward audit.

No physical producer is imported. The pure Taylor-jet subtraction module is
separately authored; archived subtraction arrays are compared, never used as
definitions. Global all-time Wronskians remain producer declarations; all
archived observation Wronskians and weighted sums are reconstructed.
"""
from fractions import Fraction
import hashlib
import importlib.util
import math
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
NAMES=('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
LD,CD=np.longdouble,np.clongdouble
SUBTRACTION_MODULE=None

def require(condition,message):
    if not condition: raise RuntimeError(message)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def near_array(actual,expected,label):
    actual,expected=np.asarray(actual),np.asarray(expected)
    require(actual.shape==expected.shape,label+': shape')
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)),label+': nonfinite')
    residual=np.max(np.abs(actual-expected),initial=LD(0))
    scale=np.max(np.abs(expected),initial=LD(0))
    allowance=LD('1e-12')+512*np.finfo(LD).eps*max(LD(1),scale)
    require(residual<=allowance,label+': raw mismatch, residual='+str(residual)+' allowance='+str(allowance))
    return float(residual)

def near_json(actual,expected,label):
    """JSON binary64 conversion/arithmetic allowance, separate from raw LD.

    Producers serialize JSON scalars through float(), and primary jet metadata
    also uses short binary64 sums/products. Sixteen binary64 epsilon units
    allow those metadata conversions/operations. This comparator never sets a
    physical cross/refinement/Ward gate and never compares raw LD arrays.
    """
    actual,expected=np.asarray(actual,dtype=LD),np.asarray(expected,dtype=LD)
    require(actual.shape==expected.shape,label+': JSON shape')
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)),label+': JSON nonfinite')
    scale=np.max(np.abs(expected),initial=LD(0));residual=np.max(np.abs(actual-expected),initial=LD(0))
    allowance=LD('1e-12')+16*np.finfo(np.float64).eps*max(LD(1),scale)
    require(residual<=allowance,label+': JSON conversion mismatch')
    return float(residual)

def require_real(array,label):
    require(array.dtype==np.dtype(LD) and np.all(np.isfinite(array)),label+': real extended dtype/nonfinite')

def require_complex(array,label):
    require(array.dtype==np.dtype(CD) and np.all(np.isfinite(array)),label+': complex extended dtype/nonfinite')

def jets(times,source):
    """Exact integer recurrence; no producer/source module import."""
    require(source in ('positive_B','signed_uB'),'Unregistered source id')
    times=np.asarray(times,dtype=LD);z=times+4;inside=np.abs(z)<1
    u,D=z[inside],1-z[inside]**2
    poly=[1];output=np.zeros((6,)+times.shape,dtype=LD)
    for n in range(6):
        value=np.zeros_like(u)
        for coefficient in reversed(poly): value=value*u+coefficient
        output[n][inside]=np.exp(1-1/D)*value/D**(2*n)
        new=[0]*(len(poly)+4)
        for i,coefficient in enumerate(poly):
            if i:
                for j,factor in ((0,1),(2,-2),(4,1)): new[i-1+j]+=i*coefficient*factor
            new[i+1]+=(4*n-2)*coefficient;new[i+3]-=4*n*coefficient
        while len(new)>1 and new[-1]==0: new.pop()
        poly=new
    if source=='signed_uB':
        original=output.copy()
        for n in range(6): output[n]=z*original[n]+(n*original[n-1] if n else 0)
    return output

def forcing(times,hjet):
    L=-LD(1)/np.asarray(times,dtype=LD)
    return np.stack([sum(math.comb(n,j)*(4*math.factorial(j+1)*L**(j+2)*hjet[n-j]
                                       -2*math.factorial(j)*L**(j+1)*hjet[n-j+1])
                         for j in range(n+1))-hjet[n+2] for n in range(4)])

def pure_subtractions(eta,k,hjet):
    global SUBTRACTION_MODULE
    if SUBTRACTION_MODULE is None:
        path=ROOT/'theory/metric_wkb.py'
        require(path.is_file(),'Registered independent pure metric WKB module missing')
        name='registered_metric_audit_wkb'
        spec=importlib.util.spec_from_file_location(name,path)
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
        SUBTRACTION_MODULE=module
    return SUBTRACTION_MODULE.metric_subtractions(eta,k,hjet)

BASELINE_COEFFICIENTS={
 'Q':(16,23,21,15,5),
 'R':(384,493,-204,-1188,-1940,-1990,-868,308,420,105),
 'P':(384,1536,2944,3039,252,-5740,-15260,-19250,-8652,3612,4620,1155)}

def stable_baselines(k,omega,a):
    """Independent power-sum evaluation of exact rationalized differences.

    Polynomial coefficients are derived/checkable from the generic Taylor
    subtraction in verify_raw_metric_baselines.py. The producer uses Horner;
    this audit groups an explicit polynomial sum and never imports it.
    """
    ratio=k/omega;difference=2*a*a/(omega+k)
    polynomial={name:sum(LD(c)*ratio**n for n,c in enumerate(coefficients)) for name,coefficients in BASELINE_COEFFICIENTS.items()}
    return difference**3*polynomial['Q']/(32*k*omega**3),difference**4*polynomial['R']/(1024*k*omega**2),-difference**4*polynomial['P']/(3072*k*omega**2)

def literal_baseline_check(literal,stable,bare,subtraction,label):
    require(literal.shape==stable.shape==bare.shape==subtraction.shape,label+': shape')
    require(all(np.all(np.isfinite(v)) for v in (literal,stable,bare,subtraction)),label+': nonfinite')
    allowance=4096*np.finfo(LD).eps*np.maximum(LD(1),np.abs(bare)+np.abs(subtraction))
    gap=np.abs(literal-stable)
    require(np.all(gap<=allowance),label+': literal subtraction outside leading-term ulp allowance')
    return float(np.max(gap,initial=LD(0)))

def operator_terms(v,vp,dv,dvp,k,L,M,h0,h1,epsilon):
    """Pure bilinear algebra, also testable on arbitrary synthetic vectors."""
    D=vp-L*v;dDmode=dvp-L*dv;dDoperator=-epsilon*h1*v;dD=dDmode+dDoperator
    B=np.real(np.conj(v)*v);B_D=np.real(np.conj(D)*D)
    dB=2*np.real(np.conj(v)*dv)/epsilon
    dDmode2=2*np.real(np.conj(D)*dDmode)/epsilon
    dDoperator2=2*np.real(np.conj(D)*dDoperator)/epsilon
    baseR=(B_D+(k*k+M)*B)/2;baseP=(B_D-(k*k/3+M)*B)/2
    massR=M*h0*B;massP=-massR
    modeR=(dDmode2+(k*k+M)*dB)/2;modeP=(dDmode2-(k*k/3+M)*dB)/2
    operatorR=dDoperator2/2;operatorP=operatorR
    return {'D':D,'dDmode':dDmode,'dDoperator':dDoperator,'dD':dD,
            'B':B,'dB':dB,'baseR':baseR,'baseP':baseP,
            'massR':massR,'massP':massP,'modeR':modeR,'modeP':modeP,
            'operatorR':operatorR,'operatorP':operatorP,
            'deltaR':modeR+operatorR+massR,'deltaP':modeP+operatorP+massP}

def baseline_primitives(times,cutoffs):
    """Closed physical finite-band baselines for every saved history time.

    These rational primitives are pure action-matched baseline identities;
    observation weighted sums still use separately reconstructed raw modes.
    """
    a=-LD(1)/np.asarray(times,dtype=LD)[:,None];K=np.asarray(cutoffs,dtype=LD)[None,:]
    v=K/np.sqrt(K*K+2*a*a);pi=np.arccos(LD(-1))
    Q=(v*v/(2*(1+v))-v**3/48-v**5/16)/(2*pi*pi)
    rpoly=175*v**9+350*v**8-740*v**7-1830*v**6-438*v**5+954*v**4+1852*v**3+2750*v**2-545*v-2880
    ppoly=525*v**11+1050*v**10-1960*v**9-4970*v**8-470*v**7+4030*v**6+3072*v**5+2114*v**4+417*v**3-1280*v**2-1920*v-960
    R=-v*v*rpoly/(7680*pi*pi*(v+1)**2);P=v*v*ppoly/(7680*pi*pi*(v+1)**2)
    return np.stack((Q,R,P),axis=2)

def reconstructed(eta,k,u,w,hjet,epsilon):
    """Use actual unconstrained complex variations, with metric prefactors."""
    eta=LD(str(eta));a=L=-LD(1)/eta;M=2*a*a
    sub=pure_subtractions(eta,k,hjet)
    v=np.exp(-CD(1j)*k*eta)/np.sqrt(2*k);vp=-CD(1j)*k*v
    dv=v*u;dvp=v*(w-CD(1j)*k*u)
    terms=operator_terms(v,vp,dv,dvp,k,L,M,hjet[0],hjet[1],epsilon)
    D,dDmode,dDoperator,dD=(terms[n] for n in ('D','dDmode','dDoperator','dD'))
    B,dB=(terms[n] for n in ('B','dB'))
    dBprime=2*np.real(np.conj(vp)*dv+np.conj(v)*dvp)/epsilon
    g=forcing(eta,hjet)[0]
    vpp=-k*k*v;dvpp=-k*k*dv-epsilon*g*v
    dBsecond=2*np.real(np.conj(vpp)*dv+2*np.conj(vp)*dvp+np.conj(v)*dvpp)/epsilon
    baseR,baseP,massR,massP,modeR,modeP,operatorR,operatorP,deltaR,deltaP=(terms[n] for n in ('baseR','baseP','massR','massP','modeR','modeP','operatorR','operatorP','deltaR','deltaP'))
    S=sub['S']['value'];S1=sub['S']['first'];S2=sub['S']['second']
    dS=sub['S']['delta'];dS1=sub['S']['delta_first'];dS2=sub['S']['delta_second']
    Qliteral=B-S;Rliteral=baseR-sub['R']['value'];Pliteral=baseP-sub['P']['value']
    Qc0,R0,P0=stable_baselines(k,np.sqrt(k*k+M),a)
    literal_baseline_check(Qliteral,Qc0,B,S,'Variance baseline literal/stable')
    literal_baseline_check(Rliteral,R0,baseR,sub['R']['value'],'Density baseline literal/stable')
    literal_baseline_check(Pliteral,P0,baseP,sub['P']['value'],'Pressure baseline literal/stable')
    q=dB-dS-2*hjet[0]*Qc0
    q1=dBprime-dS1-2*hjet[1]*Qc0+2*hjet[0]*S1
    q2=dBsecond-dS2-2*hjet[2]*Qc0+4*hjet[1]*S1+2*hjet[0]*S2
    pi=np.arccos(LD(-1));measure=k*k/(2*pi*pi)
    data={'q':measure*q,'q_prime':measure*q1,'q_second':measure*q2,
          'rho':measure*(deltaR-sub['R']['delta']-4*hjet[0]*R0),
          'p':measure*(deltaP-sub['P']['delta']-4*hjet[0]*P0),
          'Q0':measure*Qc0/a**2,'rho0':measure*R0/a**4,'p0':measure*P0/a**4}
    data['current']=data['q']
    linearW=(dv*np.conj(vp)+v*np.conj(dvp)-dvp*np.conj(v)-vp*np.conj(dv))/epsilon
    canonicalEnergy=(2*np.real(np.conj(vp)*dvp)+k*k*(epsilon*dB))/epsilon
    W=float(np.max(np.abs(2*u.real-w.imag/k))/epsilon)
    physical=float(np.max(np.abs(linearW)))
    archive={'v':v,'v_prime':vp,'delta_v':dv,'delta_v_prime':dvp,'D_v':D,'delta_D_v':dD,
             'delta_D_v_mode':dDmode,'delta_D_v_operator':dDoperator,'bare_variance0':B,
             'bare_R0':baseR,'bare_P0':baseP,'bare_delta_R_mode':modeR,'bare_delta_P_mode':modeP,
             'bare_delta_R_operator':operatorR,'bare_delta_P_operator':operatorP,
             'bare_delta_R_mass':massR,'bare_delta_P_mass':massP,
             'bare_delta_R_scale':-4*hjet[0]*baseR,'bare_delta_P_scale':-4*hjet[0]*baseP,
             'bare_delta_Q_scale':-2*hjet[0]*B,
             'sub_delta_R_scale':-4*hjet[0]*sub['R']['value'],
             'sub_delta_P_scale':-4*hjet[0]*sub['P']['value'],
             'sub_delta_Q_scale':-2*hjet[0]*S,
             'bare_delta_R':deltaR,'bare_delta_P':deltaP,
             'bare_delta_A':dB,'bare_delta_A_prime':dBprime,'bare_delta_A_second':dBsecond,
             'canonical_energy_variation':canonicalEnergy,'linear_wronskian':linearW,
             'baseline_Q_combined':Qc0,'baseline_R_combined':R0,'baseline_P_combined':P0,
             'baseline_Q_literal':Qliteral,'baseline_R_literal':Rliteral,'baseline_P_literal':Pliteral,
             'source_jet':hjet.copy(),'forcing_jet':forcing(eta,hjet)}
    # All generic Taylor quantities are reconstructed independently.
    flat={}
    for name in ('R0','R2','R4','P0','P2','P4','S0','S2','R','P','S','W2','W4','J2','J4'):
        fields=sub[name];flat[name]=fields['value'];flat['delta'+name]=fields['delta']
        if name in ('S','W2'):
            flat[name+'_prime']=fields['first'];flat[name+'_second']=fields['second']
            flat['delta'+name+'_prime']=fields['delta_first'];flat['delta'+name+'_second']=fields['delta_second']
    archive.update({'sub_'+key:value for key,value in flat.items()})
    for n in range(3):
        archive[f'C_{n}']=np.asarray(2*math.factorial(n+1)*a**(n+2),dtype=LD)
    archive['delta_C_0']=np.asarray(hjet[2]+2*L*hjet[1],dtype=LD)
    archive['delta_C_1']=np.asarray(hjet[3]+2*L*hjet[2]+2*L*L*hjet[1],dtype=LD)
    archive['delta_C_2']=np.asarray(hjet[4]+2*L*hjet[3]+4*L*L*hjet[2]+4*L**3*hjet[1],dtype=LD)
    # Same saved modes, deliberately altered definitions; no new physical run.
    mutations={'missing_density_metric_D_operator':-measure*operatorR,
               'missing_pressure_metric_D_operator':-measure*operatorP,
               'missing_density_metric_mass_operator':-measure*massR,
               'missing_pressure_metric_mass_operator':-measure*massP,
               'missing_density_metric_scale_contact':4*measure*hjet[0]*R0,
               'missing_pressure_metric_scale_contact':4*measure*hjet[0]*P0,
               'missing_variance_metric_scale_contact':2*measure*hjet[0]*Qc0,
               'missing_density_metric_W4':measure*sub['R4']['delta'],
               'missing_pressure_metric_W4':measure*sub['P4']['delta']}
    return data,archive,W,physical,mutations,sub

def simpson(values,end,dt):
    require(end>0 and end%2==0 and end<len(values),'Non-Simpson observation endpoint')
    return dt/3*(values[0]+values[end]+4*np.sum(values[1:end:2],axis=0,dtype=LD)+2*np.sum(values[2:end:2],axis=0,dtype=LD))

def archive_path(directory,record):
    relative=Path(record['path']);path=(Path(directory)/relative).resolve()
    require(not relative.is_absolute() and '..' not in relative.parts and path.is_relative_to(Path(directory).resolve()),'Archive path escape')
    require(sha(path)==record['sha256'],'Archive hash mismatch')
    return path

def audit_run(run,directory,experiment):
    source,setting=run['source'],run['setting'];config=experiment['independent']['settings'][setting]
    dt=LD(Fraction(config['dt']).numerator)/Fraction(config['dt']).denominator
    width=LD(Fraction(config['momentum_panel_width']).numerator)/Fraction(config['momentum_panel_width']).denominator
    steps=int((LD(str(experiment['final_eta']))-LD(str(experiment['initial_eta'])))/dt)
    order=experiment['independent']['momentum_order'];nodes=int(LD(max(experiment['cutoffs']))/width)*order
    require(order==16 and run['dt']==float(dt) and run['momentum_panel_width']==float(width),'Registered run spacings')
    require(run['time_steps']==steps and run['momentum_nodes']==nodes and run['time_gauss_nodes']==8,'Registered run counts')
    require(run['local_time_momentum_pairs']==steps*8*nodes and run['direct_stress_history_points']==steps+1,'Declared work inventory')
    require(np.finfo(LD).nmant>=63,'Validator extended precision unavailable')
    path=archive_path(directory,run['archive'])
    report={'source':source,'setting':setting,'archive_sha256':sha(path),'reconstructed_observation_points':0,'ward_endpoints':[],'raw_mutations':{}}
    common={'k','momentum_weights','history_eta','history_values','history_ledger_integrand','history_baseline_contact','history_source_jet','history_forcing_jet','history_geometry','history_quantity_names','history_geometry_names','history_cutoffs','observation_eta'}
    text={'history_quantity_names','history_geometry_names','history_cutoffs'}
    complex_stems={'u','w','v','v_prime','delta_v','delta_v_prime','D_v','delta_D_v','delta_D_v_mode','delta_D_v_operator','linear_wronskian'}
    with np.load(path,allow_pickle=False) as arc:
        require(common<=set(arc.files),'Raw history arrays missing')
        require(tuple(arc['history_quantity_names'].tolist())==NAMES and arc['history_cutoffs'].tolist()==experiment['cutoffs'],'Raw quantity names/cutoffs')
        require(tuple(arc['history_geometry_names'].tolist())==('a0','L','L_prime','L_second','L_third','M','M_prime','M_second','M_third','M_fourth'),'Geometry names')
        for name in set(arc.files)-text:
            stem=name.rsplit('_',1)[0] if name.rsplit('_',1)[-1].isdigit() else name
            (require_complex if stem in complex_stems else require_real)(arc[name],name)
        k,wk=arc['k'],arc['momentum_weights'];gx,gw=np.polynomial.legendre.leggauss(16)
        panels=int(LD(max(experiment['cutoffs']))/width)
        expected_k=(np.arange(panels,dtype=LD)[:,None]*width+(gx.astype(LD)[None,:]+1)*width/2).ravel()
        expected_w=np.broadcast_to(gw.astype(LD)[None,:]*width/2,(panels,16)).ravel()
        require(k.shape==wk.shape==(nodes,),'Momentum array shapes')
        near_array(k,expected_k,'Momentum nodes');near_array(wk,expected_w,'Momentum weights')
        time=arc['history_eta'];history=arc['history_values'];hj=arc['history_source_jet'];gj=arc['history_forcing_jet']
        require(time.shape==(steps+1,) and history.shape==(steps+1,3,9) and hj.shape==(steps+1,6) and gj.shape==(steps+1,4),'History shapes')
        near_array(time,LD(str(experiment['initial_eta']))+np.arange(steps+1,dtype=LD)*dt,'History time grid')
        near_array(hj,jets(time,source).T,'Complete source history');near_array(gj,forcing(time,hj.T).T,'Complete forcing history')
        near_array(history[:,:,[5,6,7]],baseline_primitives(time,experiment['cutoffs']),'Complete finite-band baseline histories')
        a=-1/time;M=2*a*a
        geometry=np.stack((a,a,a*a,2*a**3,6*a**4,M,2*M*a,6*M*a*a,24*M*a**3,120*M*a**4),axis=1)
        near_array(arc['history_geometry'],geometry,'Geometry history')
        baseline=3*hj[:,1,None]*a[:,None]**4*(history[:,:,6]+history[:,:,7])
        ledger=a[:,None]*(history[:,:,3]-3*history[:,:,4])-baseline
        near_array(arc['history_baseline_contact'],baseline,'Metric Ward baseline contact')
        near_array(arc['history_ledger_integrand'],ledger,'Direct-history Ward integrand')
        initial=time<=LD(str(experiment['source']['support'][0]));response_indices=[0,1,2,3,4,8]
        require(np.count_nonzero(history[initial][:,:,response_indices])==0,'Initial direct-history response')
        require([r['eta'] for r in run['rows']]==experiment['observation_eta'],'Raw observation grid')
        near_array(arc['observation_eta'],np.array(experiment['observation_eta'],dtype=LD),'Observation eta array')
        expected_names=set(common);epsilon=LD(str(experiment['source']['epsilon']))
        for i,row in enumerate(run['rows']):
            eta=row['eta'];hjet=arc[f'source_jet_{i}'];u,w=arc[f'u_{i}'],arc[f'w_{i}']
            require(row['source']==source and u.shape==w.shape==(nodes,) and hjet.shape==(6,),'Mode/source shape')
            near_array(hjet,jets(np.array([eta],dtype=LD),source)[:,0],'Observation source jets')
            near_json([row['source_jet_over_epsilon'][f'h{n}'] for n in range(6)],hjet,'JSON source jet')
            near_json([row['forcing_jet_over_epsilon'][f'g{n}'] for n in range(4)],forcing(LD(str(eta)),hjet),'JSON forcing jet')
            data,archive,W,physical,mutations,sub=reconstructed(eta,k,u,w,hjet,epsilon)
            require(W<=experiment['gates']['wronskian_over_epsilon'] and physical<=experiment['gates']['wronskian_over_epsilon'],'Saved-mode Wronskians')
            near_json([row['wronskian_scaled'],row['physical_wronskian_scaled']],[W,physical],'JSON Wronskians')
            for name,array in data.items():
                target=f'integrand_{name}_{i}';expected_names.add(target);near_array(arc[target],array,'Saved integrand '+name)
            for name,array in archive.items():
                target=f'{name}_{i}';expected_names.add(target);near_array(arc[target],array,'Saved mode/operator/subtraction '+name)
            for n in range(5):
                expected_names.add(f'omega_{n}_{i}');expected_names.add(f'delta_omega_{n}_{i}')
            # Frequency jets are recomputed from fixed-r geometry, not saved derivatives.
            omega=[np.sqrt(k*k+2*a[int((LD(str(eta))-time[0])/dt)]**2)]
            mass0=2*(-LD(1)/LD(str(eta)))**2;L=-LD(1)/LD(str(eta))
            mj=[math.factorial(n+1)*mass0*L**n for n in range(5)]
            dm=[2*sum(math.comb(n,j)*mj[j]*hjet[n-j] for j in range(n+1)) for n in range(5)]
            domega=[dm[0]/(2*omega[0])]
            for n in range(1,5):
                omega.append((mj[n]-sum(math.comb(n,j)*omega[j]*omega[n-j] for j in range(1,n)))/(2*omega[0]))
                domega.append((dm[n]-sum(math.comb(n,j)*(domega[j]*omega[n-j]+omega[j]*domega[n-j]) for j in range(1,n))-2*domega[0]*omega[n])/(2*omega[0]))
            for n in range(5):
                near_array(arc[f'omega_{n}_{i}'],omega[n],'Saved baseline frequency derivative')
                near_array(arc[f'delta_omega_{n}_{i}'],domega[n],'Saved metric frequency derivative')
            expected_names.update((f'u_{i}',f'w_{i}'))
            step=int((LD(str(eta))-time[0])/dt);integrated=simpson(ledger,step,dt)
            require([p['K'] for p in row['finite_k']]==experiment['cutoffs'],'Raw cutoff grid')
            for j,point in enumerate(row['finite_k']):
                count=int(LD(point['K'])/width)*16
                for n,name in enumerate(NAMES):
                    value=np.sum(data[name][:count]*wk[:count],dtype=LD)
                    near_json([point['values'][name]],[value],'Raw/JSON '+name)
                    near_array(history[step,j,n],value,'Raw/direct-history '+name)
                endpoint=history[step,j,3]-history[0,j,3];residual=float(abs(integrated[j]-endpoint))
                near_json([point['ward']['ledger'],point['ward']['direct_endpoint'],point['ward']['endpoint_difference']],[integrated[j],endpoint,residual],'Recorded Ward ledger')
                report['ward_endpoints'].append({'eta':eta,'K':point['K'],'residual':residual,'ledger':float(integrated[j]),'gate':experiment['gates']['ward_endpoint'] if setting=='fine' else None})
                for name,array in mutations.items():
                    residual=float(abs(np.sum(array[:count]*wk[:count],dtype=LD)))
                    report['raw_mutations'][name]=max(report['raw_mutations'].get(name,0),residual)
                # Actual Simpson witness for omitting the finite-band baseline term.
                wrong=abs(simpson(baseline[:,j],step,dt))
                report['raw_mutations']['missing_Ward_baseline_contact']=max(report['raw_mutations'].get('missing_Ward_baseline_contact',0),float(wrong))
                report['reconstructed_observation_points']+=1
            if eta<=experiment['source']['support'][0]:require(np.count_nonzero(u)==np.count_nonzero(w)==0,'Initial forced modes changed')
        require(set(arc.files)==expected_names,'Unexpected/missing archived arrays')
    return report
