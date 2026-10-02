"""Raw independent stress evidence reconstruction; no producer imports or import-time experiment.

Saved histories contain direct stresses, not full modes at each history time. We
reconstruct every saved observation from modes and independently integrate the
complete direct-history Ward ledger. Global producer Wronskian maxima remain
provenance-backed declarations; observation Wronskians are reconstructed here.
"""
from fractions import Fraction
import hashlib
from pathlib import Path
import numpy as np

NAMES = ('q','q_prime','q_second','rho','p','Q0','anomaly','current')
LD, CD = np.longdouble, np.clongdouble

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def near_array(actual, expected, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, label+': shape')
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), label+': nonfinite')
    residual = np.max(np.abs(actual-expected), initial=LD(0))
    scale = np.max(np.abs(expected), initial=LD(0))
    require(residual <= LD('1e-12')+512*np.finfo(LD).eps*max(LD(1),scale), label+': raw mismatch')
    return float(residual)

def jets(times, source):
    """Independent exact-integer recurrence, then extended arithmetic evaluation."""
    require(source in ('positive_B','signed_uB'), 'Unregistered source id')
    times = np.asarray(times, dtype=LD)
    z = times+4
    inside = np.abs(z)<1
    u, D = z[inside], 1-z[inside]**2
    poly = [1]
    output = np.zeros((6,)+times.shape, dtype=LD)
    for n in range(6):
        val = np.zeros_like(u)
        for c in reversed(poly):
            val = val*u+c
        output[n][inside] = np.exp(1-1/D)*val/D**(2*n)
        # P_(n+1)=(1-u²)² P'_n + [(4n-2)u-4nu³] P_n.
        next_poly = [0]*(len(poly)+4)
        for i,c in enumerate(poly):
            if i:
                for j,d in ((0,1),(2,-2),(4,1)):
                    next_poly[i-1+j] += i*c*d
            next_poly[i+1] += (4*n-2)*c
            next_poly[i+3] -= 4*n*c
        while len(next_poly)>1 and next_poly[-1]==0:
            next_poly.pop()
        poly = next_poly
    if source == 'signed_uB':
        base = output.copy()
        for n in range(6):
            output[n] = z*base[n]+(n*base[n-1] if n else 0)
    return output

def reconstructed(eta, k, u, w, jet, epsilon):
    """Rebuild direct minimal bilinears and complete W2/W4 variations."""
    a = L = -1/LD(str(eta)); m2=2*a*a
    s,s1,s2 = epsilon*jet[:3]
    om=np.sqrt(k*k+m2); om1=L*m2/om
    om2=3*L*L*m2/om-L*L*m2*m2/om**3
    U=-L*L/om-om2/(4*om**2)+3*om1**2/(8*om**3)
    d2=s/(2*om); d21=s1/(2*om)-s*om1/(2*om**2)
    d22=s2/(2*om)-s1*om1/om**2-s*om2/(2*om**2)+s*om1**2/om**3
    d4=-U*d2/om-d22/(4*om**2)+om2*d2/(4*om**3)+3*om1*d21/(4*om**3)-3*om1**2*d2/(4*om**4)
    beta=-k*k/3-m2; c=1-beta/om**2
    j=L*d21/om**2-2*L*om1*d2/om**3+om1*d21/(2*om**3)-3*om1**2*d2/(4*om**4)
    sub_r=(s/om-L*L*d2/om**2+j)/4
    sub_p=(c*d2-s/om+c*d4+2*beta*U*d2/om**3+s*U/om**2-L*L*d2/om**2+j)/4
    v=np.exp(-CD(1j)*k*LD(str(eta)))/np.sqrt(2*k)
    vp=-CD(1j)*k*v; dv=v*u; dvp=v*(w-CD(1j)*k*u)
    A=2*np.real(np.conj(v)*dv)
    Dv=vp-L*v; dDv=dvp-L*dv
    B=2*np.real(np.conj(Dv)*dDv)
    bare_r=(B+(k*k+m2)*A+s/(2*k))/2
    bare_p=(B-(k*k/3+m2)*A-s/(2*k))/2
    A1=w.real/k; A2=np.real(CD(2j)*k*w-s)/k
    C=-s/(4*om**3); C1=-s1/(4*om**3)+3*s*om1/(4*om**4)
    C2=-s2/(4*om**3)+3*s1*om1/(2*om**4)+3*s*om2/(4*om**4)-3*s*om1**2/om**5
    baseline_sub=1/(2*om)-U/(2*om**2)
    an=sub_r-3*sub_p-m2*C-s*baseline_sub+(C2-2*L*C1-2*L*L*C)/2
    measure=k*k/(2*np.arccos(LD(-1))**2)
    data={'q':measure*(A-C)/epsilon,'q_prime':measure*(A1-C1)/epsilon,
          'q_second':measure*(A2-C2)/epsilon,'rho':measure*(bare_r-sub_r)/epsilon,
          'p':measure*(bare_p-sub_p)/epsilon,'Q0':measure*(1/(2*k)-baseline_sub)/(a*a),
          'anomaly':measure*an/epsilon}
    data['current']=data['q']+jet[0]*data['Q0']/4
    physical=np.max(np.abs(dv*np.conj(vp)+v*np.conj(dvp)-dvp*np.conj(v)-vp*np.conj(dv)))/epsilon
    W=np.max(np.abs(2*u.real-w.imag/k))/epsilon
    # Concrete alternative operators, applied to SAME linear modes. No new run.
    controls={'missing_density_mass_operator':-measure*s/(4*k)/epsilon,
              'missing_pressure_mass_operator':measure*s/(4*k)/epsilon,
              'missing_pressure_W4':measure*c*d4/(4*epsilon),
              'improved_stress_density':-measure*(L*L*A-L*A1)/(2*epsilon),
              'improved_stress_pressure':-measure*(A2-3*L*A1+3*L*L*A)/(6*epsilon)}
    archive={'bare_rho':bare_r/epsilon,'sub_rho':sub_r/epsilon,'bare_p':bare_p/epsilon,
             'sub_p':sub_p/epsilon,'delta_W2':d2,'delta_W4':d4}
    return data,archive,float(W),float(physical),controls

def simpson(values, end, dt):
    require(end>0 and end%2==0, 'Non-Simpson observation endpoint')
    return dt/3*(values[0]+values[end]+4*np.sum(values[1:end:2],axis=0,dtype=LD)+2*np.sum(values[2:end:2],axis=0,dtype=LD))

def audit_run(run, directory, experiment):
    source, setting=run['source'],run['setting']
    config=experiment['independent']['settings'][setting]
    dt=LD(Fraction(config['dt']).numerator)/Fraction(config['dt']).denominator
    width=LD(Fraction(config['momentum_panel_width']).numerator)/Fraction(config['momentum_panel_width']).denominator
    steps=int((LD(str(experiment['final_eta']))-LD(str(experiment['initial_eta'])))/dt)
    nodes=int(LD(max(experiment['cutoffs']))/width)*experiment['independent']['momentum_order']
    require(run['dt']==float(dt) and run['momentum_panel_width']==float(width), 'Run grid spacing changed')
    require(run['time_steps']==steps and run['momentum_nodes']==nodes and run['time_gauss_nodes']==8, 'Run grid counts changed')
    require(run['local_time_momentum_pairs']==steps*8*nodes and run['direct_stress_history_points']==steps+1, 'Declared evaluation count changed')
    path=(Path(directory)/run['archive']['path']).resolve()
    require(not Path(run['archive']['path']).is_absolute() and path.is_relative_to(Path(directory).resolve()), 'Archive path escape')
    require(sha(path)==run['archive']['sha256'], 'Archive hash mismatch')
    require(np.finfo(LD).nmant>=63, 'Validator extended precision unavailable')
    report={'source':source,'setting':setting,'archive_sha256':sha(path),'reconstructed_observation_points':0,'ward_endpoints':[],'raw_mutations':{}}
    with np.load(path,allow_pickle=False) as arc:
        expected={'k','momentum_weights','history_eta','history_values','history_ledger_integrand','history_source_jet','history_d_prime','history_quantity_names','history_cutoffs'}
        for i in range(len(experiment['observation_eta'])):
            expected.update({f'{name}_{i}' for name in ('u','w','source_jet','bare_rho','sub_rho','bare_p','sub_p','delta_W2','delta_W4')})
            expected.update({f'integrand_{name}_{i}' for name in NAMES})
        require(set(arc.files)==expected, 'Unexpected or missing raw arrays')
        require(tuple(arc['history_quantity_names'].tolist())==NAMES and arc['history_cutoffs'].tolist()==experiment['cutoffs'], 'Raw names/cutoffs changed')
        for name in expected-{'history_quantity_names','history_cutoffs'}:
            arr=arc[name]
            dtype=CD if name.startswith(('u_','w_')) else LD
            require(arr.dtype==np.dtype(dtype) and np.all(np.isfinite(arr)), 'Raw dtype/nonfinite: '+name)
        k,wk=arc['k'],arc['momentum_weights']
        gx,gw=np.polynomial.legendre.leggauss(16)
        panels=int(LD(max(experiment['cutoffs']))/width)
        ek=(np.arange(panels,dtype=LD)[:,None]*width+(gx.astype(LD)[None,:]+1)*width/2).ravel()
        ew=np.broadcast_to(gw.astype(LD)[None,:]*width/2,(panels,16)).ravel()
        near_array(k,ek,'Momentum nodes');near_array(wk,ew,'Momentum weights')
        require(k.shape==wk.shape==(nodes,), 'Momentum raw shape')
        time=arc['history_eta']; history=arc['history_values']; hj=arc['history_source_jet']
        expected_time=LD(str(experiment['initial_eta']))+np.arange(steps+1,dtype=LD)*dt
        near_array(time,expected_time,'History time grid')
        require(history.shape==(steps+1,3,8) and hj.shape==(steps+1,6), 'History shape')
        near_array(hj,jets(time,source).T,'Complete history source jets')
        a=-1/time; dp=LD(str(experiment['source']['epsilon']))*(hj[:,1]-2*a*hj[:,0])/(a*a)
        near_array(arc['history_d_prime'],dp,'Physical mass derivative history')
        ledger=a[:,None]**2*history[:,:,5]*(hj[:,1,None]-2*a[:,None]*hj[:,0,None])/2+a[:,None]*history[:,:,3]-3*a[:,None]*history[:,:,4]
        near_array(arc['history_ledger_integrand'],ledger,'Direct-history Ward integrand')
        initial=time<=LD(str(experiment['source']['support'][0]))
        require(np.count_nonzero(history[initial][:,:,[0,1,2,3,4,6,7]])==0, 'Causality or initial-state history altered')
        require([r['eta'] for r in run['rows']]==experiment['observation_eta'], 'Raw observation grid')
        eps=LD(str(experiment['source']['epsilon']))
        for i,row in enumerate(run['rows']):
            eta=row['eta']; jet=arc[f'source_jet_{i}']; u,w=arc[f'u_{i}'],arc[f'w_{i}']
            require(row['source']==source and u.shape==w.shape==(nodes,) and jet.shape==(6,), 'Mode/source shape')
            near_array(jet,jets(np.array([eta],dtype=LD),source)[:,0],'Observation source jets')
            near_array(np.array([row['source_jet_over_epsilon'][f'f{n}'] for n in range(6)],dtype=LD),jet,'JSON source jet')
            calculated,archived,W,physical,mutations=reconstructed(eta,k,u,w,jet,eps)
            require(W<=experiment['gates']['wronskian_over_epsilon'] and physical<=experiment['gates']['wronskian_over_epsilon'],'Reconstructed Wronskian')
            near_array(np.array([row['wronskian_scaled'],row['physical_wronskian_scaled']]),np.array([W,physical]),'JSON Wronskian')
            for name,array in calculated.items():near_array(arc[f'integrand_{name}_{i}'],array,'Saved integrand '+name)
            for name,array in archived.items():near_array(arc[f'{name}_{i}'],array,'Saved bare/subtraction '+name)
            step=int((LD(str(eta))-time[0])/dt)
            integrated_ward=simpson(ledger,step,dt)
            for j,point in enumerate(row['finite_k']):
                cutoff=experiment['cutoffs'][j]; require(point['K']==cutoff, 'Raw cutoff')
                count=int(LD(cutoff)/width)*16
                for n,name in enumerate(NAMES):
                    value=np.sum(calculated[name][:count]*wk[:count],dtype=LD)
                    near_array(np.array([point['values'][name],history[step,j,n]],dtype=LD),np.array([value,value],dtype=LD),'Raw/JSON/history '+name)
                endpoint=history[step,j,3]-history[0,j,3]
                residual=float(abs(integrated_ward[j]-endpoint))
                near_array(np.array([point['ward']['ledger'],point['ward']['direct_endpoint'],point['ward']['endpoint_difference']]),np.array([float(integrated_ward[j]),float(endpoint),residual]),'Recorded Ward ledger')
                report['ward_endpoints'].append({'eta':eta,'K':cutoff,'residual':residual,'gate':experiment['gates']['ward_endpoint'] if setting=='fine' else None})
                for name,array in mutations.items():
                    difference=float(abs(np.sum(array[:count]*wk[:count],dtype=LD)))
                    report['raw_mutations'][name]=max(report['raw_mutations'].get(name,0),difference)
                report['reconstructed_observation_points']+=1
            if eta<=experiment['source']['support'][0]:require(np.count_nonzero(u)==np.count_nonzero(w)==0,'Initial forced modes nonzero')
    return report
