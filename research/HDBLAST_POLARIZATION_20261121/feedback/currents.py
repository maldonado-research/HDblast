"""Evolve homogeneous two-polarization daughter currents and their leaky field response."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.optimize import brentq


HERE=Path(__file__).resolve().parent
ORBITS=HERE.parent/'orbits'
pairs=runpy.run_path(str(ORBITS/'pairs.py'))
coarse=runpy.run_path(str(ORBITS/'coarse.py'))
K=2.54


def integrate(inputs,label,case,spec,power,steps,seed,feedback,exe,tmp):
    h,f,m=inputs[label]
    table=pairs['field_table'](h,f,m,2,2,spec[2],spec[3],steps)
    basis=table.copy()
    a1,a2,phase=spec[4:7]
    table[:,4:6]=a1*basis[:,4:6]+a2*np.cos(phase)*basis[:,6:8]
    table[:,6:8]=a2*np.sin(phase)*basis[:,6:8]
    tags=pairs['seeds'](h,m,power,seed)
    # Positions drop out in the invariant frame. The second column supplies the birth radius.
    tags[:,1]=np.interp(tags[:,0],h['t'],h['a'])
    with tempfile.TemporaryDirectory(dir=tmp) as path:
        path=Path(path);table.tofile(path/'fields');tags.tofile(path/'births')
        cmd=[str(exe),str(path/'fields'),str(path/'births'),str(m['D']),str(int(feedback)),str(path/'raw'),str(path/'history')]
        values=np.array(subprocess.check_output(cmd,text=True).split(),dtype=float)
        raw=np.fromfile(path/'raw',dtype='=f8').reshape(-1,8)
        trace=np.fromfile(path/'history',dtype='=f8').reshape(-1,10)
    if len(raw)!=2*len(tags) or len(trace)!=steps or not np.isfinite(raw).all() or not np.isfinite(trace).all():
        raise ValueError('Incomplete response trace')
    gas,injection,expansion,axial_work,supplied_work,induced_work,field,loss,iterations=values
    matter_error=(gas-injection-expansion-axial_work-supplied_work-induced_work)/gas
    field_error=(field+loss+induced_work)/max(abs(induced_work),1.e-35)
    if max(abs(matter_error),abs(field_error))>2e-9:raise ValueError(('energy ledger',matter_error,field_error))
    if np.max(abs(trace[:,0]-table[1:,0]))>1e-12:raise ValueError('Unmatched clocks')
    conformal=np.linalg.norm(trace[:,1:5],axis=1)
    if np.max(conformal-trace[:,7])>2.e-12:raise ValueError('Integrated speed-of-light current bound')
    if not np.isclose(trace[-1,8],m['D'],rtol=1e-13,atol=0.):raise ValueError('Charge weights')
    rhs=np.linalg.norm(table[1:,4:],axis=1);a=table[1:,1];induced=trace[:,1:5]/a[:,None]**2
    peak=np.max(np.linalg.norm(induced,axis=1))
    m0=pairs['pressure'](h,m,basis,spec,raw)
    c=coarse['compare'](h,m,raw)
    stats=dict(history=label,case=case,feedback=feedback,steps=steps,pairs=len(tags),seed=seed,
                launch=[float(a1),float(a2),float(phase)],emitted_energy_over_original=float(a1**2+a2**2),
                local_lateral_pressure_over_thermal=m0['pooled_lateral_pressure_over_thermal'],
                per_sign_pressure_eigenvalues_over_thermal=np.linalg.eigvalsh(m0['per_sign_pooled_lateral_covariance_over_target']).tolist(),
                energy_over_drifting_thermal=m0['energy_over_drifting_thermal'],ks_transverse=m0['ks_transverse'],
                coarse_cell_total_variation=c['coarse_total_variation'],
                final_induced_physical_E=induced[-1,:2].tolist(),final_induced_physical_B=induced[-1,2:].tolist(),
                peak_induced_physical_norm=float(peak),peak_supplied_transverse_physical_norm=float(np.max(rhs)),
                peak_induced_over_peak_supplied=float(peak/max(np.max(rhs),np.finfo(float).tiny)),
                final_induced_conformal_norm=float(conformal[-1]),final_sampled_particle_current_bound=float(trace[-1,7]),
                final_transverse_comoving_current=trace[-1,5:7].tolist(),
                ledger=dict(gas=float(gas),injected=float(injection),expansion=float(expansion),axial_work=float(axial_work),
                            supplied_transverse_work=float(supplied_work),induced_transverse_work=float(induced_work),
                            induced_field_energy=float(field),lost_to_open_boundary=float(loss),
                            particle_residual_over_gas=float(matter_error),field_residual_over_induced_work=float(field_error)),
                largest_implicit_iteration_count=int(iterations))
    return stats,raw[:,5:]


def isolated_pair_control(exe,tmp):
    """An initially drifting neutral pair gives a known signed source in both modes."""
    dt=1.e-7;mass=.15;momentum=.8
    direction=np.array([.6,.8]);charge_step=1.e7*dt
    tab=np.array([[0.,1.,0.,0.,0.,0.,0.,0.],[dt,1.,0.,0.,0.,0.,0.,0.]])
    tag=np.array([[0.,1.,0.,0.,0.,*(momentum*direction),0.]])
    tab.tofile(tmp/'test_fields');tag.tofile(tmp/'test_seed')
    def equation(avg):
        final=momentum+charge_step*avg
        velocity=(momentum+final)/(np.hypot(1.,momentum)+np.hypot(1.,final))
        return avg+charge_step*mass*velocity/(2.*K)
    midpoint=brentq(equation,-charge_step*mass/(2.*K),0.,xtol=1.e-16)
    result=subprocess.check_output([str(exe),str(tmp/'test_fields'),str(tmp/'test_seed'),str(mass),'1',
                                    str(tmp/'test_out'),str(tmp/'test_history')],text=True)
    row=np.fromfile(tmp/'test_history',dtype='=f8')
    raw=np.fromfile(tmp/'test_out',dtype='=f8').reshape(-1,8)
    # Current kick at the center, then the last half-step's exact curl and split leakage.
    e=2.*midpoint*np.exp(-dt/K)*np.cos(dt)*direction
    b=2.*midpoint*np.exp(-.5*dt/K)*np.sin(dt)*direction
    p=(momentum+charge_step*midpoint)*direction
    error=max(np.max(abs(row[1:3]-e)),np.max(abs(row[3:5]-b)),
              np.max(abs(raw[:,5:7]-raw[:,0,None]*p)))
    if not result or error>5.e-14:raise ArithmeticError('Two-axis source control')
    return dict(max_absolute_error=float(error),final_conformal_E=e.tolist(),final_conformal_B=b.tolist())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('history','fields','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--jobs',type=int,default=1)
    args=p.parse_args()
    inputs=pairs['read_inputs'](args.history,args.fields)
    tasks=[];seed=pairs['SCRAMBLE']+1
    for label in ('dedicated_half','shared_half','dedicated_native','shared_native'):
        matched=pairs['uniform_pair'](pairs['FIT'][label.split('_')[0]+'_half'][1])
        cases={'two_pressure_fit':matched,
               'equal_early_and_late_4pct':(2,2,'early','late_birth',.04/np.sqrt(2),.04/np.sqrt(2),np.pi/2),
               'equal_early_quadratures_4pct':(2,2,'early','early',.04/np.sqrt(2),.04/np.sqrt(2),np.pi/2)}
        if label.endswith('half'):
            settings=[(12,512,seed),(11,256,seed+18),(12,256,seed)]
        else:settings=[(11,256,seed+18)]
        for name,spec in cases.items():
            for power,n,s in settings:
                for feedback in (False,True):tasks.append((label,name,spec,power,n,s,feedback))
    with tempfile.TemporaryDirectory(prefix='homogeneous-currents-') as folder:
        tmp=Path(folder);exe=tmp/'homogeneous'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'homogeneous.cpp'),'-o',str(exe)],check=True)
        control=isolated_pair_control(exe,tmp)
        with ThreadPoolExecutor(args.jobs) as pool:
            outputs=list(pool.map(lambda item:integrate(inputs,*item,exe,tmp),tasks))
    rows=[]
    for i in range(0,len(outputs),2):
        passive,p0=outputs[i];active,p1=outputs[i+1]
        if passive['feedback'] or not active['feedback'] or tasks[i][:-1]!=tasks[i+1][:-1]:
            raise RuntimeError('Paired trajectory comparison')
        cosine_product=np.einsum('ij,ij->i',p0,p1)
        angles=np.arctan2(np.linalg.norm(np.cross(p0,p1),axis=1),cosine_product)
        theta=inputs[active['history']][2]['theta']
        lateral=np.linalg.norm((p1-p0)[:,:2],axis=1)/theta
        active['paired_deflection']=dict(median_angle_degrees=float(np.degrees(np.median(angles))),
                                        percentile90_angle_degrees=float(np.degrees(np.quantile(angles,.9))),
                                        rms_delta_lateral_momentum_over_theta=float(np.sqrt(np.mean(lateral**2))))
        rows.extend([passive,active])
    out=dict(runs=rows,one_pair_control=control)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False,sort_keys=True)+'\n')


if __name__=='__main__':main()
