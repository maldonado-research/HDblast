"""Drive coexact Maxwell probes through the recorded round five-dimensional collar."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np


HERE=Path(__file__).resolve().parent
LAMBDAS=(2,5,10,20,50,103,258,515,1030,7)
BASECOLS=('t','a','K','E','B','I','D','phi')


def check_fixed_power_launch():
    """Check the source construction in a Hopf orthonormal frame."""
    knots,weights=np.polynomial.legendre.leggauss(32)
    u=(knots+1.)/2.
    chi=2.*np.pi*u
    for kappa in (2,5,7,10,20):
        amplitude=.2/np.sqrt(kappa-1.)
        c=np.sqrt(kappa-1.)*u**((kappa-2.)/2.)*np.exp(1j*(kappa-2)*chi)
        pattern=np.column_stack([np.zeros_like(u),c.real,-c.imag])
        norm2=np.einsum('ij,ij->i',pattern,pattern)
        emission=amplitude*pattern
        emission[:,0]=np.sqrt(1.-amplitude**2*norm2)
        if not np.allclose(np.einsum('ij,ij->i',emission,emission),1.,rtol=0.,atol=1.e-14):
            raise ValueError('Nonuniform pointwise source power')
        projection=.5*np.dot(weights,np.einsum('ij,ij->i',emission,pattern))
        if not np.isclose(projection,amplitude,rtol=1.e-13,atol=0.):
            raise ValueError('Transverse launch projection')


def build(prior, temp, lambdas, gauge_control=False):
    base=prior/'prior'
    helpers=runpy.run_path(str(prior/'feedback/feedback.py'))
    source=runpy.run_path(str(base/'envelope/build.py'))['transform'](
        (base/'coupled/joined_forward.cpp').read_text(),centered=True)
    source=runpy.run_path(str(prior/'boost/evolve.py'))['make_source'](
        source,helpers,(prior/'feedback/angular.hpp').read_text())
    def change(old,new):
        nonlocal source
        if source.count(old)!=1:
            raise ValueError('Unexpected ancestor source at '+old)
        source=source.replace(old,new)
    source='#include <fstream>\n'+source
    n=len(lambdas)+1
    modes=','.join(map(str,lambdas+(50,)))
    change('constexpr int P = 4;',f'''constexpr int P = 4, probe_count={n}, unforced={n-1};
constexpr double probe_lambda[probe_count]={{{modes}}};
constexpr int channels=3+3*probe_count;''')
    change('    COUNT\n','''    PROBE_C, PROBE_MOM=PROBE_C+probe_count, PROBE_NET=PROBE_MOM+probe_count,
    PROBE_DRIVER=PROBE_NET+probe_count, COUNT=PROBE_DRIVER+probe_count
''')
    change('return COUNT + 3 * n + channel;', 'return COUNT + channels*n + channel;')
    if source.count('COUNT + 3 * cells * P')!=2: raise ValueError('Mesh allocation')
    source=source.replace('COUNT + 3 * cells * P','COUNT + channels*cells*P')
    change('channel < 3;', 'channel < channels;')
    change('COUNT + 3 * (P * from + j) + channel','COUNT + channels * (P * from + j) + channel')
    change('    void rhs(const vector<double>',(HERE/'mode_rhs.hpp').read_text()+'\n    void rhs(const vector<double>')
    change('            outgoing = field[P - 1].o;', '''            probe_cell(y,dy,cell,geom,emitted,w.k);
            outgoing = field[P - 1].o;''')
    change('        dy[A_W] = a * h;', '''        probe_wall(y,dy,w,a,h);
        dy[A_W] = a * h;''')
    change('if (argc != 34)', 'if (argc != 36)')
    change('int collision_mode = 0;', 'double branching=1.;\nint collision_mode = 0;')
    change('collision_mode=std::stoi(argv[33]);', '''collision_mode=std::stoi(argv[33]);
    branching=std::stod(argv[34]); std::ofstream trace(argv[35],std::ios::binary);
    trace.exceptions(std::ios::badbit|std::ios::failbit);
    double last_probe_print=-1.;''')
    change('return proper_density*boosted_rate(', 'return branching*proper_density*boosted_rate(')
    anchor='        const double margin=gas.intrinsic_random-.9*wnew.e;'
    change(anchor,'''        if(t<.020 || t>6.24 || t-last_probe_print>=.0004) {
            const double record[]={t+dt,a,wnew.k,wnew.x,wnew.b,y[s.index(0,I)],gas.comoving_rest,y[PHI]};
            trace.write(reinterpret_cast<const char*>(record),sizeof(record));
            for(int m=0;m<probe_count;++m) {
                const double mode[]={y[PROBE_MOM+m]/wnew.k,probe_lambda[m]*y[PROBE_C+m]/a,
                    y[s.index(0,3+m)],y[PROBE_NET+m],y[PROBE_DRIVER+m]};
                trace.write(reinterpret_cast<const char*>(mode),sizeof(mode));
            }
            last_probe_print=t;
        }
'''+anchor)
    if gauge_control:
        change('constexpr int channels=', 'int probe_gauge_mode=0;\nconstexpr int channels=')
        change('struct Diag {','''Values probe_material(Values w, double psi) {
    if (probe_gauge_mode==1 && psi>=fall_phi) { w.k=prepare_gauge+width*cap; w.dk=0.; }
    if (probe_gauge_mode==2) { w.k=width; w.dk=0.; }
    return w;
}
struct Diag {''')
        change('probe_wall(y,dy,w,a,h);','probe_wall(y,dy,probe_material(w,y[PHI]),a,h);')
        change('probe_cell(y,dy,cell,geom,emitted,w.k);','probe_cell(y,dy,cell,geom,emitted,probe_material(w,y[PHI]).k);')
        change('if (argc != 36)', 'if (argc != 37)')
        change('double last_probe_print=-1.;', 'double last_probe_print=-1.; probe_gauge_mode=std::stoi(argv[36]);')
        change('t+dt,a,wnew.k,wnew.x', 't+dt,a,probe_material(wnew,y[PHI]).k,wnew.x')
        change('y[PROBE_MOM+m]/wnew.k,probe_lambda', 'y[PROBE_MOM+m]/probe_material(wnew,y[PHI]).k,probe_lambda')
    src=temp/'probe.cpp'
    src.write_text(source)
    exe=temp/'probe'
    subprocess.run(['g++','-std=c++17','-O3','-I',str(prior/'boost'),str(src),'-o',str(exe)],check=True)
    return exe


def solve_case(exe, prior, history, outdir, key, name, mesh_factor, temp, lambdas, trace_folder, gauge_mode=0):
    n=len(lambdas)+1
    meta=json.loads((history/'markers.json').read_text())
    spec=meta['runs'][key]
    with np.load(history/(key+'.npz')) as tab:
        times=tab['t']; oldE=tab['X']; oldB=tab['B']; ab=tab['a'][0]
    config=json.loads((prior/'prior/envelope/ENVELOPE.json').read_text())['parameters']['mass_4e12']
    path=temp/(name+'.bin')
    wrapper=temp/name
    args=['2',repr(meta['action']),str(prior/'boost/BOOST_TABLE.dat'),'1',repr(spec['channel_fraction']),str(path)]
    if gauge_mode:
        args.append(str(gauge_mode))
    wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(args)+'\n')
    wrapper.chmod(0o755)
    end=spec['markers']['0.05']['time']+.0004
    cfg=dict(config,cfl=config['cfl']*(.5 if key.endswith('half') else 1.),sample=0,end=end,
             early_cells=config['early_cells']*mesh_factor,late_cells=config['late_cells']*mesh_factor)
    run=runpy.run_path(str(prior/'prior/coupled/joined.py'))['run']
    try: run(wrapper,cfg)
    except subprocess.CalledProcessError as exc: raise RuntimeError(exc.stderr) from exc
    raw=np.fromfile(path,dtype='=f8').reshape(-1,len(BASECOLS)+n*5)
    if not np.isfinite(raw).all(): raise ValueError('Nonfinite probe')
    data=raw[:,:len(BASECOLS)].T
    t,a,gauge,E0,B0,I0,D,phi=data
    E,B,I,net,driver=raw[:,len(BASECOLS):].reshape(-1,n,5).transpose(2,1,0)
    def lin_each(values,times):
        return np.array([np.interp(times,t,row) for row in values])
    actual_birth=t[np.flatnonzero(D>0)[0]-1]
    pre=t<=min(actual_birth,times[0])
    tnodes=times[times<=end]
    if mesh_factor==1:
        max_mean_error=max(np.max(abs(np.interp(tnodes,t,E0)-oldE[:len(tnodes)])),
                           np.max(abs(np.interp(tnodes,t,B0)-oldB[:len(tnodes)])))
        if max_mean_error>5e-6: raise ValueError(('Mean history reproduction',max_mean_error))
    else: max_mean_error=None
    if trace_folder:
        traces=outdir/trace_folder
        traces.mkdir(exist_ok=True)
        np.savez_compressed(traces/(name+'.npz'),t=tnodes,k=np.array(lambdas)/ab,kappa=lambdas,
            E_cos=lin_each(E[:-1],tnodes), B_sin=lin_each(B[:-1],tnodes), I_cos=lin_each(I[:-1],tnodes))
    wall_energy=.5*a**4*gauge*(E*E+B*B)
    scale=np.max(wall_energy+abs(net)+abs(driver),axis=1)
    error=np.max(abs(wall_energy-net-driver),axis=1)
    erel=error[:-1]/scale[:-1]
    pre_error=max(np.max(abs(E[0,pre]-E0[pre])),np.max(abs(B[0,pre]-B0[pre])),
                  np.max(abs(I[0,pre]-I0[pre])))
    if gauge_mode==0 and pre_error>1e-8 and mesh_factor==1: raise ValueError(('lowest mode',pre_error))
    if max(erel)>3e-5 or np.any(E[-1]!=0) or np.any(B[-1]!=0) or np.any(I[-1]!=0):
        raise ValueError(('Wall power or unforced mode',max(erel)))
    def sample(values,when): return float(np.interp(when,t,values))
    measures={}
    for mi,kappa in enumerate(lambdas):
        mode={'k_birth':float(kappa/ab),'energy_balance_relative':float(erel[mi]),
              'peak_incident_wall':float(np.max(abs(I[mi]))),'peak_electric':float(np.max(abs(E[mi]))),
              'peak_magnetic':float(np.max(abs(B[mi])))}
        for label,when in [('birth',times[0]),('mesh_birth',actual_birth)]+[(w,s['time']) for w,s in spec['markers'].items()]:
            mode[label]={field:sample(values,when) for field,values in {
                'E':E[mi],'B':B[mi],'I':I[mi],'net_incoming_work':net[mi],'gauge_work':driver[mi],
                'comoving_wall_energy':wall_energy[mi]}.items()}
        measures[str(kappa)]=mode
    return dict(reference_birth_left_time=float(times[0]),actual_birth_left_time=float(actual_birth),
        original_wall_history_max_abs_error=max_mean_error,lowest_mode_prebirth_max_abs_error=float(pre_error),
        max_absolute_unforced_field=float(np.max(abs(raw[:,-5:]))),
        launch='unit coefficient of the reference emitted characteristic in a normalized independent harmonic',
        modes=measures)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--spectrum',action='store_true',help='Scan integer curl eigenvalues 2 to 40 on finer time and radial meshes.')
    p.add_argument('--gauge-controls',action='store_true',help='Rerun with the final gauge release suppressed, or with fixed gauge from launch.')
    p.add_argument('--cases',nargs='+',choices=['dedicated_native','dedicated_half','shared_native','shared_half','dedicated_spatial'],default=None)
    opt=p.parse_args(); prior=opt.prior.resolve(); out=opt.output_dir.resolve();out.mkdir(exist_ok=True,parents=True)
    if opt.spectrum and opt.gauge_controls: p.error('Choose one scan')
    check_fixed_power_launch()
    cases=opt.cases or ['dedicated_native','dedicated_half','shared_native','shared_half','dedicated_spatial']
    lambdas=tuple(range(2,41)) if opt.spectrum else LAMBDAS
    active=[(name,0) for name in cases]
    if opt.gauge_controls:
        cases=opt.cases or ['dedicated_half','dedicated_spatial']
        active=[(name,mode) for name in cases for mode in (1,2)]
    with tempfile.TemporaryDirectory(prefix='hdb-coexact-') as wd:
        temp=Path(wd); exe=build(prior,temp,lambdas,opt.gauge_controls)
        def one(task):
            name,mode=task
            key='dedicated_native' if name=='dedicated_spatial' else name
            spatial=name=='dedicated_spatial'
            if mode: name+=('_hold_release' if mode==1 else '_fixed')
            print('Replaying',name,flush=True)
            value=solve_case(exe,prior,opt.history,out,key,name,2 if spatial else 1,temp,lambdas,
                             None if mode else ('spectrum_traces' if opt.spectrum else 'traces'), mode)
            print('Completed',name,flush=True)
            return name,value
        with ThreadPoolExecutor(opt.jobs) as pool: results=dict(pool.map(one,active))
    product='GAUGE.json' if opt.gauge_controls else ('SPECTRUM.json' if opt.spectrum else 'PROBE.json')
    (out/product).write_text(json.dumps(results,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__': main()
