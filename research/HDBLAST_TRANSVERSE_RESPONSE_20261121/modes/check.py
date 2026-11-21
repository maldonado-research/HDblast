"""Check the transverse solver against an independent constant-coefficient matrix exponential."""
import argparse
import json
import math
from pathlib import Path
import struct
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parent
Q=1.e7
K=2.54


def frozen(exe,folder,density,velocity,gyro,physical_k,duration,steps):
    # Choose a=1 and one uniformly populated pair cohort present right at t=0.
    gamma=1/math.sqrt(1-velocity**2)
    path=folder/"frozen.bin"
    with path.open("wb") as f:
        f.write(struct.pack("=iiid",1,steps,0,duration/steps))
        np.tile([1.,gyro*gamma/Q,0.],(2*steps+1,1)).astype("=f8").tofile(f)
        f.write(struct.pack("=dddii",density,-velocity*gamma,0.,0,0))
    line=subprocess.check_output([str(exe),str(path),str(physical_k),"1"],text=True)
    vals=np.fromstring(line,sep=" ")
    # Variables are physical E/B followed by delta p/(m Q) and delta(number)/Q.
    A=np.zeros((8,8)); by,ez,bz,ey,px,py,pz,d=range(8)
    A[by,ez]=A[bz,ey]=physical_k
    A[ez,by]=A[ey,bz]=-physical_k
    s=Q**2*density/K
    A[ez,pz]=-s/gamma**3;A[ez,d]=-s*velocity;A[ey,py]=-s/gamma
    A[px,by]=-velocity;A[px,py]=gyro
    A[py,ey]=1.;A[py,px]=-gyro;A[pz,ez]=1.;A[d,px]=physical_k/gamma
    sol=expm(duration*A)
    measured=vals[-16:].reshape(4,4)
    field_error=float(np.linalg.norm(measured-sol[:4,:4],2)/max(1.,np.linalg.norm(sol[:4,:4],2)))
    gram=vals[9:25].reshape(4,4)
    exact=(sol[[px,py],:4].T@sol[[px,py],:4])*Q**2
    gram_error=float(np.linalg.norm(gram-exact,2)/max(1.,np.linalg.norm(exact,2)))
    return {"field_matrix_relative_error":field_error,"momentum_gram_relative_error":gram_error,
            "max_real_eigenvalue":float(max(np.linalg.eigvals(A).real)),"steps":steps}


def free_stream_bound(state,marker,k):
    """Absolute upper on the free unscattered Fourier phase from finite injected transverse speed."""
    use=state['t']<=marker['time']
    # Linear interpolation on each archived segment attains endpoint extrema.
    aa=np.append(state['a'][use],marker['a'])
    ee=np.append(state['X'][use],marker['X'])
    if not (np.all(ee<0) and np.all(np.diff(aa)>=0)): raise ValueError('Bound assumes expanding and single-signed field')
    acc=Q*min(abs(ee))*min(aa)/max(aa)
    b=math.sqrt(1.05**2-1)
    duration=marker['time']-state['t'][0]
    integral=min(duration,b/acc)+math.asinh(max(0.,acc*duration-b))/acc
    return float(k*b/min(aa)*integral)


def check_three_dimensional_bounds():
    """Verify the generator estimates used to include finite pitch and arbitrary k."""
    rng=np.random.default_rng(81927)
    rotate=np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,0.]])
    guide_fraction=0.;jac_error=0.;force=0.;skew_error=0.
    b=math.sqrt(1.05**2-1.)
    for i,axial in enumerate(np.r_[0.,np.geomspace(.01,8.e4,512)]):
        angle=rng.uniform(-np.pi,np.pi)
        perp=b*rng.uniform(.01,1.)
        p=np.array([perp*np.cos(angle),perp*np.sin(angle),axial*(-1 if i%2 else 1)])
        gamma=np.hypot(1.,np.linalg.norm(p));v=p/gamma
        derivative=(np.eye(3)-np.outer(v,v))/gamma
        gen=rotate@derivative
        budget=perp/gamma**2
        guide_fraction=max(guide_fraction,max(np.linalg.eigvalsh((gen+gen.T)/2))/budget)
        # Centered force difference, selected mid-energy region to avoid cancellation.
        if axial<=100.:
            h=1e-5
            def lorentz(u):return rotate@u/np.hypot(1.,np.linalg.norm(u))
            fd=np.column_stack([(lorentz(p+e*h)-lorentz(p-e*h))/(2*h) for e in np.eye(3)])
            jac_error=max(jac_error,np.linalg.norm(fd-gen,2)/np.linalg.norm(gen,2))
        cross=lambda a: np.array([[0.,-a[2],a[1]],[a[2],0.,-a[0]],[-a[1],a[0],0.]])
        force=max(force,np.linalg.norm(np.concatenate((np.eye(3),cross(v)),axis=1),2))
        curl=cross(rng.normal(size=3))
        maxwell=np.block([[np.zeros((3,3)),1j*curl],[-1j*curl,np.zeros((3,3))]])
        skew_error=max(skew_error,np.linalg.norm(maxwell+maxwell.conj().T))
    if guide_fraction>1+1e-12 or jac_error>1e-8 or force>math.sqrt(2)+1.e-12 or skew_error>1e-12:
        raise ValueError("Three-dimensional linear generator failed a norm check")
    return dict(probes=513,largest_symmetric_guide_eigenvalue_over_bound=guide_fraction,
                guide_jacobian_relative_error=jac_error,max_lorentz_operator_norm=force,
                max_vacuum_skew_hermitian_defect=skew_error)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--history',type=Path,default=ROOT.parent/'history')
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_text())
    phases={}
    for name,row in meta['runs'].items():
        with np.load(args.history/(name+'.npz')) as s: phases[name]=free_stream_bound(s,row['markers']['0.01'],32768)
    tests={}
    with tempfile.TemporaryDirectory() as f:
        d=Path(f);exe=d/'modes'
        subprocess.run(['g++','-std=c++17','-O3',str(ROOT/'evolve.cpp'),'-o',str(exe)],check=True)
        for name,wp,v,gyro,k,t in [('vacuum',1e-10,.98,5,62,.1),('plasma_at_rest',20.,0.,6.,12.,.08),
                ('growing_no_guide',20.,.95,0.,25.,.14),('growing_guide',20.,.95,7.,30.,.13),
                ('strong_guide',20.,.90,44.,30.,.13)]:
            w=K*(wp/Q)**2
            tests[name]=[frozen(exe,d,w,v,gyro,k,t,n) for n in (256,1024)]
            fine=tests[name][-1]
            if max(fine['field_matrix_relative_error'],fine['momentum_gram_relative_error'])>2e-5:
                raise ValueError((name,fine))
    if tests['growing_no_guide'][-1]['max_real_eigenvalue']<=.01:raise ValueError('Expected positive mode')
    if tests['strong_guide'][-1]['max_real_eigenvalue']>1e-4:raise ValueError('Expected stable constant background')
    args.output.write_text(json.dumps({'constant_background':tests,'free_stream_phase_ceiling_k32768':phases,
                                     'three_dimensional_generator':check_three_dimensional_bounds()},
                                      indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__': main()
