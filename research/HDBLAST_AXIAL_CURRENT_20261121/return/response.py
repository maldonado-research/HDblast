"""Return characteristic for three time-dependent physical electric coefficients at the wall."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent


def prepare(geometry, work):
    with np.load(geometry) as file:
        g={key:file[key] for key in file.files}
    time=g['t']
    np.column_stack((time,g['geo'].reshape(time.size,-1))).astype('=f8').tofile(work/'geometry')
    exe=work/'collar'
    subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off',str(HERE/'collar.cpp'),'-o',str(exe)],check=True)
    return g


def calculate(g, work, t, e, refine=1, cfl=.06):
    t,e=np.asarray(t),np.asarray(e)
    if e.shape!=(len(t),3) or not np.isfinite(e).all() or np.any(np.diff(t)<=0):
        raise ValueError("Expected three electric components on increasing times")
    np.column_stack((t,e)).astype('=f8').tofile(work/'drive')
    cmd=[work/'collar',work/'geometry',work/'drive',work/'result',g['geo'].shape[1],refine,float(g['depth']),cfl]
    subprocess.run(list(map(str,cmd)),check=True)
    raw=np.fromfile(work/'result',dtype='=f8').reshape(len(t),3,7)
    if not np.isfinite(raw).all():
        raise ValueError('Nonfinite collar response')
    I,out,energy,wall,tube,metric,numerical=raw.transpose(2,0,1)
    residual=energy-wall+tube-metric+numerical
    scale=(np.abs(wall)+tube+np.abs(metric)+numerical+energy).max(axis=0)
    peak=np.max(np.abs(e),axis=0)
    ratio=lambda a,b: np.divide(a,b,out=np.zeros_like(a),where=b>0)
    summary=dict(cells=refine*g['geo'].shape[1],cfl=cfl,
        peak_wall_change=peak.tolist(),peak_return_I=np.max(np.abs(I),axis=0).tolist(),
        peak_return_over_peak_wall=ratio(np.max(np.abs(I),axis=0),peak).tolist(),
        peak_bulk_balance_residual_over_scale=ratio(np.max(abs(residual),axis=0),scale).tolist(),
        end_bulk_energy=energy[-1].tolist(),wall_to_bulk=wall[-1].tolist(),through_tube=tube[-1].tolist(),
        work_from_background=metric[-1].tolist(),dg_dissipation=numerical[-1].tolist())
    return I,summary


def continuum_cap(g,duration):
    """Bound the continuous characteristic equations for this tabulated metric."""
    data=g['geo']
    # On each spatial element every metric component is a cubic in the nodal model.
    # The Bernstein coefficients enclose all interior values, including overshoots.
    xs=(np.array([-1.,-np.sqrt(.2),np.sqrt(.2),1.])+1)/2
    basis=np.column_stack(((1-xs)**3,3*xs*(1-xs)**2,3*xs**2*(1-xs),xs**3))
    coefs=np.einsum('ij,tcjk->tcik',np.linalg.inv(basis),data)
    lo=coefs.min(axis=(0,1,2)); hi=coefs.max(axis=(0,1,2))
    margin=64*np.finfo(float).eps*np.maximum(1,np.maximum(abs(lo),abs(hi)))
    lo-=margin; hi+=margin
    rmin=lo[0];rmax=hi[0];nmin=lo[3];nmax=hi[3]
    if rmin<=0 or nmin<=0:raise ValueError('No positive causal speeds guaranteed')
    tau=float(g['depth'])/nmin
    grad,expansion=[nmax*max(abs(lo[j]),abs(hi[j]))/rmin for j in (1,2)]
    nx=max(abs(lo[4]),abs(hi[4]))
    # Bound N(rx+u)/(2r) separately to retain cancellation in the inward source.
    mixed=.5*nmax*max(abs(lo[1]+lo[2]),abs(hi[1]+hi[2]))/rmin
    diag=nx+1.5*(grad+expansion)
    potential=2.*nmax/rmin**3 * duration*nmax*rmax
    matrix=tau*(np.array([[diag,mixed],[.5*(grad+expansion),diag]])+potential)
    matrix[1,0]+=1.  # outgoing wall datum is E - I
    radius=max(abs(np.linalg.eigvals(matrix)))
    if radius>=1:raise ValueError("Insufficient bound")
    cap=np.linalg.solve(np.eye(2)-matrix,[0.,1.])
    return dict(max_one_way_time=tau,min_one_way_time=float(g['depth'])/nmax,
                I_peak_over_E_peak_upper=float(cap[0]),O_peak_over_E_peak_upper=float(cap[1]),
                duration=duration,bounding_matrix=matrix.tolist(),spectral_radius=float(radius),
                metric_lower=lo.tolist(),metric_upper=hi.tolist())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--snapshot',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--drive-output',type=Path)
    p.add_argument('--refine',type=int,default=1)
    p.add_argument('--cfl',type=float,default=.06)
    opt=p.parse_args()
    with tempfile.TemporaryDirectory() as folder:
        work=Path(folder)
        g=prepare(opt.geometry,work)
        if opt.snapshot:
            with np.load(opt.snapshot) as file:
                t=file['t'];e=file['delta_E']
        else:
            t=g['t'];e=np.zeros((len(t),3));e[:,2]=g['drive']
        inc,report=calculate(g,work,t,e,opt.refine,opt.cfl)
        report['characteristic_cap']=continuum_cap(g,t[-1]-t[0])
        if not opt.snapshot:
            report['manufactured_max_error']=float(np.max(abs(inc[:,2]-g['reference'])))
            report['manufactured_peak_reference']=float(np.max(abs(g['reference'])))
        if opt.drive_output:
            np.savez_compressed(opt.drive_output,t=t,delta_incoming=inc)
        opt.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
