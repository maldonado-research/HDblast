#!/usr/bin/env python3
"""Read-only archive diagnostics and exact finite-difference controls; NumPy only.

No evolution run. Archive path can be provided as the first argument.
"""
from fractions import Fraction as F
from pathlib import Path
from io import BytesIO
import hashlib,json,sys,zipfile
import numpy as np

archive=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921.zip')
prefix='HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/'
checks=[]
def record(name,ok,**extra):
    checks.append(dict(name=name,pass_check=bool(ok),**extra))
    if not ok:raise RuntimeError(name)

# Coordinate spacing h=1; f(xi)=xi^3 has f'(0)=f''(0)=0.
values={k:F(k**3) for k in range(-3,1)}
values[1]=values[-1];values[2]=values[-2]
d1=(values[-2]-8*values[-1]+8*values[1]-values[2])/12
d2=(-values[-2]+16*values[-1]-30*values[0]+16*values[1]-values[2])/12
record('even-plus-linear closure enforces derivative',d1==0)
record('even-plus-linear closure has cubic second-derivative error',d2==F(-4,3),value=str(d2),scaled_error='-(4/3)*h for f(z)=z^3')

# Mapped V2 constructs the third ghost by linear extrapolation rather than
# consistent even extension or polynomial extrapolation.
q={k:F(k*k) for k in range(-3,3)};q[3]=2*q[2]-q[1]
delta6=q[-3]-6*q[-2]+15*q[-1]-20*q[0]+15*q[1]-6*q[2]+q[3]
record('V2 third ghost creates KO forcing on quadratic',delta6==-2,value=str(delta6),Q_over_epsilon=str(delta6/64))
q[3]=F(9)
record('consistent polynomial third ghost removes quadratic KO forcing',q[-3]-6*q[-2]+15*q[-1]-20*q[0]+15*q[1]-6*q[2]+q[3]==0)

# Shell identity with arbitrary nonzero example values.
eb,sig,sig1,Bt,pt=map(F,[2,3,5,7,11])
Az=Bz=eb*sig/6;pz=-eb*sig1/2;Atz=eb*(sig*Bt+sig1*pt)/6
M=-3*Atz+3*Az*Bt-pt*pz
wrongM=3*Az*Bt-pt*pz
record('correct velocity slope cancels shell M',M==0)
record('zero velocity slope produces false nonzero shell M',wrongM!=0,false_residual=str(wrongM))

with zipfile.ZipFile(archive) as zf:
    solver_bytes=zf.read(prefix+'rolloff5d.py')
    data_name=prefix+'referee_numerics/runs/ref_t1e3_dz2.5e-4_timeseries.npz'
    raw=zf.read(data_name)
    with np.load(BytesIO(raw),allow_pickle=False) as data:
        zz=data['z'];rho=data['rho'];rec=data['rec'];rb=float(rho[-1])
    record('archived radial nodes ordered',bool(np.all(np.diff(zz)>0)))
    record('archived warp factor positive',bool(np.all(rho>0)))
    table=[]
    for t in [.01,.1,.5,1.,2.,2.5,3.]:
        rr=float(np.interp(-t,zz,rho))
        gain=float(np.exp(-3*t)*(rb/rr)**3)
        table.append(dict(time=t,ray_position=-t,rho=rr,Cplus_unweighted_gain=gain))
    record('constraint gain near old reported front exceeds ten thousand',table[-2]['Cplus_unweighted_gain']>10000)
    static_phi=float(rec[0,1])
    # Registered-v2 seed supplied in its log; this check is arithmetic and
    # the numeric value is explicitly a rounded log datum, not recomputed seed.
    seed_from_v2_log=8.979e-5
    record('registered dc1e-4 automatic growth-fit window is empty',30*seed_from_v2_log>2e-3,
           lower_threshold=30*seed_from_v2_log,upper_threshold=.002)

result={'status':'PASS','scope':'Archive diagnostics and exact stencil controls; no evolution replay or interval bounds.',
        'check_count':len(checks),'checks':checks,'archive':str(archive),
        'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
        'solver_sha256':hashlib.sha256(solver_bytes).hexdigest(),
        'background_archive_member':data_name,'background_member_sha256':hashlib.sha256(raw).hexdigest(),
        'rho_b':rb,'constraint_transport_amplification':table,
        'transport_scope':'Static-background formula evaluated on archived numerical rho with linear interpolation; explanatory diagnostic, not proof of the cause of the original run failure.',
        'self_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('ARCHIVED_NUMERICS_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'check_count':result['check_count'],'gain_at_t2p5':table[-2]['Cplus_unweighted_gain']}))
