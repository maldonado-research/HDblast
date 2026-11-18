"""Conformally coupled massive spectator on fixed B1 histories. Requires NumPy."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np

PIN = '9a7069e44c379e68555511e82fe841e0f01a72b23e77fa5ec69caa1373beaf70'
RHO_B = 7.838285538073243
G, PHI_STAR = 100., 0.5
B_COUPLING = 6.49e-6
CUTS = (0.09, 0.12, 0.15)
SELECTED = (0., 10., 18., 25., 50., 100.)
BANDS = ((5., 18.), (5., 80.))
# Both Gauss nodes on [0,1], and weights for the two commutator-free exponentials.
MINUS, PLUS = (0.5-math.sqrt(3)/6, 0.5+math.sqrt(3)/6)
FRONT, BACK = (0.25+math.sqrt(3)/6, 0.25-math.sqrt(3)/6)


class Wall:
    def __init__(self, archive, group, grid, case):
        folder = f'{group}/{grid}/'
        raw = archive.read(folder+case+'.csv')
        meta = json.loads(archive.read(folder+'metadata.json'))
        if hashlib.sha256(raw).hexdigest() != meta['cases'][case]['sha256']:
            raise ValueError('Wall source digest mismatch')
        rows = list(csv.DictReader(io.StringIO(raw.decode())))
        self.t = np.array([float(r['T']) for r in rows])
        self.fields = np.array([[float(r[k]) for k in ('A_wall','B_wall','phi_wall')] for r in rows])
        self.slopes = np.array([[float(r[k]) for k in ('A_T_wall','B_T_wall','phi_T_wall')] for r in rows])
        self.a_initial = self.fields[0,0]
        if not np.all(np.diff(self.t)>0) or self.t[0]!=0 or not np.isfinite(self.fields).all():
            raise ValueError('Invalid time history')

    def sample(self, j, u):
        dt = self.t[j+1]-self.t[j]
        y0, y1 = self.fields[j:j+2]
        d0, d1 = self.slopes[j:j+2]
        sec = (y1-y0)/dt
        c2 = 3*sec-2*d0-d1
        c3 = -2*sec+d0+d1
        y = y0+dt*u*(d0+u*(c2+u*c3))
        dy = d0+u*(2*c2+u*3*c3)
        a = math.exp(y[0]-self.a_initial)
        lapse = math.exp(y[1])/RHO_B
        ell = lapse/a
        mu = G*abs(y[2]-PHI_STAR)
        mu_t = G*math.copysign(1.,y[2]-PHI_STAR)*dy[2]
        a_eta = a*dy[0]/ell
        M = a*mu
        M_eta = a_eta*mu+a*mu_t/ell
        return State(a,lapse,ell,M,M_eta,a_eta,mu)


class State:
    __slots__ = ('a','lapse','ell','M','M_eta','a_eta','mu')
    def __init__(self,a,lapse,ell,M,M_eta,a_eta,mu):
        self.a,self.lapse,self.ell,self.M,self.M_eta,self.a_eta,self.mu=a,lapse,ell,M,M_eta,a_eta,mu


def wave_number_grid():
    points=list(SELECTED)
    rules=[]
    for lo,hi in BANDS:
        for order in (16,32):
            nodes, weights=np.polynomial.legendre.leggauss(order)
            k=0.5*((hi-lo)*nodes+(hi+lo))
            segment=slice(len(points),len(points)+len(k))
            points.extend(k)
            rules.append((lo,hi,order,segment,(hi-lo)*weights/2))
    return np.array(points),rules


def rotate(matrix, alpha, beta, duration):
    freq=np.sqrt(alpha*beta)
    cs=np.cos(freq*duration)
    sn=np.sin(freq*duration)/freq
    a,b,c,d=matrix
    return (cs*a+alpha*sn*c, cs*b+alpha*sn*d,
            cs*c-beta*sn*a, cs*d-beta*sn*b)


def propagate(segment, k, matrix, duration, midpoint):
    sm,sp=segment
    qm=k*k+sm.M**2
    qp=k*k+sp.M**2
    if midpoint:
        return rotate(matrix,sm.ell,sm.ell*qm,duration)
    # For time ordering the first exponential weights the earlier Gauss node more.
    for early,late in ((FRONT,BACK),(BACK,FRONT)):
        alpha=early*sm.ell+late*sp.ell
        beta=early*sm.ell*qm+late*sp.ell*qp
        matrix=rotate(matrix,alpha,beta,duration)
    return matrix


def evolve(wall,k,method='cf4'):
    zero=np.zeros_like(k)
    one=np.ones_like(k)
    matrix=(one,zero,zero,one)
    initial=wall.sample(0,0.)
    result={}
    s=eta=0.
    cut_index=0
    emax=0.
    for j,dt in enumerate(np.diff(wall.t)):
        u=0.
        while u < 1:
            def quadrature(v):
                span=v-u
                return (wall.sample(j,u+MINUS*span),wall.sample(j,u+PLUS*span))
            sm,sp=quadrature(1.)
            ds=(1-u)*dt*(sm.lapse+sp.lapse)/2
            at_cut=s+ds>=CUTS[cut_index]-2e-14
            v=1.
            if at_cut and s+ds>CUTS[cut_index]+2e-14:
                a,b=u,1.
                for _ in range(46):
                    test=(a+b)/2
                    mm,pp=quadrature(test)
                    s_test=s+(test-u)*dt*(mm.lapse+pp.lapse)/2
                    if s_test < CUTS[cut_index]:
                        a=test
                    else:
                        b=test
                v=(a+b)/2
                sm,sp=quadrature(v)
            dur=dt*(v-u)
            if method=='midpoint':
                center=wall.sample(j,(u+v)/2)
                matrix=propagate((center,center),k,matrix,dur,True)
            else:
                matrix=propagate((sm,sp),k,matrix,dur,False)
            s+=dur*(sm.lapse+sp.lapse)/2
            eta+=dur*(sm.ell+sp.ell)/2
            state=wall.sample(j,v)
            emax=max(emax,*(math.hypot(18./x.a,x.mu) for x in (sm,sp,state,initial)))
            if at_cut:
                result[CUTS[cut_index]]=(state,matrix,eta,wall.t[j]+v*dt,emax)
                cut_index+=1
                if cut_index==len(CUTS):
                    return initial,result
            u=v
    raise ValueError('Incomplete requested proper-time cuts')


def frequency(k, state):
    w=np.sqrt(k*k+state.M**2)
    return w, state.M*state.M_eta/(2*w*w)


def oscillator(k,initial,state,matrix,seed):
    w0,d0=frequency(k,initial)
    v0=1/np.sqrt(2*w0)
    p0=(-1j*w0-(d0 if seed=='wkb' else 0))*v0
    a,b,c,d=matrix
    return a*v0+b*p0,c*v0+d*p0


def energy_data(k,pair,quiet,seed):
    s,m,si=pair
    sq,mq,sqi=quiet
    v,p=oscillator(k,si,s,m,seed)
    vq,pq=oscillator(k,sqi,sq,mq,seed)
    w,D=frequency(k,s)
    wq,Dq=frequency(k,sq)
    def occupation(v,p,w,D):
        return np.abs(w*v-1j*(p+D*v))**2/(2*w)
    n=occupation(v,p,w,0.)
    nq=occupation(vq,pq,wq,0.)
    wkb=occupation(v,p,w,D)
    wkbq=occupation(vq,pq,wq,Dq)
    relative=np.abs(vq*p-pq*v)**2
    minimax=relative/(2*(np.sqrt(1+relative)+1))
    field_variance=np.abs(v/s.a)**2
    quiet_variance=np.abs(vq/sq.a)**2
    chi,Pi=v/s.a,s.a*p-s.a_eta*v
    chiq,Piq=vq/sq.a,sq.a*pq-sq.a_eta*vq
    relative_field=np.abs(chiq*Pi-Piq*chi)**2
    # A conditional number relative to the quiet end-point frequency, common to both.
    shared_occupation_difference=occupation(v,p,wq,0.)-nq
    qvol=sq.a**4
    vol=s.a**4
    energy=n*w/vol
    energy_quiet=nq*wq/qvol
    vacuum_delta=(w/vol-wq/qvol)/2
    return dict(relative=relative,relative_field=relative_field,minimax_common_occupation=minimax,
                field_variance=field_variance,field_variance_quiet=quiet_variance,
                delta_field_variance=field_variance-quiet_variance,n=n,n_quiet=nq,
                delta_n=n-nq,wkb_n=wkb,wkb_n_quiet=wkbq,delta_wkb_n=wkb-wkbq,
                shared_delta_n=shared_occupation_difference,
                energy_excitation=energy,energy_excitation_quiet=energy_quiet,
                delta_energy_excitation=energy-energy_quiet,
                delta_energy_zero_point=vacuum_delta,
                delta_energy_total=energy-energy_quiet+vacuum_delta,
                wronskian_error=abs(-2*np.imag(np.conjugate(v)*p)-1))


def compare(k,rules,run,quiet):
    i,sol=run
    iq,solq=quiet
    if max(abs(i.a-iq.a),abs(i.M-iq.M),abs(i.M_eta-iq.M_eta))>1e-12:
        raise ValueError('Initial states differ')
    result={}
    for t in CUTS:
        st,mat,eta,T,emax=sol[t]
        qt,mt,etaq,Tq,_=solq[t]
        a,b,c,d=mat
        vals={'T':T,'quiet_T':Tq,'conformal_clock':eta,'quiet_conformal_clock':etaq,
              'max_screened_physical_frequency_at_sampled_nodes':emax,
              'symplectic_max_error':float(np.max(np.abs(a*d-b*c-1)))}
        for seed in ('instant','wkb'):
            arr=energy_data(k,(st,mat,i),(qt,mt,iq),seed)
            selected={}
            for j,k0 in enumerate(SELECTED):
                selected[str(k0)]={name:float(vec[j]) for name,vec in arr.items()}
            # Preserve the complex canonical mode amplitudes for independent covariance checks.
            for suffix,initial,state,matrix in (('pulse',i,st,mat),('quiet',iq,qt,mt)):
                v,p=oscillator(k,initial,state,matrix,seed)
                for j,k0 in enumerate(SELECTED):
                    selected[str(k0)][suffix+'_mode']=[float(v[j].real),float(v[j].imag),float(p[j].real),float(p[j].imag)]
            bands={}
            for lo,hi,n,sl,ws in rules:
                names=('relative','relative_field','minimax_common_occupation','field_variance',
                       'field_variance_quiet','delta_field_variance','n','n_quiet','delta_n','delta_wkb_n',
                       'energy_excitation','energy_excitation_quiet','delta_energy_excitation',
                       'delta_energy_zero_point','delta_energy_total')
                bands[f'{lo:g}-{hi:g}/{n}']={key:float(np.sum(ws*k[sl]**2*arr[key][sl])/(2*np.pi**2))
                                               for key in names}
            vals[seed]={'selected':selected,'band':bands,'max_wronskian_error':float(np.max(arr['wronskian_error']))}
        result[str(t)]=vals
    return result


def calculate(path):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise ValueError('Unexpected published pulse archive')
    k,rules=wave_number_grid()
    results={}
    with zipfile.ZipFile(path) as archive:
        for grid in ('g1','g2','g4','g8'):
            quiet=evolve(Wall(archive,'runs',grid,'quiet'),k)
            cases=[('plus','runs','delayed'),('minus','runs','delayed_negative')]
            if grid in ('g4','g8'):
                cases.extend((('half_plus','half_runs','delayed_positive'),('half_minus','half_runs','delayed_negative')))
            if grid=='g8':
                cases.append(('without_remote_controls','runs','uncompensated'))
            for label,group,case in cases:
                r=evolve(Wall(archive,group,grid,case),k)
                results[f'{grid}/{label}']=compare(k,rules,r,quiet)
            if grid=='g8':
                mq=evolve(Wall(archive,'runs',grid,'quiet'),k,'midpoint')
                mp=evolve(Wall(archive,'runs',grid,'delayed'),k,'midpoint')
                results['g8/midpoint_plus']=compare(k,rules,mp,mq)
    return {'input_zip_sha256':PIN,'model':{'G':G,'phi_star':PHI_STAR,'rho_b':RHO_B,
            'M5_over_H0':B_COUPLING**(-1/3),'selected_k':SELECTED,'bands':BANDS,'proper_cuts':CUTS},
            'results':results}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('wall_zip',type=Path)
    d=calculate(p.parse_args().wall_zip)
    print(json.dumps(d,indent=2,allow_nan=False))
