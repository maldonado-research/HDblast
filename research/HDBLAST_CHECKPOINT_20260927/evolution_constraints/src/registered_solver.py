"""Floating-point audit solver for the registered shell; not an interval certificate.

Same continuum model as Chat13. Analytic stretched grid, quintic Hermite ghost
extension, potential derivatives cached, and independent constraint diagnostics.
The fixed stretched grid is intended for linear/small-displacement validation;
it does not certify the resolution of an outgoing nonlinear domain wall.
"""
from pathlib import Path
import argparse, json, math, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq
from scipy import sparse
from scipy.sparse.linalg import eigs

IP=1.0357712571566784
C=2/IP-4/3
UC=np.array([-1/6,4/3,-5/3,-4/9,17/18,0.,-2/27])

def W(p): return 1-p+p**3/3
def U(p): return .5*(p*p-1)**2-(2/3)*W(p)**2
def derivatives(p):
    out=[]; c=UC.copy()
    for j in range(7):
        out.append(np.polynomial.polynomial.polyval(p,c))
        c=np.arange(1,len(c))*c[1:]
    return out

def stable_potential(psi):
    p=psi-1; w=5/3-psi*psi+psi**3/3; wp=psi*(psi-2)
    return .5*wp*wp-(2/3)*w*w,wp*(2*p-(4/3)*w)

def shell_background(tdet, dc=0., rtol=2e-12):
    """Regular cone integration in shifted scalar psi=phi+1; dimensionless units."""
    c=C+dc; y0=1e-4
    def sigma(p):return 2*W(p)+tdet*(1+c*p)
    def sigma1(p):return 2*(p*p-1)+tdet*c
    def integrate(x, dense=False):
        psi_h=10**x[0]; ph=psi_h-1; yb=x[1]
        u,u1=stable_potential(psi_h)
        u2=4*ph*ph+2*(ph*ph-1)-(4/3)*((ph*ph-1)**2+2*ph*W(ph))
        k4=u1*(u2/280+u/630)
        initial=[y0-u*y0**3/36,psi_h+u1*y0*y0/10+k4*y0**4,u1*y0/5+4*k4*y0**3,0.]
        def rhs(y,v):
            rho,psi,py,_=v; pot,pot1=stable_potential(psi)
            rpy=np.sqrt(1+rho*rho*(py*py/12-pot/6))
            return [rpy,py,pot1-4*rpy*py/rho,1/rho]
        sol=solve_ivp(rhs,(y0,yb),initial,method='DOP853',rtol=rtol,atol=[1e-13,1e-21,1e-21,1e-13],dense_output=dense)
        if not sol.success:raise RuntimeError(sol.message)
        return sol
    def residual(x):
        sol=integrate(x);rho,psi,py,_=sol.y[:,-1];p=psi-1;pot,_=stable_potential(psi)
        rpy=np.sqrt(1+rho*rho*(py*py/12-pot/6))
        return [rpy/rho-sigma(p)/6,py+sigma1(p)/2]
    guess=[math.log10(.2207*tdet**1.8),8.2282+.9*math.log(.001/tdet)]
    found=root(residual,guess,tol=1e-9,options={'eps':1e-9})
    res=np.array(residual(found.x))
    if np.max(np.abs(res))>2e-10:raise RuntimeError('Shell root failed: '+repr(res))
    sol=integrate(found.x,True); end=sol.y[:,-1]; z_end=end[3]
    def sample(z):
        targets=np.asarray(z)+z_end
        if targets.min()<0:raise ValueError('Requested point before cone table')
        ix=np.clip(np.searchsorted(sol.y[3],targets)-1,0,len(sol.t)-2)
        ys=np.array([brentq(lambda y:sol.sol(y)[3]-zz,sol.t[i],sol.t[i+1],xtol=2e-14) for zz,i in zip(targets,ix)])
        vals=sol.sol(ys);rho,psi,py=vals[:3];p=psi-1;pot,_=stable_potential(psi)
        return dict(rho=rho,phi=p,psi=psi,phiz=rho*py,Hc=np.sqrt(1+rho*rho*(py*py/12-pot/6)),lnrho=np.log(rho))
    return sample,dict(tdet=tdet,dc=dc,psi_h=float(10**found.x[0]),yb=float(found.x[1]),rho_b=float(end[0]),phi_b=float(end[1]-1),root_residual=res.tolist(),rtol=rtol)

def hermite_ghosts(order=5):
    # Degree-five polynomial through values at x=0,-1,-2,-3,-4 and derivative at zero.
    V=np.zeros((6,6)); xx=-np.arange(5.)
    V[:5]=xx[:,None]**np.arange(6)[None,:];V[5,1]=1
    return np.array([np.arange(1.,4.)[i]**np.arange(6) for i in range(3)])@np.linalg.inv(V)

def grid(hmin,L,stretch=.05):
    alpha=hmin/stretch;N=math.ceil(math.log1p(L/stretch)/alpha)
    j=np.arange(N,-1,-1,dtype=float)
    z=-stretch*np.expm1(alpha*j);zp=hmin*np.exp(alpha*j);zpp=-alpha*zp
    return z,zp,zpp

def matrices(zp,zpp):
    n=len(zp);gh=hermite_ghosts();E=sparse.lil_matrix((n+4,n)); gv=np.zeros(n+4)
    E[2:n+2]=sparse.eye(n)
    for k in range(2):
        for j in range(5):E[n+2+k,n-1-j]=gh[k,j]
        gv[n+2+k]=gh[k,5]*zp[-1]
    E=E.tocsr();D1x=sum(c*E[k:k+n] for k,c in enumerate([1,-8,0,8,-1]))/12
    D2x=sum(c*E[k:k+n] for k,c in enumerate([-1,16,-30,16,-1]))/12
    g1=sum(c*gv[k:k+n] for k,c in enumerate([1,-8,0,8,-1]))/12
    g2=sum(c*gv[k:k+n] for k,c in enumerate([-1,16,-30,16,-1]))/12
    D1=sparse.diags(1/zp)@D1x
    D2=sparse.diags(1/zp**2)@D2x-sparse.diags(zpp/zp**3)@D1x
    return D1.tocsr(),D2.tocsr(),g1/zp,g2/zp**2-g1*zpp/zp**3,gh

class Solver:
    def __init__(self,tdet=.001,hmin=2e-4,L=6,stretch=.05,dc=0.,ko=0.):
        self.z,self.zp,self.zpp=grid(hmin,L,stretch);self.n=len(self.z)
        self.D1,self.D2,self.g1,self.g2,self.ghost=matrices(self.zp,self.zpp)
        samp,self.meta=shell_background(tdet);self.bg=samp(self.z)
        self.rho=self.bg['rho'];self.hc=self.bg['Hc'];self.phz=self.bg['phiz'];self.phi=self.bg['phi']
        self.rb=self.rho[-1];self.pb=self.phi[-1];self.rho2=self.rho**2
        self.tdet=tdet;self.ko=ko;self.hmin=hmin;self.L=L;self.stretch=stretch
        self.pot=derivatives(self.phi);self.s0=2*W(self.pb)+tdet*(1+C*self.pb)
        self.s10=2*(self.pb*self.pb-1)+tdet*C;self.s20=4*self.pb
        self.state0=np.zeros((6,self.n));self.seed_meta=None
        if dc:
            sample,self.seed_meta=shell_background(tdet,dc);seed=sample(self.z)
            self.state0[0]=seed['lnrho']-self.bg['lnrho'];self.state0[1]=self.state0[0]
            self.state0[2]=seed['psi']-self.bg['psi']
            # C2 quintic smoothstep taper, outside the near-shell causal window.
            x=np.clip((self.z+.99*L)/(.14*L),0,1);w=x**3*(10-15*x+6*x*x)
            self.state0[:3]*=w
        self.state0[:,:2]=0
    def bc(self,b,f,pb,pf):
        eb=np.exp(b);ds=self.s10*f+2*self.pb*f*f+(2/3)*f**3;ds1=self.s20*f+2*f*f
        ga=self.rb*(self.s0*np.expm1(b)+eb*ds)/6
        gf=-self.rb*(self.s10*np.expm1(b)+eb*ds1)/2
        gpa=self.rb*eb*((self.s0+ds)*pb+(self.s10+ds1)*pf)/6
        gpf=-self.rb*eb*((self.s10+ds1)*pb+4*(self.pb+f)*pf)/2
        return ga,gf,gpa,gpf
    def dissipation(self,v,g):
        gh=self.ghost[:,:5]@v[-1:-6:-1]+self.ghost[:,5]*self.zp[-1]*g
        ex=np.concatenate((np.zeros(3,dtype=v.dtype),v,gh))
        return self.ko*sum(c*ex[i:i+self.n] for i,c in enumerate([1,-6,15,-20,15,-6,1]))/(64*self.zp)
    def rhs(self,state):
        a,b,f,pa,pb,pf=state;ga,gf,gpa,gpf=self.bc(b[-1],f[-1],pb[-1],pf[-1])
        deriv1=self.D1@state[:3].T;deriv2=self.D2@state[:3].T
        deriv1+=self.g1[:,None]*np.array([ga,ga,gf]);deriv2+=self.g2[:,None]*np.array([ga,ga,gf])
        az,bz,fz=deriv1.T;azz,bzz,fzz=deriv2.T
        du=np.zeros_like(f);du1=np.zeros_like(f)
        for j in range(6,0,-1):du=(du+self.pot[j]/math.factorial(j))*f
        for j in range(5,0,-1):du1=(du1+self.pot[j+1]/math.factorial(j))*f
        eb=np.expm1(2*b);su=self.rho2*(eb*(self.pot[0]+du)+du);su1=self.rho2*(eb*(self.pot[1]+du1)+du1)
        kin=2*pa+pa*pa;grad=2*self.hc*az+az*az
        out=np.empty_like(state);out[:3]=state[3:]
        out[3]=azz-3*kin+3*grad+(2/3)*su
        out[4]=bzz+3*kin-3*grad-.5*pf*pf+self.phz*fz+.5*fz*fz-su/3
        out[5]=fzz-3*(1+pa)*pf+3*(self.hc*fz+az*self.phz+az*fz)-su1
        if self.ko:
            out[3]+=self.dissipation(pa,gpa);out[4]+=self.dissipation(pb,gpa);out[5]+=self.dissipation(pf,gpf)
        out[:,:2]=0
        return out
    def linear_matrix(self):
        n=self.n;zero=sparse.csr_matrix((n,n));I=sparse.eye(n,format='csr')
        # Spatial derivatives of each field, including the linear coupled junction.
        E=sparse.csr_matrix(([1.],([0],[n-1])),shape=(1,n))
        ga=[0,self.rb*self.s0/6,self.rb*self.s10/6];gf=[0,-self.rb*self.s10/2,-self.rb*self.s20/2]
        D=[]
        for K,g in [(self.D1,self.g1),(self.D2,self.g2)]:
            D.append([[((K if i==j else zero)+sparse.csr_matrix(g[:,None])@E*(ga if i<2 else gf)[j]) for j in range(3)] for i in range(3)])
        d1,d2=D;hc=sparse.diags(self.hc);pz=sparse.diags(self.phz);r2=sparse.diags(self.rho2)
        uv=[sparse.diags(v) for v in self.pot[:3]]
        rows=[[zero,zero,zero,I,zero,zero],[zero,zero,zero,zero,I,zero],[zero,zero,zero,zero,zero,I]]
        r=[]
        for j in range(3):r.append(d2[0][j]+6*hc@d1[0][j]+((4/3)*r2@uv[0] if j==1 else (2/3)*r2@uv[1] if j==2 else zero))
        rows.append(r+[-6*I,zero,zero])
        r=[]
        for j in range(3):r.append(d2[1][j]-6*hc@d1[0][j]+pz@d1[2][j]-((2/3)*r2@uv[0] if j==1 else (1/3)*r2@uv[1] if j==2 else zero))
        rows.append(r+[6*I,zero,zero])
        r=[]
        for j in range(3):r.append(d2[2][j]+3*hc@d1[2][j]+3*pz@d1[0][j]-(2*r2@uv[1] if j==1 else r2@uv[2] if j==2 else zero))
        rows.append(r+[zero,zero,-3*I])
        M=sparse.bmat(rows,format='lil')
        for block in range(6):
            for i in range(2):M[block*n+i,:]=0
        return M.tocsc()
    def diagnostics(self,state,t):
        a,b,f,pa,pb,pf=state;ga,gf,gpa,gpf=self.bc(b[-1],f[-1],pb[-1],pf[-1])
        dz=self.D1@state[:3].T+self.g1[:,None]*np.array([ga,ga,gf]);az,bz,fz=dz.T
        azz=self.D2@a+self.g2*ga;paz=self.D1@pa+self.g1*gpa
        # Endpoint pa_z must be a genuinely independent one-sided stencil, not inserted boundary data.
        one=np.array([25/12,-4,3,-4/3,1/4])/self.zp[-1]
        independent_paz=np.dot(one,pa[-1:-6:-1]);paz[-1]=independent_paz
        du=sum(self.pot[j]*f**j/math.factorial(j) for j in range(1,7));su=self.rho2*(np.expm1(2*b)*(self.pot[0]+du)+du)
        M=-3*paz-3*(1+pa)*(az-bz)+3*(self.hc+az)*pb-pf*(self.phz+fz)
        H=-2*su+12*pa+6*pa*pa+6*(1+pa)*pb-18*self.hc*az+6*self.hc*bz-12*az*az+6*az*bz-6*azz-pf*pf-2*self.phz*fz-fz*fz
        weight=(self.rho/self.rb)**3*np.exp(np.clip(3*(t+a),-700,700))
        cs=weight[None,:]*np.array([H+2*M,H-2*M])
        # Scale by a positive background magnitude; also retain unnormalized residuals.
        scale=1+6*self.hc*self.hc+self.phz*self.phz
        mask=(self.z>-.8*self.L)&(np.arange(self.n)>5)
        return dict(momentum_max=float(np.max(np.abs(M[mask])/scale[mask])),hamiltonian_max=float(np.max(np.abs(H[mask])/scale[mask])),
            boundary_velocity_slope_error=float(independent_paz-gpa),weighted_characteristic_max=float(np.max(np.abs(cs[:,mask])))),H,M,cs

def spectral_run(args):
    s=Solver(args.tdet,args.hmin,args.L,args.stretch);M=s.linear_matrix()
    rng=np.random.default_rng(4801);v=rng.normal(size=(6,s.n));v[:,:2]=0
    numeric=s.rhs(1j*1e-24*v).imag/1e-24
    jacerr=float(np.max(np.abs(M@v.ravel()-numeric.ravel()))/(1+np.max(np.abs(numeric))))
    if jacerr>2e-12:raise RuntimeError('Linear matrix/RHS mismatch '+str(jacerr))
    vals,vecs=eigs(M,k=6,sigma=args.shift,tol=5e-11,maxiter=3000)
    order=np.argsort(np.abs(vals-args.shift));vals=vals[order];vecs=vecs[:,order]
    modes=[]
    for k,ev in enumerate(vals):
        v=vecs[:,k];res=np.linalg.norm(M@v-ev*v)/np.linalg.norm(v)
        modes.append(dict(real=float(ev.real),imag=float(ev.imag),absolute_eigen_residual=float(res)))
    chosen=vecs[:,0].real.reshape(6,s.n);chosen/=chosen[2,-1]
    diag,hh,mm,cc=s.diagnostics(chosen*1e-7,0.)
    out=dict(status='FLOATING_SPECTRAL_CONTROL',background=s.meta,hmin=args.hmin,stretch=args.stretch,L=args.L,nodes=s.n,
        matrix_rhs_relative_error=jacerr,eigenvalues=modes,selected_diagnostics_at_phi_amplitude_1e_minus7=diag,
        scope='Discrete eigenvalues are convergence diagnostics, not proof of physical mode uniqueness or nonlinear fate.')
    np.savez_compressed(args.output+'.npz',mode=chosen,z=s.z,eigenvalue=vals[0],background=np.array([s.rho,s.phi,s.hc,s.phz]))
    Path(args.output+'.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

def evolve(args):
    s=Solver(args.tdet,args.hmin,args.L,args.stretch,dc=args.dc,ko=args.ko);v=s.state0
    if args.mode:
        with np.load(args.mode) as f:
            if np.max(np.abs(f['z']-s.z))>1e-12:raise ValueError('Mode grid mismatch')
            v=np.array(f['mode'])*args.amplitude
    dt=.4*s.zp.min();steps=math.ceil(args.tf/dt);dt=args.tf/steps
    records=[];snaps={};start=time.monotonic();tau=0.;every=max(1,round(.025/dt));reason='tf'
    for i in range(steps+1):
        t=i*dt
        if i%every==0 or i==steps:
            diag,H,M,cc=s.diagnostics(v,t);phi_b=s.pb+v[2,-1];hj=(1+v[3,-1])*np.exp(-v[1,-1])
            row=dict(time=t,H0tau=tau,phi_b=phi_b,scalar_deviation=v[2,-1],H_over_H0=hj,**diag);records.append(row)
            if len(records)%20==1:print(json.dumps({'progress':row,'elapsed':time.monotonic()-start}),flush=True)
            if len(records)%40==1 or i==steps:
                snaps['state_'+str(i)]=v.copy();snaps['constraints_'+str(i)]=np.array([H,M,*cc])
            if not np.isfinite(v).all() or np.max(np.abs(v[2]))>args.max_f:reason='amplitude_or_finite_guard';break
        if i==steps:break
        oldtau=np.exp(v[1,-1]);k1=s.rhs(v);k2=s.rhs(v+.5*dt*k1);k3=s.rhs(v+.5*dt*k2);k4=s.rhs(v+dt*k3)
        v+=dt*(k1+2*k2+2*k3+k4)/6;tau+=dt*.5*(oldtau+np.exp(v[1,-1]))
    result=dict(status='FLOATING_EVOLUTION_CONTROL',background=s.meta,seed_background=s.seed_meta,hmin=args.hmin,L=args.L,stretch=args.stretch,ko=args.ko,
        dc=args.dc,mode_seed=args.mode,amplitude=args.amplitude,t_end=records[-1]['time'],stop_reason=reason,runtime_seconds=time.monotonic()-start,records=records)
    np.savez_compressed(args.output+'.npz',z=s.z,final_state=v,background=np.array([s.rho,s.phi,s.hc,s.phz]),**snaps)
    Path(args.output+'.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:x for k,x in result.items() if k!='records'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['spectrum','evolve']);p.add_argument('--tdet',type=float,default=.001)
    p.add_argument('--hmin',type=float,default=2e-4);p.add_argument('--stretch',type=float,default=.05);p.add_argument('--L',type=float,default=6)
    p.add_argument('--shift',type=float,default=1.6572);p.add_argument('--output',required=True);p.add_argument('--tf',type=float,default=4)
    p.add_argument('--dc',type=float,default=1e-7);p.add_argument('--ko',type=float,default=0.);p.add_argument('--mode')
    p.add_argument('--amplitude',type=float,default=1e-7);p.add_argument('--max-f',type=float,default=.02)
    a=p.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    spectral_run(a) if a.operation=='spectrum' else evolve(a)
