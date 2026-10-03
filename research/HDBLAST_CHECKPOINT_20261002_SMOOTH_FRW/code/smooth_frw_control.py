"""Prospective prescribed-background FRW control; execution requires a protocol hash.

Physical complex modes are evolved from exact static incoming data.  Complete
positive-reference fourth-order stress and second-order field-square terms are
expanded as in the archived rational-jet verifier.  No physical backreaction is
solved.  See the registered protocol for the finite-action convention and gates.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import uuid

import numpy as np
import scipy
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.special import roots_legendre

LD = np.longdouble
ORDER = 5
R = 4.0
END = 1.25
CUTOFFS = (24.0, 48.0, 96.0)
OBS = ("rho", "p", "Q")


class Jet:
    """Ordinary Taylor coefficients through order five, evaluated in long double."""

    # NumPy-left arithmetic must dispatch to reflected Jet operators, not form
    # an object array containing separate Jet instances.
    __array_priority__ = 1000

    def __init__(self, coefficients):
        self.c = np.asarray(coefficients, dtype=LD)
        if self.c.shape[0] != ORDER + 1:
            raise ValueError("Wrong jet order")

    @classmethod
    def constant(cls, value):
        value = np.asarray(value, dtype=LD)
        c = np.zeros((ORDER + 1,) + value.shape, dtype=LD)
        c[0] = value
        return cls(c)

    @staticmethod
    def asjet(value):
        return value if isinstance(value, Jet) else Jet.constant(value)

    def __add__(self, other):
        other = self.asjet(other)
        shape = np.broadcast_shapes(self.c.shape[1:], other.c.shape[1:])
        a = self.c.reshape((ORDER + 1,) + (1,) * (len(shape)-self.c.ndim+1) + self.c.shape[1:])
        b = other.c.reshape((ORDER + 1,) + (1,) * (len(shape)-other.c.ndim+1) + other.c.shape[1:])
        return Jet(a+b)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.c)

    def __sub__(self, other):
        return self + -self.asjet(other)

    def __rsub__(self, other):
        return self.asjet(other) + -self

    def __mul__(self, other):
        other = self.asjet(other)
        return Jet(np.stack([sum(self.c[j]*other.c[n-j] for j in range(n+1))
                             for n in range(ORDER+1)]))

    __rmul__ = __mul__

    def inv(self):
        out = [1/self.c[0]]
        for n in range(1, ORDER+1):
            out.append(-sum(self.c[j]*out[n-j] for j in range(1,n+1))/self.c[0])
        return Jet(np.stack(out))

    def __truediv__(self, other):
        return self * self.asjet(other).inv()

    def __rtruediv__(self, other):
        return self.asjet(other)*self.inv()

    def __pow__(self, power):
        if not isinstance(power, int):
            raise TypeError("Jet powers must be integers")
        if power < 0:
            return self.inv()**(-power)
        out, base = Jet.constant(1), self
        while power:
            if power & 1:
                out = out*base
            base = base*base
            power //= 2
        return out

    def d(self):
        return Jet(np.concatenate((self.c[1:]*np.arange(1,ORDER+1).reshape(
            (ORDER,)+(1,)*(self.c.ndim-1)), np.zeros_like(self.c[:1])), axis=0))

    def exp(self):
        out = [np.exp(self.c[0])]
        for n in range(1, ORDER+1):
            out.append(sum(j*self.c[j]*out[n-j] for j in range(1,n+1))/n)
        return Jet(np.stack(out))

    def sqrt(self):
        out = [np.sqrt(self.c[0])]
        for n in range(1,ORDER+1):
            out.append((self.c[n]-sum(out[j]*out[n-j] for j in range(1,n)))/(2*out[0]))
        return Jet(np.stack(out))


def geometry(t, amplitude, static=False):
    """Fast independent analytic B,B',B'' for the mode ODE (T=1)."""
    if static or t <= 0:
        b = bp = bpp = 0.0
    elif t >= 1:
        b, bp, bpp = 1.0, 0.0, 0.0
    else:
        g = -1/t + 1/(1-t)
        small = np.exp(-abs(g))
        b = small/(1+small) if g < 0 else 1/(1+small)
        logistic_first = small/(1+small)**2
        logistic_second = logistic_first*((1-small)/(1+small) if g < 0 else (small-1)/(1+small))
        gp = 1/t**2 + 1/(1-t)**2
        gpp = -2/t**3 + 2/(1-t)**3
        bp = logistic_first*gp
        bpp = logistic_first*gpp + logistic_second*gp**2
    a = np.exp(amplitude*b)
    H = amplitude*bp
    app_over_a = amplitude*bpp + H*H
    x = R*(2*b-1)**2
    xp = 4*R*(2*b-1)*bp
    return a, H, app_over_a, x, xp


def background_jets(t, amplitude, static=False):
    t = np.asarray(t, dtype=LD)
    bc = np.zeros((ORDER+1, len(t)), dtype=LD)
    if not static:
        bc[0,t >= 1] = 1
        inside = (t > 0) & (t < 1)
        if np.any(inside):
            sc = np.zeros((ORDER+1, int(np.count_nonzero(inside))),dtype=LD)
            sc[0], sc[1] = t[inside], 1
            s = Jet(sc)
            g = 1/(1-s)-1/s
            left = g.c[0] <= 0
            bj = np.empty_like(sc)
            if np.any(left):
                eg = Jet(g.c[:,left]).exp()
                bj[:,left] = (eg/(1+eg)).c
            if np.any(~left):
                eg = (-Jet(g.c[:,~left])).exp()
                bj[:,~left] = (1/(1+eg)).c
            bc[:,inside] = bj
    b = Jet(bc)
    return (amplitude*b).exp(), R*(2*b-1)**2


def counterterms(k, t, amplitude, static=False):
    """Return fully expanded subtractions and analytic derivative of rho_sub.

    Every formula is the existing exact-rational verifier's formula.  Derivative
    jets are a subtraction consistency diagnostic, not the mode solver's source.
    """
    a0, x0 = background_jets(t, amplitude, static)
    a, x = Jet(a0.c[:,:,None]), Jet(x0.c[:,:,None])
    k = np.asarray(k,dtype=LD)[None,:]
    h = x-R
    w = (k*k+a*a*R).sqrt()
    H = a.d()/a
    D = a*a*h
    B = -k*k/3-a*a*R
    C = 1-B/w**2
    u = (D-a.d().d()/a)/(2*w)-w.d().d()/(4*w**2)+3*w.d()**2/(8*w**3)
    v = (-u**2/(2*w)-u.d().d()/(4*w**2)+w.d().d()*u/(4*w**3)
         +3*w.d()*u.d()/(4*w**3)-3*w.d()**2*u/(4*w**4))
    shared2 = H*w.d()/w**2+w.d()**2/(4*w**3)
    shared4 = (H*u.d()/w**2-2*H*w.d()*u/w**3+w.d()*u.d()/(2*w**3)
               -3*w.d()**2*u/(4*w**4))
    rho = (w/(2*a**4)+((D+H**2)/w+shared2)/(4*a**4)
           +(u**2/w-(D+H**2)*u/w**2+shared4)/(4*a**4))
    p = (k*k/(6*a**4*w)+(C*u+(H**2-D)/w+shared2)/(4*a**4)
         +(C*v+B*u**2/w**3-(H**2-D)*u/w**2+shared4)/(4*a**4))
    Q = 1/(2*a**2*w)-u/(2*a**2*w**2)
    return dict(rho=rho.c[0],p=p.c[0],Q=Q.c[0],rho_prime=rho.c[1],
                Q_prime=Q.c[1],Q_second=2*Q.c[2],
                a=a.c[0,:,0],H=H.c[0,:,0],H_prime=H.c[1,:,0],
                x=x.c[0,:,0],xp=x.c[1,:,0])


def finite_cutoff_anomaly(t,amplitude,cutoff,static=False):
    """Independent analytic finite-K trace moments from the action derivation.

    Cosmic H and its cosmic derivatives are evaluated from background jets.
    This integrates the subtraction trace remainder analytically, independently
    of the stress-mode quadrature and the WKB expressions above.
    """
    a,x = background_jets(t,amplitude,static)
    Hj = a.d()/a**2
    dj = Hj.d()/a
    ej = dj.d()/a
    fj = ej.d()/a
    zj = x.d()/a
    zzj = zj.d()/a
    H,d,e,f,z,zz = [j.c[0] for j in (Hj,dj,ej,fj,zj,zzj)]
    h=x.c[0]-R
    coefficients={
        5:R*(-18*H**2*d-6*H**2*h-9*H*e+5*H*z-3*d**2-4*d*h-f+3*h**2+zz)/16,
        7:-R**2*(-40*H**4-24*H**2*d+50*H**2*h+5*H*e+10*H*z+10*d*h+f)/32,
        9:7*R**3*(55*H**4+48*H**2*d+10*H**2*h+4*H*e+3*d**2)/64,
        11:-231*H**2*R**4*(3*H**2+d)/64,
        13:1155*H**4*R**5/256,
    }
    physical_cutoff=LD(cutoff)/a.c[0]
    v=physical_cutoff/np.sqrt(physical_cutoff**2+R)
    result=np.zeros(len(t),dtype=LD)
    for n,c in coefficients.items():
        m=(n-5)//2
        moment=LD(R)**(LD(3-n)/2)*sum((-1)**j*comb(m,j)*v**(2*j+3)/(2*j+3)
                                      for j in range(m+1))
        result+=c*moment
    return np.asarray(result/(2*LD(str(np.pi))**2),dtype=float)


def quadrature(cutoff, order):
    """Independent Gauss-Legendre panels of fixed width two in comoving k."""
    edges = np.arange(0.0,cutoff+1.0,2.0)
    if edges[-1] != cutoff:
        raise ValueError("Cutoff must be a positive multiple of two")
    nodes, weights = roots_legendre(order)
    k = ((edges[:-1,None]+edges[1:,None])/2 + nodes).ravel()
    weights = np.broadcast_to(weights,(len(edges)-1,order)).ravel().copy()
    return k,weights*k*k/(2*np.pi**2)


@dataclass(frozen=True)
class Settings:
    label: str
    order: int
    rtol: float
    atol: float
    phase_step: float


PRIMARY = Settings("primary",8,2e-11,2e-13,.2)
TIGHT = Settings("tight",8,2e-13,2e-15,.1)
QUADRATURE = Settings("quadrature",12,2e-13,2e-15,.1)
TAIL_STEP = Settings("tail_step",8,2e-13,2e-15,.05)


def evolve(k, t, amplitude, settings, static=False):
    """Physical complex u,u' plus independently integrated bare energy work."""
    k = np.asarray(k,dtype=float)
    omega = np.sqrt(k*k+R)
    u0 = 1/np.sqrt(2*omega)+0j
    v0 = -1j*omega*u0
    n = len(k)
    initial = np.concatenate((u0,v0,np.zeros(n,dtype=complex)))

    def rhs(eta,state):
        a,H,app_a,x,xp = geometry(eta,amplitude,static)
        u,v = state[:n],state[n:2*n]
        absu2 = abs(u)**2
        shifted_v2 = abs(v-H*u)**2
        work = H*(-shifted_v2+(k*k+2*a*a*x)*absu2)+a*a*xp*absu2/2
        return np.concatenate((v,-(k*k+a*a*x-app_a)*u,work.astype(complex)))

    high_frequency = np.sqrt(float(np.max(k))**2+np.exp(2*amplitude)*R)
    max_step = min(1/80,settings.phase_step/high_frequency)
    sol = solve_ivp(rhs,(float(t[0]),float(t[-1])),initial,method="DOP853",
                    t_eval=t,rtol=settings.rtol,atol=settings.atol,max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    if not np.all(np.isfinite(sol.y)):
        raise FloatingPointError("Nonfinite physical modes")
    return sol.y[:n].T,sol.y[n:2*n].T,sol.y[2*n:].real.T,dict(
        nfev=sol.nfev,rtol=settings.rtol,atol=settings.atol,max_step=max_step,
        solver="complex scipy.integrate.solve_ivp DOP853")


def integrate_modes(values,weights,mask=None):
    if mask is None:
        mask = np.ones(len(weights),dtype=bool)
    return np.sum(values[:,mask]*np.asarray(weights[mask],dtype=LD)[None,:],
                  axis=1,dtype=LD).astype(float)


def finite_difference_five(values,spacing,derivative=1):
    out = np.empty_like(values)
    if derivative==1:
        out[2:-2] = (values[:-4]-8*values[1:-3]+8*values[3:-1]-values[4:])/(12*spacing)
    elif derivative==2:
        out[2:-2] = (-values[:-4]+16*values[1:-3]-30*values[2:-2]+16*values[3:-1]-values[4:])/(12*spacing**2)
    else:
        raise ValueError("Only first and second derivatives implemented")
    for index in (0,1,len(values)-2,len(values)-1):
        start = max(0,min(index-2,len(values)-5))
        offsets = np.arange(start,start+5)-index
        matrix = np.vstack([offsets.astype(float)**j for j in range(5)])
        rhs = np.zeros(5)
        rhs[derivative]=1 if derivative==1 else 2
        out[index] = np.linalg.solve(matrix,rhs) @ values[start:start+5]/spacing**derivative
    return out


def finite_json(path,payload):
    path.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")


def difference_gate(first,second,atol=2e-5,rtol=.005):
    difference = float(np.max(abs(np.asarray(first)-np.asarray(second))))
    scale = float(max(np.max(abs(first)),np.max(abs(second))))
    limit = atol+rtol*scale
    return dict(absolute_difference=difference,signal_scale=scale,
                normalized_difference=difference/max(scale,atol),limit=limit,
                passed=bool(difference<=limit))


def exchange_diagnostics(t,integrated,a,H,xp):
    rho,p,Q = [np.asarray(integrated[key]) for key in OBS]
    conformal_energy = a**4*rho
    geometry_work = a**4*H*(rho-3*p)
    source_work = a**4*xp*Q/2
    total_work = geometry_work+source_work
    ledger = cumulative_simpson(total_work,x=t,initial=0)
    ledger_residual = conformal_energy-conformal_energy[0]-ledger
    ledger_scale = float(max(np.max(abs(conformal_energy-conformal_energy[0])),
                             np.max(abs(ledger))))
    ledger_limit = 1e-5+.005*ledger_scale
    derivative = finite_difference_five(rho,t[1]-t[0])
    expansion = 3*H*(rho+p)
    source = xp*Q/2
    local = derivative+expansion-source
    local_scale = float(max(np.max(abs(derivative)),np.max(abs(expansion)),np.max(abs(source))))
    local_limit = 2e-4+.01*local_scale
    omission = conformal_energy-conformal_energy[0]-cumulative_simpson(geometry_work,x=t,initial=0)
    return dict(
        ledger_residual_max=float(np.max(abs(ledger_residual))),ledger_scale=ledger_scale,
        ledger_limit=ledger_limit,ledger_passed=bool(np.max(abs(ledger_residual))<=ledger_limit),
        local_residual_max=float(np.max(abs(local))),local_scale=local_scale,
        local_limit=local_limit,local_passed=bool(np.max(abs(local))<=local_limit),
        source_omission_ledger_residual_max=float(np.max(abs(omission))),
        source_work_cumulative_max=float(np.max(abs(cumulative_simpson(source_work,x=t,initial=0)))),
        time_nodes=len(t)),dict(conformal_energy=conformal_energy,ledger=ledger,
        geometry_work=geometry_work,source_work=source_work,local_exchange_residual=local,
        ledger_residual=ledger_residual)


def trace_diagnostics(t,integrated,a,H,x,anomaly):
    trace=-integrated["rho"]+3*integrated["p"]
    direct_kinetic=(integrated["Q_second"]+2*H*integrated["Q_prime"])/(2*a*a)
    direct=direct_kinetic-x*integrated["Q"]+anomaly
    direct_gate=difference_gate(trace,direct,atol=2e-6,rtol=5e-6)
    sampled=[]
    arrays={"trace":trace,"trace_direct":direct,"trace_anomaly":anomaly}
    for stride in (1,2):
        Qt=integrated["Q"][::stride]
        spacing=(t[1]-t[0])*stride
        Qp=finite_difference_five(Qt,spacing)
        Qpp=finite_difference_five(Qt,spacing,derivative=2)
        rhs=(Qpp+2*H[::stride]*Qp)/(2*a[::stride]**2)-x[::stride]*Qt+anomaly[::stride]
        gate=difference_gate(trace[::stride],rhs,atol=2e-4,rtol=.01)
        sampled.append(dict(time_nodes=len(Qt),**gate))
        arrays[f"trace_sampled_{len(Qt)}"]=rhs
    return dict(direct_derivatives=direct_gate,sampled_derivatives=sampled,
                direct_scope="Physical equation Q derivatives plus subtraction jets, compared with independently integrated analytic finite-K trace moments."),arrays


def run_case(output,amplitude,settings,cutoff=96.0,static=False):
    label = f"A{amplitude:g}_{settings.label}_K{cutoff:g}"+("_static" if static else "")
    directory = output/label
    directory.mkdir()
    finite_json(directory/"started.json",dict(label=label,amplitude=amplitude,cutoff=cutoff,
        quadrature_order=settings.order,static=static,created=datetime.now(timezone.utc).isoformat()))
    t = np.linspace(0,END,801)
    k,weights = quadrature(cutoff,settings.order)
    u,v,ode_work,solver = evolve(k,t,amplitude,settings,static)
    # Preserve completed physical evolution even if a later subtraction or gate fails.
    np.savez_compressed(directory/"physical_modes.npz",t=t,k=k,weights=weights,u=u,v=v,
                        bare_conformal_work=ode_work)
    rho,p,Q,ct_rho,Q_prime,Q_second = [np.empty((len(t),len(k)),dtype=LD) for _ in range(6)]
    a,H,x,xp = [np.empty(len(t),dtype=float) for _ in range(4)]
    ct_exchange_max = 0.0
    for start in range(0,len(t),16):
        sl = slice(start,min(start+16,len(t)))
        ct = counterterms(k,t[sl],amplitude,static)
        aj,Hj,xj = ct["a"][:,None],ct["H"][:,None],ct["x"][:,None]
        uj,vj = u[sl].astype(np.clongdouble),v[sl].astype(np.clongdouble)
        absu2,shiftedv2 = abs(uj)**2,abs(vj-Hj*uj)**2
        rho[sl] = (shiftedv2+(np.asarray(k,dtype=LD)[None,:]**2+aj**2*xj)*absu2)/(2*aj**4)-ct["rho"]
        p[sl] = (shiftedv2-(np.asarray(k,dtype=LD)[None,:]**2/3+aj**2*xj)*absu2)/(2*aj**4)-ct["p"]
        Q[sl] = absu2/aj**2-ct["Q"]
        uvreal=np.real(vj*np.conj(uj))
        Hprime=ct["H_prime"][:,None]
        physical_vprime=-(np.asarray(k,dtype=LD)[None,:]**2+aj**2*xj-Hprime-Hj**2)*uj
        Q_prime[sl]=(2*uvreal-2*Hj*absu2)/aj**2-ct["Q_prime"]
        Q_second[sl]=(2*abs(vj)**2+2*np.real(physical_vprime*np.conj(uj))
                      -8*Hj*uvreal+(4*Hj**2-2*Hprime)*absu2)/aj**2-ct["Q_second"]
        ct_rho[sl] = ct["rho"]
        for dest,key in ((a,"a"),(H,"H"),(x,"x"),(xp,"xp")):
            dest[sl] = ct[key]
        ct_exchange = ct["rho_prime"]+3*ct["H"][:,None]*(ct["rho"]+ct["p"])-ct["xp"][:,None]*ct["Q"]/2
        ct_exchange_max = max(ct_exchange_max,float(np.max(abs(ct_exchange))))
    wronskian = (u*np.conj(v)-v*np.conj(u))/1j
    wronskian_error = float(np.max(abs(wronskian-1)))
    wr_limit = 2e-9 if settings.label == "primary" else 2e-10
    modes = dict(rho=rho,p=p,Q=Q)
    arrays = dict(t=t,k=k,weights=weights,a=a,H=H,x=x,xp=xp,
                  wronskian_error=np.abs(wronskian-1).astype(float))
    diagnostics = dict(label=label,amplitude=amplitude,static=static,cutoff=cutoff,mode_count=len(k),
        solver=solver,wronskian_max=wronskian_error,wronskian_limit=wr_limit,
        wronskian_passed=wronskian_error<=wr_limit,counterterm_exchange_pointwise_max=ct_exchange_max,
        finite=True,cutoffs={})
    integrated_by_cutoff = {}
    for K in tuple(c for c in CUTOFFS if c<=cutoff) + ((192.0,) if cutoff==192 else ()):
        mask = k < K
        integrated = {key:integrate_modes(values,weights,mask) for key,values in modes.items()}
        integrated["Q_prime"]=integrate_modes(Q_prime,weights,mask)
        integrated["Q_second"]=integrate_modes(Q_second,weights,mask)
        integrated_by_cutoff[str(int(K))] = integrated
        exchange,exchange_arrays = exchange_diagnostics(t,integrated,a,H,xp)
        coarse, _ = exchange_diagnostics(t[::2],{key:value[::2] for key,value in integrated.items()},
                                         a[::2],H[::2],xp[::2])
        ec = np.asarray(a,dtype=LD)[:,None]**4*ct_rho
        ode_ren_ledger = integrate_modes(np.asarray(ode_work,dtype=LD)-(ec-ec[0]),weights,mask)
        energy = a**4*integrated["rho"]
        ode_residual = energy-energy[0]-ode_ren_ledger
        ode_scale = float(max(np.max(abs(energy-energy[0])),np.max(abs(ode_ren_ledger))))
        ode_limit = 1e-5+.005*ode_scale
        # Future static decomposition retains coherent p and Q separately.
        future = t>=1 if not static else np.ones(len(t),dtype=bool)
        af = a[-1]
        wf = np.sqrt(k*k+af*af*R)
        beta = (wf[None,:]*u[future]-1j*v[future])/np.sqrt(2*wf[None,:])
        occupation = abs(beta)**2
        occ_rho = integrate_modes(occupation*wf[None,:]/af**4,weights,mask)
        occ_p = integrate_modes(occupation*k[None,:]**2/(3*af**4*wf[None,:]),weights,mask)
        occ_Q = integrate_modes(occupation/(af**2*wf[None,:]),weights,mask)
        spectral_gate = difference_gate(integrated["rho"][future],occ_rho,atol=2e-5,rtol=.005)
        trace,trace_arrays=trace_diagnostics(t,integrated,a,H,x,finite_cutoff_anomaly(t,amplitude,K,static))
        diagnostics["cutoffs"][str(int(K))] = dict(exchange_fine=exchange,exchange_coarse=coarse,
            trace=trace,
            ode_bare_work_residual_max=float(np.max(abs(ode_residual))),
            ode_bare_work_limit=ode_limit,ode_bare_work_passed=bool(np.max(abs(ode_residual))<=ode_limit),
            ode_bare_work_scope="Independent accumulated bare work versus energy; subtraction cancels algebraically here. Renormalized Simpson ledger separately tests subtraction/current pairing.",
            future_spectral_energy=spectral_gate,
            future_coherent_pressure_max=float(np.max(abs(integrated["p"][future]-occ_p))),
            future_coherent_Q_max=float(np.max(abs(integrated["Q"][future]-occ_Q))),
            maxima={key:float(np.max(abs(integrated[key]))) for key in OBS},
            crossing={key:float(integrated[key][320]) for key in OBS},
            endpoint={key:float(integrated[key][-1]) for key in OBS})
        for key,value in integrated.items():
            arrays[f"K{int(K)}_{key}"] = value
        for key,value in exchange_arrays.items():
            arrays[f"K{int(K)}_{key}"] = value
        for key,value in trace_arrays.items():
            arrays[f"K{int(K)}_{key}"]=value
        arrays[f"K{int(K)}_ode_ren_ledger"] = ode_ren_ledger
        arrays[f"K{int(K)}_future_occ_rho"] = occ_rho
        arrays[f"K{int(K)}_future_occ_p"] = occ_p
        arrays[f"K{int(K)}_future_occ_Q"] = occ_Q
    for key,value in modes.items():
        arrays[f"mode_{key}"] = np.asarray(value,dtype=float)
    if not all(np.all(np.isfinite(value)) for value in arrays.values()):
        raise FloatingPointError("Nonfinite result; failed directory preserved")
    np.savez_compressed(directory/"arrays.npz",**arrays)
    finite_json(directory/"diagnostics.json",diagnostics)
    print(json.dumps(dict(completed=label,wronskian_max=wronskian_error,
         maxima=diagnostics["cutoffs"][str(int(cutoff))]["maxima"]),allow_nan=False),flush=True)
    return dict(diagnostics=diagnostics,integrated=integrated_by_cutoff,t=t,a=a,H=H,xp=xp)


def compare_observables(first,second):
    return {key:difference_gate(first[key],second[key]) for key in OBS}


def all_passed(gates):
    return all(value["passed"] for value in gates.values())


def execute(output,protocol,protocol_hash):
    actual = hashlib.sha256(protocol.read_bytes()).hexdigest()
    if actual != protocol_hash:
        raise ValueError("Protocol hash mismatch; refuse execution")
    output.mkdir(parents=True,exist_ok=False)
    provenance = dict(protocol=str(protocol),protocol_sha256=actual,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        longdouble_epsilon=float(np.finfo(LD).eps),created=datetime.now(timezone.utc).isoformat(),
        scope="Prescribed C-infinity background control only; no shell change, backreaction, decay, or thermal radiation.",
        finite_action="Positive-reference subtraction baseline; interpretation subject to registered finite-action matching document.")
    finite_json(output/"provenance.json",provenance)
    final = dict(provenance=provenance,amplitudes={},status="IN_PROGRESS")
    finite_json(output/"summary.json",final)
    for amplitude in (0.0,.2):
        runs = {setting.label:run_case(output,amplitude,setting) for setting in (PRIMARY,TIGHT,QUADRATURE)}
        primary,tight,quad = [runs[label] for label in ("primary","tight","quadrature")]
        cutoff_gates = compare_observables(tight["integrated"]["48"],tight["integrated"]["96"])
        coarse_cutoff = compare_observables(tight["integrated"]["24"],tight["integrated"]["48"])
        solver_gates = compare_observables(primary["integrated"]["96"],tight["integrated"]["96"])
        quad_gates = compare_observables(tight["integrated"]["96"],quad["integrated"]["96"])
        selected = quad
        contingency = None
        if not all_passed(cutoff_gates):
            higher = run_case(output,amplitude,TIGHT,cutoff=192.0)
            higher_quad = run_case(output,amplitude,QUADRATURE,cutoff=192.0)
            higher_step = run_case(output,amplitude,TAIL_STEP,cutoff=192.0)
            cutoff_gates = compare_observables(higher["integrated"]["96"],higher["integrated"]["192"])
            quad_gates = compare_observables(higher["integrated"]["192"],higher_quad["integrated"]["192"])
            higher_solver_gates = compare_observables(higher["integrated"]["192"],higher_step["integrated"]["192"])
            continuity = compare_observables(tight["integrated"]["96"],higher["integrated"]["96"])
            selected = higher_quad
            contingency = dict(trigger="K48-to-K96 cutoff gate failed",continuity_at_K96=continuity,
                               cutoff_K96_to_K192=cutoff_gates,quadrature_at_K192=quad_gates,
                               solver_step_refinement_at_K192=higher_solver_gates,
                               all_extended_wronskians_passed=all(run["diagnostics"]["wronskian_passed"]
                                   for run in (higher,higher_quad,higher_step)))
        key = str(int(selected["diagnostics"]["cutoff"]))
        chosen = selected["diagnostics"]["cutoffs"][key]
        t=selected["t"]
        vals=selected["integrated"][key]
        fine,fine_arrays=exchange_diagnostics(t,vals,selected["a"],selected["H"],selected["xp"])
        coarse,coarse_arrays=exchange_diagnostics(t[::2],{name:value[::2] for name,value in vals.items()},
                                                 selected["a"][::2],selected["H"][::2],selected["xp"][::2])
        time_gate=difference_gate(fine_arrays["ledger"][::2],coarse_arrays["ledger"],atol=1e-5,rtol=.005)
        checks = [all_passed(cutoff_gates),all_passed(solver_gates),all_passed(quad_gates),time_gate["passed"],
                  chosen["exchange_fine"]["ledger_passed"],chosen["exchange_fine"]["local_passed"],
                  chosen["ode_bare_work_passed"],chosen["future_spectral_energy"]["passed"],
                  chosen["trace"]["direct_derivatives"]["passed"],
                  chosen["trace"]["sampled_derivatives"][0]["passed"],
                  all(run["diagnostics"]["wronskian_passed"] for run in runs.values()),
                  selected["diagnostics"]["wronskian_passed"]]
        if contingency is not None:
            checks.append(all_passed(contingency["continuity_at_K96"]))
            checks.append(all_passed(contingency["solver_step_refinement_at_K192"]))
            checks.append(contingency["all_extended_wronskians_passed"])
        initial_fine_cutoff=compare_observables(tight["integrated"]["48"],tight["integrated"]["96"])
        ratios={name:(initial_fine_cutoff[name]["absolute_difference"]/coarse_cutoff[name]["absolute_difference"]
                 if coarse_cutoff[name]["absolute_difference"]>1e-30 else None) for name in OBS}
        final["amplitudes"][str(amplitude)] = dict(status="PASS" if all(checks) else "FAIL",
            selected=selected["diagnostics"]["label"],initial_cutoff_K24_to_K48=coarse_cutoff,
            initial_cutoff_K48_to_K96=initial_fine_cutoff,
            cutoff_increment_shrink_ratios=ratios,
            final_cutoff_gate=cutoff_gates,solver_refinement=solver_gates,quadrature_refinement=quad_gates,
            time_refinement=time_gate,exchange_fine=fine,exchange_coarse=coarse,
            contingency=contingency,selected_diagnostics=chosen)
        finite_json(output/"summary.json",final)
    static=run_case(output,0.0,TIGHT,static=True)
    static_max=static["diagnostics"]["cutoffs"]["96"]["maxima"]
    static_ok=all(value<=2e-5 for value in static_max.values()) and static["diagnostics"]["wronskian_passed"]
    final["static_negative_control"]=dict(passed=static_ok,absolute_limit=2e-5,maxima=static_max,
        interpretation="Exact physical and subtracted rho,p,Q vanish; finite residual diagnoses solver/cancellation floor.")
    passed=all(value["status"]=="PASS" for value in final["amplitudes"].values()) and static_ok
    final["status"]="PASS" if passed else "FAIL"
    finite_json(output/"summary.json",final)
    print("HDBLAST_SMOOTH_FRW_JSON="+json.dumps(final,allow_nan=False),flush=True)
    return 0 if passed else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute",action="store_true",required=True,
                        help="Run only after prospective registration commit and authorization")
    parser.add_argument("--protocol",type=Path,required=True)
    parser.add_argument("--protocol-sha256",required=True)
    parser.add_argument("--output",type=Path,default=None)
    args=parser.parse_args()
    output=args.output or Path("runs")/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
    try:
        return execute(output,args.protocol,args.protocol_sha256)
    except Exception as exc:
        if output.exists():
            finite_json(output/"failure.json",dict(exception=type(exc).__name__,message=str(exc),
                                                  failed_at=datetime.now(timezone.utc).isoformat()))
        raise


if __name__ == "__main__":
    raise SystemExit(main())
