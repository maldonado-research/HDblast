#!/usr/bin/env python3
"""Coherent ordinary first-field third jets from authenticated actual raw caps.

This adds coarse magnitude enclosures, not derivatives of interval selections.
The retained actual lower two-jets are returned unchanged. All assertions about
actual functions retain the same smooth background/master/gauge PDE premises.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path

GEOM='HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
GEOM_MANIFEST='3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
SELECTED={
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json':'d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e',
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/checks/check_master_fourth.py':'6962e8c9bdca416ad072029df1d84ae666467e11656f0078f3e5ba4033659ed0',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json':'94579087d39d9c849e9264d40bd1d71d64892a668fe3e2e08aade9ca4296b873',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json':'63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_actual_gauge_fourth.py':'9d88529967c6da55684779a67c828aa15062052ccd8307e188c7211963cf305b',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_original_bulk_fx.py':'30def8af98798c18a5f48d73294945da5522d5c3fb11269b1d8226cde70676ee',
}
KEYS=tuple((i,j) for i in range(5) for j in range(5-i))
THIRD_SUFFIXES=('yyy','yyt','ytt','ttt')


def need(ok,why):
    if not ok: raise ValueError(why)


def sha(raw): return hashlib.sha256(raw).hexdigest()


def authenticate(repo):
    source=Path(repo).resolve()/'research'
    cap=source/GEOM
    manifest=cap/'MANIFEST.json';raw=manifest.read_bytes()
    need(sha(raw)==GEOM_MANIFEST,'geometric manifest authentication')
    seen={'MANIFEST.json'};total=0
    for row in json.loads(raw)['files']:
        rel=Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'safe source member')
        seen.add(str(rel));p=cap/rel
        need(p.is_file() and not p.is_symlink(),'source member type')
        data=p.read_bytes()
        need(len(data)==row['bytes'] and sha(data)==row['sha256'],'source member authentication '+str(rel))
        total+=len(data)
    actual={str(p.relative_to(cap)) for p in cap.rglob('*') if p.is_file()}
    need(actual==seen and not any(p.is_symlink() for p in cap.rglob('*')),'closed geometric source tree')
    records={GEOM+'/MANIFEST.json':{'sha256':GEOM_MANIFEST,'members':len(seen)-1,'bytes_excluding_manifest':total}}
    for rel,pin in SELECTED.items():
        p=source/rel;data=p.read_bytes()
        need(p.is_file() and not p.is_symlink() and sha(data)==pin,'selected higher-jet source '+rel)
        records[rel]={'sha256':pin,'bytes':len(data)}
    return records


def physical_through_three(j,rho,a,ay,ell,mutation='none'):
    """Ordinary physical derivatives of the same smooth composite profile."""
    L=1/(ell*rho.val());av=a.val()
    third=-L**3*j.derivative(0,3)-3*av*L**2*j.derivative(0,2)+(ay-av**2)*L*j.derivative(0,1)
    if mutation=='third_radial_sign':third=L**3*j.derivative(0,3)-3*av*L**2*j.derivative(0,2)+(ay-av**2)*L*j.derivative(0,1)
    if mutation=='omit_middle_radial_term':third=-L**3*j.derivative(0,3)+(ay-av**2)*L*j.derivative(0,1)
    if mutation=='omit_background_radial_term':third=-L**3*j.derivative(0,3)-3*av*L**2*j.derivative(0,2)
    if mutation=='wrong_background_radial_sign':third=-L**3*j.derivative(0,3)-3*av*L**2*j.derivative(0,2)+(ay+av**2)*L*j.derivative(0,1)
    out={'':j.val(),'y':-L*j.derivative(0,1),'t':j.derivative(1,0)/ell,
         'yy':L**2*j.derivative(0,2)+av*L*j.derivative(0,1),
         'yt':-L*j.derivative(1,1)/ell,'tt':j.derivative(2,0)/ell**2,
         'yyy':third,
         'yyt':(L**2*j.derivative(1,2)+av*L*j.derivative(1,1))/ell,
         'ytt':-L*j.derivative(2,1)/ell**2,'ttt':j.derivative(3,0)/ell**3}
    if mutation=='omit_yyt_rho_term':out['yyt']=L**2*j.derivative(1,2)/ell
    if mutation=='wrong_ytt_sign':out['ytt']=L*j.derivative(2,1)/ell**2
    if mutation=='lose_third_time_ell':out['ttt']=j.derivative(3,0)/ell**2
    return out


def taylor_type(interval):
    class T:
        """Raw actual normalized jets through four; output is used only to three."""
        def __init__(self,c=None):self.c={k:interval.cv((c or {}).get(k,0)) for k in KEYS}
        @staticmethod
        def cv(v):return v if isinstance(v,T) else T({(0,0):v})
        def __add__(self,v):
            v=T.cv(v);return T({k:self.c[k]+v.c[k] for k in KEYS})
        __radd__=__add__
        def __neg__(self):return T({k:-v for k,v in self.c.items()})
        def __sub__(self,v):return self+-T.cv(v)
        def __rsub__(self,v):return T.cv(v)+-self
        def __mul__(self,v):
            v=T.cv(v);return T({(i,j):sum((self.c[a,b]*v.c[i-a,j-b] for a in range(i+1) for b in range(j+1)),interval(0)) for i,j in KEYS})
        __rmul__=__mul__
        def inv(self):
            out={(0,0):self.c[0,0].inv()}
            for n in range(1,5):
                for i in range(n+1):
                    j=n-i;out[i,j]=-out[0,0]*sum((self.c[a,b]*out[i-a,j-b] for a in range(i+1) for b in range(j+1) if a+b),interval(0))
            return T(out)
        def __truediv__(self,v):return self*T.cv(v).inv()
        def __rtruediv__(self,v):return T.cv(v)*self.inv()
        def __pow__(self,n):
            if n<0:return self.inv()**(-n)
            out=T.cv(1)
            for _ in range(n):out=out*self
            return out
        def dx(self):return T({(i,j):(i+1)*self.c.get((i+1,j),interval(0)) for i,j in KEYS})
        def dr(self):return T({(i,j):(j+1)*self.c.get((i,j+1),interval(0)) for i,j in KEYS})
        def derivative(self,i,j):return math.factorial(i)*math.factorial(j)*self.c[i,j]
        def val(self):return self.c[0,0]
    return T


class Adapter:
    """Retained lower actual enclosures plus uniform actual third-jet magnitudes."""
    def __init__(self,repo):
        self.repo=Path(repo).resolve();self.before=authenticate(self.repo)
        nr=self.repo/'research'/GEOM/'numeric'
        spec=importlib.util.spec_from_file_location('authenticated_original_actual_jets',nr/'actual_first_jets.py')
        old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
        self.old=old;self.core=old.core;self.I=old.I;self.ell=old.ELL;self.T=taylor_type(self.I)
        self.data,_,_=old.inputs()
        m4=json.loads((self.repo/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json').read_text())
        g4=json.loads((self.repo/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json').read_text())
        need(m4['status']=='PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS','actual master4 status')
        need(g4['status']=='PASS_CONDITIONAL_ACTUAL_FIRST_GAUGE_CONTINUOUS_FOURTH_JETS','actual gauge4 status')
        need(self.ell==F(1,10**6) and m4['domain']=='0<=x<=3, r>=0, x+r<=8; x=t/ell, r=z/ell, ell=1e-6'
             and g4['domain']=='0<=x<=3,r>=0,x+r<=8; x=t/ell,r=z/ell,ell=1e-6','same actual scaled triangle and units')
        need(g4['input_sha256']['master4']==SELECTED['HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json'],'gauge/master common dependency')
        for rel,pin in m4['input_sha256'].items():
            need(sha((self.repo/rel).read_bytes())==pin,'actual master input join '+rel)
        for role,rel in [('coefficients','GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'),('gauge_premise','GAUGE_WAVE_REGULARITY_PREMISE.md')]:
            need(sha((nr/'inputs'/rel).read_bytes())==g4['input_sha256'][role],'actual gauge input join '+role)
        master=json.loads((nr/'inputs/HIGHER_JETS_NORMAL.json').read_text())['exact']['interior']
        self.master_caps=[F(master['U']),F(master['first']),max(map(F,master['second'].values())),
                          max(map(F,master['third'].values())),max(map(F,m4['exact']['all_scaled_order_four_mixed_caps'].values()))]
        self.gauge_caps={tuple(map(int,k.split(','))):list(map(F,v)) for k,v in g4['actual_mixed_scaled_derivative_caps'].items()}
        need(set(self.gauge_caps)==set(KEYS) and all(len(v)==3 and min(v)>=0 for v in self.gauge_caps.values()),'all 45 actual gauge mixed caps')
        need(len(self.master_caps)==5 and min(self.master_caps)>0,'all actual master order caps')
        self.global_third,self.coefficient_boxes=self._global_third()
        self.readback()

    def raw_jets(self):
        I,T=self.I,self.T
        W=T({(i,j):I(-self.master_caps[i+j],self.master_caps[i+j])/(math.factorial(i)*math.factorial(j)) for i,j in KEYS})
        U=[T({(i,j):I(-self.gauge_caps[i,j][k],self.gauge_caps[i,j][k])/(math.factorial(i)*math.factorial(j)) for i,j in KEYS}) for k in range(3)]
        return W,U

    def _global_third(self,unused_fourth_sentinel=0):
        I,T,ell=self.I,self.T,self.ell
        rr,eta,pp,aa,hh,RR,BB,kk,_,_=self.core.background(I(0,8),self.data,'none')
        def radial(j):return T({**{(0,n):ell**n*j.c[n] for n in range(4)},(0,4):unused_fourth_sentinel})
        rho,p,a,R,B,k=map(radial,(rr,pp,aa,RR,BB,kk))
        ct,ht,vt=self.old.temporal(I(0,3))
        c,H,vartheta=(T({**t.c,(4,0):unused_fourth_sentinel}) for t in (ct,ht,vt))
        W,U=self.raw_jets()
        amps,scalar=self.old.frame_amplitudes(W,U,rho,p,a,R,B,k,c,H,vartheta)
        ay=-1/rho.val()**2-p.val()**2/3
        out={name+suffix:value for name,j in {**amps,'s':scalar}.items()
             for suffix,value in physical_through_three(j,rho,a,ay,ell).items() if suffix in THIRD_SUFFIXES}
        need(len(out)==32,'all 28 metric and four scalar third jets')
        coeff={name:{str(i)+','+str(j):v.c[i,j].row() for i,j in KEYS if i+j<=3}
               for name,v in [('rho',rho),('p',p),('a',a),('R',R),('B',B),('k',k),('c',c),('H',H),('vartheta',vartheta)]}
        return out,coeff

    def actual_jet(self,x,r,clip_causal=False):
        """Return all 70 metric and 10 scalar ordinary physical derivative boxes.

        Lower two-jets retain the existing continuous cell enclosures. Third
        derivatives are coarser uniform magnitude enclosures on the same actual
        physical triangle; they are not differentiated lower interval boxes.
        """
        self.readback()
        m,s,bg=self.old.actual_jet(x,r,clip_causal=clip_causal)
        m.update({k:v for k,v in self.global_third.items() if not k.startswith('s')})
        s.update({k:v for k,v in self.global_third.items() if k.startswith('s')})
        need(len(m)==70 and len(s)==10,'complete ordinary third-jet interface')
        self.readback()
        return m,s,bg

    def readback(self):
        need(authenticate(self.repo)==self.before,'actual jet source closure changed')


def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {str(k):encode(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [encode(x) for x in v]
    if hasattr(v,'row'):return v.row()
    return v
