#!/usr/bin/env python3
"""Independent exact first natural-field binding; no scientific imports.

The tensor calculation uses flat dS slicing and arbitrary scalar coordinate
jets.  A separate global-slicing calculation binds the actual l=2 master.
All arithmetic and decisions are rational; no validation depends on assert.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


class Failure(Exception):
    pass


def need(ok, message):
    if not ok:
        raise Failure(message)


class P:
    def __init__(self, terms=None):
        self.d = {k:F(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def cv(v):
        return v if isinstance(v,P) else P({():F(v)})
    @staticmethod
    def var(name):
        return P({((name,1),):F(1)})
    def __add__(self,b):
        b=P.cv(b); d=dict(self.d)
        for k,v in b.d.items():
            d[k]=d.get(k,F(0))+v
            if not d[k]: del d[k]
        return P(d)
    __radd__=__add__
    def __neg__(self): return P({k:-v for k,v in self.d.items()})
    def __sub__(self,b): return self+-P.cv(b)
    def __rsub__(self,b): return P.cv(b)+-self
    def __mul__(self,b):
        b=P.cv(b); d={}
        for k,v in self.d.items():
            for l,w in b.d.items():
                powers=dict(k)
                for n,e in l: powers[n]=powers.get(n,0)+e
                key=tuple(sorted((n,e) for n,e in powers.items() if e))
                d[key]=d.get(key,F(0))+v*w
        return P(d)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            need(len(self.d)==1,'inverse must be a monomial')
            (key,value),=self.d.items()
            return P({tuple((x,e*n) for x,e in key):value**n})
        out=P.cv(1)
        for _ in range(n): out=out*self
        return out
    def __truediv__(self,b): return self*P.cv(b)**-1
    def __rtruediv__(self,b): return P.cv(b)*self**-1
    def derivative(self,rule):
        out=P()
        for key,value in self.d.items():
            for name,e in key:
                d=rule(name)
                if not d.d: continue
                powers=dict(key);powers[name]-=1
                base=P({tuple(sorted((n,q) for n,q in powers.items() if q)):value*e})
                out+=base*d
        return out
    def substitute(self,name,value):
        out=P()
        for key,co in self.d.items():
            powers=dict(key);n=powers.pop(name,0)
            out+=P({tuple(sorted(powers.items())):co})*value**n
        return out
    def serial(self):
        return [{'coefficient':str(v),'powers':dict(k)} for k,v in sorted(self.d.items())]


rho,a,p,V,V1,V2,V3,s,alpha,theta,H,ell = map(P.var,
    ('rho','a','p','V','V1','V2','V3','s','alpha','theta','H','ell'))
ZERO=P()


def psi(*indices):
    return P.var('Psi:'+','.join(map(str,sorted(indices))))


def w(i=0,j=0): return P.var('W:%d,%d'%(i,j))


def background_rule(name,axis,phase,mutation):
    if axis==0:
        rules={'rho':a*rho,'a':-rho**-2-p*p/3,
               'p':V1-4*a*p,'V':p*V1,'V1':p*V2,'V2':p*V3,
               'alpha':alpha*(V1/p-F(5,2)*a)}
        if mutation=='alpha_sign': rules['alpha']=-rules['alpha']
        if name in rules: return rules[name]
    if axis==1:
        if phase=='tensor' and name=='s': return s
        if phase=='separated':
            if name=='H': return 1-H*H
            if name=='theta':
                return ZERO if mutation=='theta_frozen' else -F(3,2)*H*theta
    return ZERO


def D(expr,axis,phase='tensor',mutation='none'):
    def rule(name):
        if name.startswith('Psi:'):
            old=name[4:]
            inds=[] if not old else list(map(int,old.split(',')))
            return psi(*(inds+[axis]))
        if name.startswith('W:'):
            i,j=map(int,name[2:].split(','))
            if axis==0:
                sign=1 if mutation=='radial_sign' else -1
                return sign*w(i,j+1)/(ell*rho)
            if axis==1:return w(i+1,j)/ell
            return ZERO
        return background_rule(name,axis,phase,mutation)
    return P.cv(expr).derivative(rule)


def hamiltonian(expr,mutation='none'):
    return expr if mutation=='no_hamiltonian' else expr.substitute('V',6/rho**2+p*p/2-6*a*a)


CHECKS=[]


def equal(name,left,right=0,ham=False,mutation='none'):
    remainder=P.cv(left)-right
    if ham: remainder=hamiltonian(remainder,mutation)
    if remainder.d:
        raise Failure(name+' | '+json.dumps(remainder.serial(),sort_keys=True,separators=(',',':')))
    CHECKS.append(name)


def tensor_calculation(mutation):
    """Full coordinate connection in dy²+rho²(-dt²+exp(2t)dx_i²)."""
    n=5;g=[P.cv(1),-rho*rho]+[rho*rho*s*s]*3
    gi=[1/v for v in g]; q=psi()
    # Build the background connection from metric derivatives rather than
    # copying a table of perturbation Ricci coefficients.
    G={}
    for up in range(n):
        for b in range(n):
            for c in range(n):
                v=gi[up]*(
                    (D(g[up],b) if up==c else 0)
                    +(D(g[up],c) if up==b else 0)
                    -(D(g[b],up) if b==c else 0))/2
                if v.d:G[up,b,c]=v
    def gam(up,b,c): return G.get((up,b,c),ZERO)
    hv=[-4*q]+[2*q*g[i] for i in range(1,n)]
    def h(b,c):return hv[b] if b==c else ZERO
    def dh(d,b,c):
        return D(h(b,c),d)-sum((gam(k,d,b)*h(k,c)+gam(k,d,c)*h(b,k) for k in range(n)),ZERO)
    C={}
    for up in range(n):
        for b in range(n):
            for c in range(n):
                v=gi[up]*(dh(b,c,up)+dh(c,b,up)-dh(up,b,c))/2
                if v.d:C[up,b,c]=v
    def dc(up,b,c):return C.get((up,b,c),ZERO)
    trace=[sum((dc(c,b,c) for c in range(n)),ZERO) for b in range(n)]
    ric={}
    for b in range(n):
        for c in range(b,n):
            v=sum((D(dc(k,b,c),k) for k in range(n)),ZERO)-D(trace[b],c)
            # nabla_k C^k_bc - nabla_c C^k_bk, with every slot.
            for j in range(n):
                for k in range(n):
                    v+=gam(k,k,j)*dc(j,b,c)-gam(j,k,b)*dc(k,j,c)-gam(j,k,c)*dc(k,b,j)
                v+=gam(j,c,b)*trace[j]
            ric[b,c]=v
    box=lambda z:-D(D(z,1),1)-3*D(z,1)+sum((D(D(z,i),i)/s**2 for i in range(2,n)),ZERO)
    expected_yy=-4*D(D(q,0),0)-16*a*D(q,0)+2*box(q)/rho**2
    expected_trace=-box(q)-rho*rho*(D(D(q,0),0)+10*a*D(q,0)+6*(D(a,0)+4*a*a)*q)
    for (b,c),value in ric.items():
        target=expected_yy if b==c==0 else (
            -3*D(D(q,0)+2*a*q,c) if b==0 else (
                g[b]/rho**2*expected_trace if b==c else ZERO))
        equal('arbitrary_Ricci_%d%d'%(b,c),value,target)
    sign=1 if mutation=='chi_sign' else -1
    Pmaster=rho*rho*q
    chi=sign*3*D(Pmaster,0)/(rho*rho*p)
    M=D(D(Pmaster,0),0)-(2*a+2*D(p,0)/p)*D(Pmaster,0)+(box(Pmaster)+4*Pmaster)/rho**2-F(2,3)*p*p*Pmaster
    for (b,c),value in ric.items():
        v=value-(p*D(chi,c) if b==0 else 0)-(p*D(chi,b) if c==0 else 0)
        if b==c:
            if mutation!='drop_Vprime':v-=F(2,3)*V1*chi*g[b]
            if mutation!='drop_Vh':v-=F(2,3)*V*hv[b]
        target=2*M/rho**2 if b==c==0 else (-g[b]/rho**2*M if b==c else ZERO)
        if b==c==0 and mutation in ('none','no_hamiltonian'):
            equal('explicit_background_Hamiltonian_witness',v-target,
                  16*q*(a*a-rho**-2-p*p/12+V/6))
        equal('natural_E_%d%d'%(b,c),v,target,ham=True,mutation=mutation)
    # Direct scalar variation including the nonharmonic divergence contraction.
    scalar_box=sum((gi[b]*(D(D(chi,b),b)-sum((gam(k,b,b)*D(chi,k) for k in range(n)),ZERO)) for b in range(n)),ZERO)
    Hess=[D(p,0)]+[a*p*g[b] for b in range(1,n)]
    scalar=scalar_box-sum((gi[b]*gi[b]*hv[b]*Hess[b] for b in range(n)),ZERO)-V2*chi
    divergence=sum((gi[b]*dc(0,b,b) for b in range(n)),ZERO)
    if mutation!='drop_DC':scalar-=divergence*p
    retained= D(D(chi,0),0)+4*a*D(chi,0)+box(chi)/rho**2-V2*chi+6*p*D(q,0)+4*V1*q
    equal('direct_scalar_variation',scalar,retained)
    equal('natural_S_master_derivative',scalar,-3/(rho*rho*p)*(D(M,0)+2*a*M),ham=True,mutation=mutation)
    return {'dimension':n,'Ricci_rows':len(ric),'Einstein_rows':len(ric),
            'arbitrary_scalar_coordinate_jets':True,'DC_y':divergence.serial(),
            'maximum_master_order_in_scalar':3,
            'Hamiltonian_witness':'deltaE_yy-2 M/rho^2 = 16 Psi (a^2-rho^-2-p^2/12+V/6) before Hamiltonian reduction'}


def separated_calculation(mutation):
    dy=lambda z:D(z,0,'separated',mutation)
    dt=lambda z:D(z,1,'separated',mutation)
    k=rho*(F(5,2)*a-V1/p)
    kval=-k if mutation=='k_sign' else k
    pm=alpha*theta*w();q=pm/rho**2
    coefB=3*alpha/(rho**3*p)
    actual_chi=coefB*theta*(w(0,1)/ell+kval*w())
    equal('actual_chi_constraint',actual_chi,-3*dy(pm)/(rho*rho*p))
    equal('actual_alpha_log_derivative',dy(alpha)/alpha,V1/p-F(5,2)*a)
    equal('actual_k_log_inward_derivative',-rho*dy(alpha)/alpha,k)
    QQ=2*(V1/p)**2-10*a*V1/p+F(35,4)*a*a-V2-p*p/6
    if mutation=='Q_a2':QQ+=a*a/4
    if mutation=='Q_p2':QQ+=p*p/3
    if mutation=='Q_V2':QQ+=2*V2
    # Radial transformation is checked separately, preserving +5/2.
    radial=dy(dy(pm))-(2*a+2*dy(p)/p)*dy(pm)-F(2,3)*p*p*pm
    radial_target=alpha*theta/rho**2*(w(0,2)/ell**2)+alpha*theta*(-QQ+F(5,2)/rho**2)*w()
    if mutation=='lose_radial_constant':radial_target-=F(5,2)*pm/rho**2
    equal('radial_alpha_Q_plus_5_over_2',radial,radial_target)
    eigen=-6 if mutation=='angular_eigenvalue' else -8
    box=lambda z:-dt(dt(z))-3*H*dt(z)+eigen*(1-H*H)*z
    master_const=0 if mutation=='missing_master_4' else 4
    temporal=box(theta*w())+master_const*theta*w()
    temporal_target=theta*(-w(2,0)/ell**2+(-F(5,2)+F(35,4)*H*H)*w())
    if mutation=='lose_temporal_constant':temporal_target+=F(5,2)*theta*w()
    equal('temporal_l2_minus_5_over_2',temporal,temporal_target)
    potential=rho*rho*QQ-F(35,4)*H*H
    ellpower=1 if mutation=='residual_ell' else 2
    residual=w(2,0)-w(0,2)+ell**ellpower*potential*w()
    M=radial+(box(pm)+master_const*pm)/rho**2
    equal('separated_master_conversion',M,-alpha*theta/(ell**2*rho**2)*residual)
    # Explicit inward-coordinate product rule is independent of differentiating
    # the factored expression.  All third master jets remain free.
    hratio=V1/p
    hratio_y=V2-hratio*hratio+4*a*hratio
    Q_y=(4*hratio-10*a)*hratio_y+(-10*hratio+F(35,2)*a)*(-rho**-2-p*p/3)-p*V3-p*(V1-4*a*p)/3
    potential_r=-ell*rho*(2*a*rho*rho*QQ+rho*rho*Q_y)
    radial_residual=w(2,1)-w(0,3)+ell*ell*(potential_r*w()+potential*w(0,1))
    equal('master_residual_r_chain_rule',radial_residual,-ell*rho*dy(residual))
    if mutation=='scalar_one_point':radial_residual=ZERO
    scalar_ell=1 if mutation=='scalar_derivative_ell' else ell
    factor=-3*alpha*theta/(ell**2*rho**5*p)*(radial_residual/scalar_ell+k*residual)
    scalar=-3/(rho*rho*p)*(dy(M)+2*a*M)
    equal('separated_scalar_full_radial_residual',scalar,factor)
    # Registered ordinary-y conversion and its essential rho derivative.
    yy=dy(dy(w()))
    yy_target=w(0,2)/(ell*ell*rho*rho)+a*w(0,1)/(ell*rho)
    if mutation=='drop_variable_rho':yy_target-=a*w(0,1)/(ell*rho)
    equal('ordinary_second_y_conversion',yy,yy_target)
    comparator_defect=residual.substitute('W:2,0',w(0,2))
    equal('free_comparator_retains_potential_defect',comparator_defect,ell*ell*potential*w())
    need(bool(comparator_defect.d),'comparator defect must be nonzero')
    if mutation=='comparator_claim':equal('false_comparator_natural_equation',comparator_defect,0)
    return {'eigenvalue':'-8','radial_constant':'+5/2 rho^-2',
            'temporal_constant':'-5/2','constant_sum':'0',
            'master_factor':'-alpha theta/(ell^2 rho^2)',
            'scalar_factor':'-3 alpha theta/(ell^2 rho^5 p) * (Rmaster_r/ell+k Rmaster)',
            'Q':QQ.serial(),'residual_radial_derivative':radial_residual.serial(),
            'nonzero_free_comparator_defect':comparator_defect.serial()}


def closure_auth(repo,pins_path):
    pinsraw=pins_path.read_bytes()
    need(hashlib.sha256(pinsraw).hexdigest()=='b67a008b2a6f8799190ed87e0239fdbf4e59b5ee93a054ae04eaf00482d62960','independent fixed pin file authentication')
    pins=json.loads(pinsraw);snap={}
    for rel,row in sorted(pins['direct'].items()):
        path=repo/rel
        need(not path.is_symlink(),'source symlink '+rel)
        data=path.read_bytes();sha=hashlib.sha256(data).hexdigest()
        need(len(data)==row['bytes'] and sha==row['sha256'],'source authentication '+rel)
        snap[rel]=sha
    counts={}
    for root,manifest_hash in sorted(pins['capsules'].items()):
        cap=repo/root;manifest=cap/'MANIFEST.json'
        need(not cap.is_symlink(),'capsule root symlink '+root)
        raw=manifest.read_bytes()
        need(hashlib.sha256(raw).hexdigest()==manifest_hash,'manifest authentication '+root)
        entries=json.loads(raw)['files'];listed={'MANIFEST.json'}
        for row in entries:
            rel=row['path'];path=cap/rel
            need(not Path(rel).is_absolute() and '..' not in Path(rel).parts,'unsafe capsule member')
            need(not path.is_symlink(),'capsule symlink '+root+'/'+rel)
            data=path.read_bytes()
            need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'capsule authentication '+root+'/'+rel)
            listed.add(rel)
        members=list(cap.rglob('*'))
        need(not any(x.is_symlink() for x in members),'capsule tree symlink '+root)
        actual={str(x.relative_to(cap)) for x in members if x.is_file()}
        need(actual==listed,'capsule closure '+root)
        counts[root]=len(entries)
    return snap,counts


def isolated_source_bindings(repo,mutation):
    path=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    tree=ast.parse(path.read_text())
    fun=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='frame_amplitudes')
    # Isolated AST is the authentic reconstruction; it has no imported globals.
    ns={};exec(compile(ast.Module(body=[fun],type_ignores=[]),str(path),'exec'),ns)
    class Jet:
        def __init__(self,v):self.v=P.cv(v)
        def __add__(self,b):return Jet(self.v+(b.v if isinstance(b,Jet) else b))
        __radd__=__add__
        def __neg__(self):return Jet(-self.v)
        def __sub__(self,b):return self+-b
        def __rsub__(self,b):return -self+b
        def __mul__(self,b):return Jet(self.v*(b.v if isinstance(b,Jet) else b))
        __rmul__=__mul__
        def __truediv__(self,b):return Jet(self.v/(b.v if isinstance(b,Jet) else b))
        def __rtruediv__(self,b):return Jet(P.cv(b)/self.v)
        def scaled_derivative(self,axis):
            def rule(name):
                if not name.startswith('W:'):return ZERO
                i,j=map(int,name[2:].split(','))
                return w(i+1,j) if axis=='x' else w(i,j+1)
            return Jet(self.v.derivative(rule))
        def dx(self):return self.scaled_derivative('x')
        def dr(self):return self.scaled_derivative('r')
    # Set gauge jets to zero for the authenticated longitudinal coefficient join.
    ns['ELL']=Jet(ell)
    rr=alpha/rho**2;bb=3*alpha/(rho**3*p);kk=rho*(F(5,2)*a-V1/p)
    amps,scalar=ns['frame_amplitudes'](Jet(w()),[Jet(0)]*3,*map(Jet,(rho,p,a,rr,bb,kk,P.var('c'),H,theta)))
    q=rr*theta*w()
    for key,multiple in {'A':-4,'B':0,'C':0,'D':-2,'E':0,'L':2,'M':0}.items():
        equal('source_longitudinal_'+key,amps[key].v,multiple*q)
    equal('source_longitudinal_chi',scalar.v,bb*theta*(w(0,1)/ell+kk*w()))
    physical=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='physical_derivatives')
    physns={'ELL':ell}
    exec(compile(ast.Module(body=[physical],type_ignores=[]),str(path),'exec'),physns)
    class OrdinaryJet:
        def __init__(self,v):self.v=v
        def val(self):return self.v
        def derivative(self,i,j):return P.var('J:%d,%d'%(i,j))
    actual_adapter=physns['physical_derivatives'](OrdinaryJet(P.var('J:0,0')),OrdinaryJet(rho),OrdinaryJet(a))
    J=lambda i,j:P.var('J:%d,%d'%(i,j))
    expected_adapter={'':J(0,0),'y':-J(0,1)/(ell*rho),'t':J(1,0)/ell,
        'yy':J(0,2)/(ell*ell*rho*rho)+a*J(0,1)/(ell*rho),
        'yt':-J(1,1)/(ell*ell*rho),'tt':J(2,0)/(ell*ell)}
    for suffix,target in expected_adapter.items():equal('source_physical_adapter_'+(suffix or 'value'),actual_adapter[suffix],target)
    # Complete plus-sign pullback was independently bound in PR106; pin it here
    # and inspect the actual source assignment, without importing its result.
    assign=next(x for x in ast.walk(fun) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='scalar' for t in x.targets))
    expected='v - p * f if injection == \'wrong_scalar_shift\' else v + p * f'
    equal_ast=ast.dump(assign.value,include_attributes=False)==ast.dump(ast.parse(expected,mode='eval').body,include_attributes=False)
    need(equal_ast,'actual plus scalar Lie assignment')
    pullback_fun=ast.parse(ast.unparse(fun)).body[0]
    if mutation=='lie_scalar_sign':
        live_assign=next(x for x in ast.walk(pullback_fun) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='scalar' for t in x.targets))
        live_assign.value.orelse=ast.parse('v-p*f',mode='eval').body
    pullns={'ELL':Jet(ell)}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[pullback_fun],type_ignores=[])),str(path),'exec'),pullns)
    _,shifted=pullns['frame_amplitudes'](Jet(w()),[Jet(P.var('f')),Jet(0),Jet(0)],*map(Jet,(rho,p,a,rr,bb,kk,P.var('c'),H,theta)))
    equal('pure_gauge_scalar_plus_sign',shifted.v-scalar.v,p*P.var('f'))
    CHECKS.append('authenticated_plus_scalar_Lie_assignment')
    # Bind the actual polynomial to its generating superpotential, rather than
    # assigning the name V to a similarly named historical potential.
    corepath=path.with_name('interval_core.py');coretree=ast.parse(corepath.read_text())
    Uassign=next(x for x in coretree.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='U' for t in x.targets))
    coeffs=eval(compile(ast.Expression(Uassign.value),str(corepath),'eval'),{'F':F,'list':list,'map':map})
    eta=P.var('eta');superpotential=F(1,3)+eta**2+eta**3/3
    de=lambda z:z.derivative(lambda name:P.cv(1) if name=='eta' else ZERO)
    model=de(superpotential)**2/2-F(2,3)*superpotential**2
    actual=sum((co*eta**i for i,co in enumerate(coeffs)),ZERO)
    equal('actual_scalar_potential_polynomial',actual,model)
    bg=next(x for x in coretree.body if isinstance(x,ast.FunctionDef) and x.name=='background')
    rhsassign=next(x for x in ast.walk(bg) if isinstance(x,ast.Assign) and isinstance(x.value,ast.Dict) and any(isinstance(t,ast.Name) and t.id=='rhs' for t in x.targets))
    namespace={'rho':rho,'e':eta,'p':p,'a':a,'alpha':alpha,'h':V1/p,
               'F':F,'poly':lambda cs,e:V1,'UP':[None,None]}
    rhs=eval(compile(ast.Expression(rhsassign.value),str(corepath),'eval'),namespace)
    expected={'rho':-rho*rho*a,'eta':-rho*p,'p':-rho*(V1-4*a*p),
              'a':1/rho+rho*p*p/3,'alpha':rho*(F(5,2)*a-V1/p)*alpha}
    for key in expected:equal('source_background_inward_'+key,rhs[key],expected[key])
    # The final k,C,B expressions are individually isolated authenticated ASTs.
    for name,target in [('k',kk),('C',rr),('B',bb)]:
        assignment=next(x for x in bg.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in x.targets))
        value=eval(compile(ast.Expression(assignment.value),str(corepath),'eval'),namespace)
        equal('source_background_'+name,value,target)
    potpath=repo/'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/check_actual_potential.py'
    pottree=ast.parse(potpath.read_text())
    formal=next(x for x in pottree.body if isinstance(x,ast.FunctionDef) and x.name=='formal_identities')
    Qassign=next(x for x in formal.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Q' for t in x.targets))
    Qsource=eval(compile(ast.Expression(Qassign.value),str(potpath),'eval'),{'F':F,'h':V1/p,'a':a,'p':p,'u2':V2})
    Qexact=2*(V1/p)**2-10*a*V1/p+F(35,4)*a*a-V2-p*p/6
    equal('source_potential_Q',Qsource,Qexact)
    runfun=next(x for x in pottree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
    numericQ=next(x for x in runfun.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Q' for t in x.targets))
    numericpotential=next(x for x in runfun.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='potential' for t in x.targets))
    equal('source_numeric_potential_Q',eval(compile(ast.Expression(numericQ.value),str(potpath),'eval'),{'F':F,'h':V1/p,'a':a,'p':p,'u2':V2}),Qexact)
    Tsymbol=P.var('T')
    for endpoint,choose,target in [('lower',lambda lo,hi:hi,rho*rho*Qexact-F(35,4)*Tsymbol*Tsymbol),('upper',lambda lo,hi:lo,rho*rho*Qexact)]:
        value=eval(compile(ast.Expression(numericpotential.value),str(potpath),'eval'),{'F':F,'r':rho,'Q':Qexact,'T':Tsymbol,'I':choose})
        equal('source_numeric_potential_'+endpoint,value,target)
    astpins={name:hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest()
             for name,node in [('frame_amplitudes',fun),('physical_derivatives',physical),('background',bg),('U_assignment',Uassign),('formal_Q_assignment',Qassign),('numeric_Q_assignment',numericQ),('numeric_potential_assignment',numericpotential)]}
    return {'functions_evaluated':['actual_first_jets.frame_amplitudes and physical_derivatives (isolated AST only)'],
            'isolated_expression_ASTs':['interval_core.U','interval_core.background.rhs, k, C, B','check_actual_potential.formal_identities.Q'],
            'canonical_AST_sha256':astpins,
            'actual_PDE_definition':'authenticated HIGHER_JETS.md and ACTUAL_POTENTIAL_WZZ_DERIVATION.md; exact smooth solution is retained analytic premise',
            'numeric_potential_role':'rational enclosure rho^2 Q-[0,35/4 T^2], not an exact point evaluator or PDE solver',
            'scientific_modules_imported':False,'pure_gauge_cancellation':'standard exact natural covariance on E0=S0=0',
            'pullback_binding':'PR106 independently bound full plus Lie pair; authenticated complete capsule'}


CONTROL_EXPECTATIONS={
    'chi_sign':'natural_E_00','alpha_sign':'actual_chi_constraint',
    'k_sign':'actual_chi_constraint','drop_Vh':'natural_E_00',
    'drop_Vprime':'natural_E_00','drop_DC':'direct_scalar_variation',
    'no_hamiltonian':'natural_E_00','angular_eigenvalue':'temporal_l2_minus_5_over_2',
    'missing_master_4':'temporal_l2_minus_5_over_2','theta_frozen':'temporal_l2_minus_5_over_2',
    'lose_radial_constant':'radial_alpha_Q_plus_5_over_2','lose_temporal_constant':'temporal_l2_minus_5_over_2',
    'Q_a2':'radial_alpha_Q_plus_5_over_2','Q_p2':'radial_alpha_Q_plus_5_over_2','Q_V2':'radial_alpha_Q_plus_5_over_2',
    'radial_sign':'actual_chi_constraint','residual_ell':'separated_master_conversion',
    'scalar_one_point':'separated_scalar_full_radial_residual','drop_variable_rho':'ordinary_second_y_conversion',
    'lie_scalar_sign':'pure_gauge_scalar_plus_sign',
    'scalar_derivative_ell':'separated_scalar_full_radial_residual',
    'comparator_claim':'false_comparator_natural_equation'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--run-controls',action='store_true')
    ap.add_argument('--mutation',choices=['none']+list(CONTROL_EXPECTATIONS),default='none')
    args=ap.parse_args();need(not args.output.exists(),'output must be new')
    pins=Path(__file__).with_name('SOURCE_PINS.json')
    before,counts=closure_auth(args.repo_root,pins)
    tensor=tensor_calculation(args.mutation)
    separated=separated_calculation(args.mutation)
    sources=isolated_source_bindings(args.repo_root,args.mutation)
    after,again=closure_auth(args.repo_root,pins)
    need(before==after and counts==again,'source readback authentication')
    report={'status':'PASS_INDEPENDENT_EXACT_FIRST_NATURAL_REVIEW',
            'checks':CHECKS,'tensor':tensor,'separated':separated,'source_binding':sources,
            'source_sha256':before,'complete_capsule_members':counts,
            'normal_optimized_decisions':'no assert; all Fraction arithmetic',
            'scope':'Exact implementation-to-retained-proof algebra; actual smooth neighbourhood PDE and background are analytic premises; no new existence theorem'}
    if args.run_controls:
        controls={}
        with tempfile.TemporaryDirectory(prefix='sol-natural-') as tmp:
            for name,expected in CONTROL_EXPECTATIONS.items():
                out=Path(tmp)/(name+'.json')
                command=[sys.executable,'-I','-S','-B']
                if sys.flags.optimize:command.append('-O')
                command += [str(Path(__file__).resolve()),'--repo-root',str(args.repo_root),'--output',str(out),'--mutation',name]
                done=subprocess.run(command,capture_output=True,text=True)
                need(done.returncode!=0 and not out.exists(),'control failed to reject '+name)
                lines=done.stderr.strip().splitlines()
                need(bool(lines) and lines[-1].startswith('FAIL: '+expected+' | '),'control rejected for wrong reason '+name+' '+done.stderr)
                controls[name]={'rejected':True,'identity':expected,'exact_remainder':json.loads(lines[-1].split(' | ',1)[1])}
        report['corrupted_subprocess_controls']=controls
    encoded=(json.dumps(report,indent=2,sort_keys=True)+'\n').encode()
    args.output.write_bytes(encoded)
    print(report['status']+'; checks='+str(len(CHECKS)))


if __name__=='__main__':
    try: main()
    except Failure as exc:
        print('FAIL: '+str(exc),file=sys.stderr)
        sys.exit(2)
