#!/usr/bin/env python3
"""Authenticated bounded actual first-field oracle on the retained full domain.

The reflected domain is covered by the original *global* differentiated error
certificate.  This is not an extension of the sharper exterior certificate.
All imports of project code follow closed-tree authentication.  No output from
this module supplies a numerical TT response or a new physical/PDE premise.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import types
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
DEFAULT_REPO = Path('/workspace/hdblast-absolute-publish-31c')
CAPSULES = {
    'G': ('HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
          '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),
    'S': ('HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026',
          '204e647b9f68f2af169b052deabee7c5d7f4bd71195fae31cefbfc2d10265777'),
}
APEX = (F(5,4), F(1,8))
INJECTIONS = ('none', 'wrong_reflection', 'wrong_branch', 'wrong_shell_sign',
              'omit_G', 'omit_dynamic_lift')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(value):
    if hasattr(value, 'lo') and hasattr(value, 'hi'):
        return [str(value.lo), str(value.hi)]
    if isinstance(value, dict):
        return {str(key): encode(item) for key,item in value.items()}
    if isinstance(value, (list,tuple)):
        return [encode(item) for item in value]
    if isinstance(value,F):
        return str(value)
    return value


def authenticate(repo):
    """Authenticate all exact members, reject symlinks/unlisted files."""
    result = {}
    for key,(label,pin) in CAPSULES.items():
        root = repo/'research'/label
        need(root.is_dir() and not root.is_symlink(), 'ROOT_TYPE:'+key)
        manifest = root/'MANIFEST.json'
        need(manifest.is_file() and not manifest.is_symlink() and sha(manifest)==pin,
             'MANIFEST:'+key)
        seen = {'MANIFEST.json'}
        for row in json.loads(manifest.read_text())['files']:
            rel = Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
                 'MEMBER_PATH:'+key)
            seen.add(str(rel)); path = root/rel
            need(path.is_file() and not path.is_symlink(), 'MEMBER_TYPE:'+key+'/'+str(rel))
            need(path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],
                 'MEMBER_BYTES:'+key+'/'+str(rel))
        for path in root.rglob('*'):
            need(not path.is_symlink(), 'TREE_SYMLINK:'+key)
            if path.is_file():
                need(str(path.relative_to(root)) in seen, 'UNLISTED_FILE:'+key)
        result[key] = {'label':label, 'manifest_sha256':pin,
                       'members_including_manifest':len(seen)}
    return result


def load(name,path,replace=None):
    if replace:
        raw = path.read_text()
        for old,new in replace:
            need(raw.count(old)==1, 'EXACT_MUTATION_ANCHOR')
            raw=raw.replace(old,new)
        module=types.ModuleType(name); module.__file__=str(path)
        exec(compile(raw,str(path),'exec'), module.__dict__)
        return module
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class Context:
    """Bounded cell API; call finish() to reauthenticate after use.

    Actual input domain: 0<=x<=3, r>=0, x+r<=8.  Each derivative enclosure
    is independently certified for the same retained smooth solution.
    """
    def __init__(self,repo=DEFAULT_REPO,injection='none'):
        need(injection in INJECTIONS,'INJECTION_NAME')
        self.repo=Path(repo).resolve(); self.own_pin=sha(__file__)
        self.before=authenticate(self.repo)
        self.G=self.repo/'research'/CAPSULES['G'][0]
        self.S=self.repo/'research'/CAPSULES['S'][0]
        self.numeric=self.G/'numeric'
        for path in (self.numeric,self.G/'projection',self.S/'numeric'):
            sys.path.insert(0,str(path))
        changes=[]
        if injection=='wrong_reflection':
            changes=[('(-1)**(n+1)*z[n]', '(-1)**n*z[n]')]
        if injection=='wrong_branch':
            changes=[('lin=x-2*r if negative else -x', 'lin=-x if negative else -x')]
        self.j=load('actual_first_jets',self.numeric/'actual_first_jets.py',changes)
        sys.modules['actual_first_jets']=self.j
        self.shell=load('reflected_shell_values',self.numeric/'shell_source_values.py')
        self.coeff=load('reflected_source_coefficients',self.G/'projection/source_coefficients.py')
        self.dyadic=load('reflected_fast_dyadic',self.numeric/'fast_dyadic.py')
        self.selected=load('reflected_selected_magnetic',self.S/'numeric/evaluate_selected_magnetic.py')
        self.data,self.errors,self.pins=self.j.inputs()
        inherited=load('reflected_inherited_error_replay',self.numeric/'check_inherited_errors.py')
        replay=inherited.run('drop_dynamic_lift' if injection=='omit_dynamic_lift' else 'none')
        need(replay==self.errors,'EXACT_GLOBAL_C3_ERROR_REPLAY')
        self.map_path=self.G/'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
        self.maps=self.coeff.load_maps(self.map_path)
        rows=self.maps['maps']['metric_forcing']['ab']
        need(rows[5]==rows[6] and len(rows[5])==78,'ACTUAL_EQUAL_78_TERM_TT_ROWS')
        self.selected_map=json.loads((self.S/'numeric/CONSENSUS_SELECTED_MAGNETIC_MAP.json').read_text())
        self.transcription=self.selected.verify_transcription(self.selected_map)
        self.injection=injection
        self.controls=self.check_controls()

    def check_controls(self):
        """Independent expanded polynomial witnesses exercise production code."""
        j=self.j; original=j.qodd
        def manufactured(s,order,data):
            need(not (s.lo<0<s.hi),'MANUFACTURED_BRANCH_STRADDLE')
            sg=1 if s.lo>=0 else -1
            return [sg*math.factorial(4)//math.factorial(4-n)*s**(4-n)
                    for n in range(order+1)]
        # Exercise the actual parity extension separately on a negative cell.
        positive=j.qpositive
        try:
            j.qpositive=lambda r,order,data: [
                F(math.factorial(4),math.factorial(4-n))*r**(4-n)
                for n in range(order+1)]
            got=j.qodd(j.I(-F(9,8)),3,self.data)
            expected=[-F(math.factorial(4),math.factorial(4-n))*(-F(9,8))**(4-n)
                      for n in range(4)]
            need(all(v.lo<=e<=v.hi for v,e in zip(got,expected)),
                 'REFLECTED_ODD_DERIVATIVE_PARITY')
        finally:
            j.qpositive=positive
        try:
            j.qodd=manufactured
            w,p=j.comparisons(j.I(APEX[0]),j.I(APEX[1]),self.data)
            wpoly={(3,1):F(4),(1,3):F(4)}
            ppoly={(4,1):F(5,2),(3,2):F(-2),(2,3):F(5),
                   (1,4):F(-2),(0,5):F(1,2)}
            def derivative(poly,i,k):
                return sum((co*F(math.factorial(a),math.factorial(a-i))*
                            F(math.factorial(b),math.factorial(b-k))*
                            APEX[0]**(a-i)*APEX[1]**(b-k)
                            for (a,b),co in poly.items() if a>=i and b>=k),F(0))
            c=6*j.ELL/130
            for i,k in j.KEYS:
                ew=derivative(wpoly,i,k); ep=c*derivative(ppoly,i,k)
                aw,ap=w.derivative(i,k),p.derivative(i,k)
                need(aw.lo<=ew<=aw.hi and ap.lo<=ep<=ap.hi,
                     'REFLECTED_COMPARATOR_POLYNOMIAL_JET')
            lhs=p.derivative(2,0)-p.derivative(0,2)
            rhs=c*derivative(wpoly,0,1)
            need(lhs.lo<=rhs<=lhs.hi,'REFLECTED_COMPARATOR_WAVE_EQUATION')
        finally:
            j.qodd=original
        self.shell.algebra_check()
        # True adapter production path, not just the closed formula.
        raw=self.shell.raw_fields(*map(F,(2,3,5,7,10,3,7,11,23)))
        qt=self.shell.projected(raw)['ab'][5]
        g=self.shell_convert(qt,F(1,10),F(13))
        need(g==F(416,3),'PHYSICAL_TT_SHELL_SOURCE_SIGN_AND_NORMALIZATION')
        return {'negative_argument_odd_derivative_parity_orders':4,
                'reflected_expanded_polynomial_actual_comparator_jets':20,
                'reflected_comparator_wave_equation':True,
                'physical_shell_independent_unit_normal_algebra':True,
                'physical_TT_shell_manufactured_exact_G':'416/3'}

    def shell_convert(self,q,ell,rho):
        if self.injection=='omit_G':
            return 0*q
        sign=-1 if self.injection=='wrong_shell_sign' else 1
        return sign*2*ell*rho*q

    def error_caps(self,x,r):
        if r.lo>=x.hi:
            ew,eu=self.j.exterior_caps(); kind='exterior_backward_closed'
        else:
            ew=[F(v['exact_upper']) for v in self.errors['master_error_budget_by_scaled_total_order']['actual_minus_comparison']]
            eu={k:[F(v['exact_upper']) for v in row] for k,row in self.errors['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'].items()}
            kind='global_with_dynamic_lift_and_reflected_shell_terms'
        return ew,eu,kind

    def raw_jets(self,x,r):
        j=self.j;x=j.I.cv(x);r=j.I.cv(r)
        need(0<=x.lo<=x.hi<=3 and 0<=r.lo<=r.hi and x.hi+r.hi<=8,'RAW_JET_DOMAIN')
        w,p=j.comparisons(x,r,self.data);u=[p,j.T(),j.T()]
        ew,eu,kind=self.error_caps(x,r)
        for i,k in j.KEYS:
            den=math.factorial(i)*math.factorial(k)
            w.c[i,k]+=j.I(-ew[i+k],ew[i+k])/den
            for ch in range(3):
                cap=eu[f'{i},{k}'][ch];u[ch].c[i,k]+=j.I(-cap,cap)/den
        if r.lo==r.hi==0:
            for i in range(4):u[0].c[i,0]=j.I(0)
            for i in range(3):u[1].c[i,1]=j.I(0);u[2].c[i,1]=j.I(0)
        return {'master':w,'f':u[0],'b':u[1],'d':u[2]},(ew,eu,kind)

    def first(self,x,r):
        """All actual ordinary physical metric/scalar two-jets on whole cell."""
        return self.j.actual_jet(x,r)

    def bulk_tt(self,x,r):
        """Actual registered F_T cells; no nominal-as-actual substitution."""
        m,s,b=self.first(x,r); values={**m,**s,**b}
        values={key:self.dyadic.D(value) if hasattr(value,'lo') else value
                for key,value in values.items()}
        return self.coeff.evaluate_encoded(self.maps['maps']['metric_forcing']['ab'][5],values)

    def shell_tt(self,x):
        """Physical-side actual nonhomogeneous TT Neumann source at r=0."""
        j=self.j;m,s,b=self.first(x,j.I(0))
        data,_,_=j.inputs()
        eta=j.I(*data['shell']['boxes']['eta_b']['exact'])
        raw=self.shell.raw_fields(m['A'],m['D'],m['E'],m['L'],m['M'],s['s'],
                                  b['a'],b['p'],4*(1+eta))
        q=self.shell.projected(raw)['ab'][5]
        return self.shell_convert(q,j.ELL,b['rho'])

    def point(self,x=APEX[0],r=APEX[1]):
        j=self.j;m,s,b=self.first(j.I(x),j.I(r))
        raw,caps=self.raw_jets(j.I(x),j.I(r))
        # Rebuild via raw independently added derivative errors and check the
        # exact same frame reconstruction endpoint identities against first().
        rr,eta,pp,aa,hh,RR,BB,kk,_,_=j.core.background(j.I(r),self.data,'none')
        rho,p,a,R,B,k=map(j.radial,(rr,pp,aa,RR,BB,kk))
        c,H,theta=j.temporal(j.I(x))
        amps,scalar=j.frame_amplitudes(raw['master'],[raw[ch] for ch in ('f','b','d')],
                                     rho,p,a,R,B,k,c,H,theta)
        rebuilt={key+suffix:value for key,amp in amps.items()
                 for suffix,value in j.physical_derivatives(amp,rho,a).items()}
        rebuilt_s={'s'+suffix:value for suffix,value in j.physical_derivatives(scalar,rho,a).items()}
        need(encode(m)==encode(rebuilt) and encode(s)==encode(rebuilt_s),
             'RAW_ACTUAL_DERIVATIVE_RECONSTRUCTION_ENDPOINT_IDENTITY')
        symbols={**m,**s,**b,'u_x':j.I(0)}
        grouped=self.selected.supplied_groups(symbols)
        expanded={name:self.coeff.evaluate_encoded(terms,symbols)
                  for name,terms in self.selected_map['components'].items()}
        qg,qe=grouped['quadratic'],expanded['quadratic']
        need(max(qg.lo,qe.lo)<=min(qg.hi,qe.hi),'QUADRATIC_EQUIVALENT_ROUTES_OVERLAP')
        quadratic=j.I(max(qg.lo,qe.lo),min(qg.hi,qe.hi))
        return {'apex':{'x':str(x),'r':str(r)},'symbols':{**m,**s,**b},
                'raw_actual_scaled_derivatives':{name:{f'{i},{k}':jet.derivative(i,k)
                   for i,k in j.KEYS} for name,jet in raw.items()},
                'actual_error_caps':{'kind':caps[2],'master_by_total_order':caps[0],
                                    'gauge_by_mixed_order':caps[1]},
                'raw_to_physical_twojet_endpoint_identity':True,
                'quadratic':{'grouped_F':grouped['F'],'grouped_4F_over_R':qg,
                             'expanded_14_term_4F_over_R':qe,'intersection':quadratic,
                             'excludes_zero':quadratic.lo>0 or quadratic.hi<0},
                'linear_TT_response_supplied':False}

    def finish(self):
        after=authenticate(self.repo)
        need(after==self.before and sha(__file__)==self.own_pin,'SOURCE_CHANGED')
        return {'before':self.before,'after':after,
                'project_imports_follow_closed_authentication':True,
                'source_unchanged':True,'wrapper_sha256':self.own_pin}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,default=DEFAULT_REPO)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject-error',choices=INJECTIONS,default='none')
    args=parser.parse_args();out=args.output
    need(not out.exists() and not out.is_symlink(),'FRESH_OUTPUT_REQUIRED')
    need(out.parent.resolve()==HERE,'OUTPUT_MUST_USE_EXCLUSIVE_SCRATCH_ROOT')
    ctx=Context(args.repo_root,args.inject_error)
    point=ctx.point();j=ctx.j
    shell_rows=[]
    for i in range(8):
        x=j.I(F(9*i,64),F(9*(i+1),64))
        shell_rows.append({'x':x,'G_T':ctx.shell_tt(x)})
    endpoint=ctx.shell_tt(j.I(F(9,8)))
    result={'schema':'HDBLAST_REFLECTED_ACTUAL_FIRST_TWOJET_V1',
            'status':'CERTIFIED_CONDITIONAL_GLOBAL_ACTUAL_FIRST_JETS_AT_POSTCONTACT_POINT',
            'point':point,'shell_history':{'x':['0','9/8'],'whole_gap_free_cells':8,
                                         'rows':shell_rows,'endpoint_G_T':endpoint},
            'controls':ctx.controls,'selected_map_transcription':ctx.transcription,
            'source_authentication':ctx.finish(),
            'domain':{'retained_global':'0<=x<=3,r>=0,x+r<=8',
              'Pdagger_backward_region':'0<=s<=5/4,max(0,s-9/8)<=q<=11/8-s',
              'shell_history':'0<=s<=9/8; all shell jets are physical-side one-sided'},
            'limits':['Global actual C3 derivative errors remain retained; no exterior-domain relabeling.',
              'Retained smooth background/master/gauge and differentiated integral identities remain conditional premises.',
              'No numerical first refinement, TT response, finite-amplitude family, or new physics claim.']}
    with out.open('x') as handle:
        json.dump(encode(result),handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':result['status'],'output':str(out),'sha256':sha(out),
                      'quadratic_display_only':[float(x) for x in
                        (point['quadratic']['intersection'].lo,point['quadratic']['intersection'].hi)]}))


if __name__=='__main__':
    main()
