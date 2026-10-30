#!/usr/bin/env python3
"""Exact source-to-proof binding of actual physical Cauchy initial data.

This is a formal Laurent identity audit, not a test at finitely many points.
Project source is authenticated before its AST is executed and after use.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode=True
ENGINE_PIN='02d2c85ec4188149f4a3c53c62a7781d17cf0a2f8e4b53cf3da5b697689e0b5b'
INPUTS={
 'registered':('research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py','2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1'),
 'background':('research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/numeric/evaluate_initial_profiles.py','d9cb8529e47b9590ec8dc94f2f4484b06dbeb0ddee6c16015feccde3e4440240'),
 'physical_proof':('research/HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md','2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed'),
 'actual_bounds':('research/HDBLAST_COMMON_GAUGE_FORCING_20261016/bulk-producer/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json','fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0'),
}

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(data):return hashlib.sha256(data).hexdigest()
def authenticate(root):
    out={}
    for key,(rel,pin) in INPUTS.items():
        path=root/rel
        need(path.is_file() and not path.is_symlink(),'SOURCE_TYPE:'+key)
        raw=path.read_bytes();need(sha(raw)==pin,'SOURCE_PIN:'+key)
        out[key]={'path':rel,'bytes':len(raw),'sha256':pin}
    return out
def load_engine():
    path=Path(__file__).resolve().with_name('laurent.py')
    need(sha(path.read_bytes())==ENGINE_PIN,'ENGINE_PIN')
    spec=importlib.util.spec_from_file_location('physical_laurent',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def registered_prefix(source,alg,mutation=None):
    if mutation=='wrong_hamiltonian_denominator':source=source.replace('(-6, 10, 42)','(-6, 12, 42)')
    if mutation=='omit_Pi_cross_term':source=source.replace('d*u*vy+','0*d*u*vy+')
    if mutation=='wrong_scalar_velocity_sign':source=source.replace('E*(v-d*u)','E*(v+d*u)')
    if mutation=='wrong_hamiltonian_gradient':source=source.replace('36*u*u-2*v*v','32*u*u-2*v*v')
    if mutation=='wrong_y_z_derivative_scale':source=source.replace('-ir*value.dz()','-self.ell*ir*value.dz()')
    tree=ast.parse(source)
    cls=next(x for x in tree.body if isinstance(x,ast.ClassDef) and x.name=='Context')
    fn=next(x for x in cls.body if isinstance(x,ast.FunctionDef) and x.name=='_registered')
    body=[]
    for statement in fn.body:
        if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='zero' for t in statement.targets):break
        body.append(statement)
    fn=ast.FunctionDef(name='registered_prefix',args=fn.args,body=body+[ast.Return(value=ast.Call(func=ast.Name(id='locals',ctx=ast.Load()),args=[],keywords=[]))],decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    env={'F':F};exec(compile(module,'authenticated_registered_prefix','exec'),env)
    fake=SimpleNamespace(init=SimpleNamespace(poly=alg.poly,UP=alg.UP),ell=F(1,10**6))
    s=alg.symbols();bg={key:s[key] for key in ('rho','eta','p','a','alpha')}
    return env['registered_prefix'](fake,bg,s['E'],s['Q0'])

def registered_background(source,alg):
    tree=ast.parse(source)
    fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='background')
    assignments=[node for node in ast.walk(fn) if isinstance(node,ast.Assign) and isinstance(node.value,ast.Dict) and any(isinstance(t,ast.Name) and t.id=='rhs' for t in node.targets)]
    need(len(assignments)==1,'BACKGROUND_RHS_UNIQUE')
    env=dict(alg.symbols());env.update({'poly':alg.poly,'UP':alg.UP,'e':alg.symbols()['eta'],'F':F})
    env['h']=alg.poly(alg.UP[1],env['e'])/env['p']
    return eval(compile(ast.Expression(assignments[0].value),'authenticated_background_rhs','eval'),env)

def registered_potential(source):
    tree=ast.parse(source);body=[]
    for statement in tree.body:
        if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('U','UP') for t in statement.targets):body.append(statement)
        if isinstance(statement,ast.For) and any(isinstance(n,ast.Name) and n.id=='UP' for n in ast.walk(statement)):body.append(statement)
    need(len(body)==3,'POTENTIAL_ASSIGNMENTS')
    env={'F':F};exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),'authenticated_potential','exec'),env)
    return tuple(tuple(row) for row in env['UP'])

def audit(root,mutation=None):
    pins=authenticate(root);alg=load_engine();s=alg.symbols()
    potential=registered_potential((root/INPUTS['background'][0]).read_text())
    need(potential==tuple(alg.UP[:5]),'ACTUAL_POTENTIAL_POLYNOMIAL')
    potential_json=json.loads((root/INPUTS['actual_bounds'][0]).read_text())
    need(tuple(map(F,potential_json['actual_model']['U_coefficients']))==alg.POTENTIAL,'ACTUAL_POTENTIAL_REGISTERED_JSON')
    impl=registered_prefix((root/INPUTS['registered'][0]).read_text(),alg,mutation)
    phy=alg.physical_table();rho,p,a=[s[key] for key in ('rho','p','a')]
    dy=alg.dy;rows=[]
    def identity(name,left,right=0):
        left,right=alg.Expr.cv(left),alg.Expr.cv(right);diff=left-right
        need(not diff,'NONZERO_IDENTITY:'+name+':terms='+str(len(diff.terms)))
        rows.append({'name':name,'left_collected_terms':len(left.terms),'right_collected_terms':len(right.terms),'residual_collected_terms':0})
    bg=registered_background((root/INPUTS['background'][0]).read_text(),alg)
    for key in ('rho','eta','p','a','alpha'):identity('actual_background_Dz_'+key,bg[key],alg.derivative_table()[key])
    identity('actual_scalar_background_equation',dy(p)+4*a*p,phy['vp'])
    identity('actual_radial_background_equation',dy(a),-1/rho**2-p*p/3)
    identity('actual_profile_u',impl['u'],phy['u'])
    identity('actual_profile_v',impl['v'],phy['v'])
    identity('actual_profile_ut',impl['ut'],phy['ut'])
    identity('actual_profile_vt',impl['vt'],phy['vt'])
    u,v,ut,vt,uy,vy,uyy,uyt=[phy[key] for key in ('u','v','ut','vt','uy','vy','uyy','uyt')]
    identity('first_reconstruction',p*v,-3*(uy+2*a*u))
    identity('time_first_reconstruction',p*vt,-3*(uyt+2*a*ut))
    identity('exponential_time_profile',vt,s['E']*(v-3*u/(rho*p)))
    # Distinct +20 and -32 angular curvature contributions are retained before
    # collecting; cancellation uses the exact background ODE, never a point box.
    firstH=(-6*uyy-36*a*uy-(24*dy(a)+48*a*a+4*p*p)*u
            +20*u/rho**2-32*u/rho**2-2*p*vy-2*phy['vp']*v)
    identity('first_physical_Hamiltonian',firstH)
    identity('first_physical_radial_momentum',(-3*(uyt+2*a*ut)-p*vt)/rho)
    # On the spherical block, P1=K1-gamma0 tr(K1) vanishes:
    # K1_yy=-2ut/rho, three equal angular entries=ut/rho.  The
    # first scalar momentum product has Pi0=0 and no angular phi0 gradient.
    first_trace=(-2*ut/rho)+3*(ut/rho)
    identity('first_physical_angular_momentum',ut/rho-first_trace)
    identity('Hamiltonian_A_unsimplified_vs_registered',phy['A'],impl['aa'])
    identity('Hamiltonian_D_vs_registered',phy['D'],impl['dd'])
    identity('physical_angular_momentum_F_equals_3K',phy['Fmom'],3*impl['K'])
    identity('physical_radial_momentum_G_minus_Fy_equals_pPi',phy['Gmom']-dy(phy['Fmom']),p*impl['Pi'])
    for j,(lam,sector) in enumerate(zip((0,8,24),phy['sectors'])):
        f=impl['fs'][j];zeta=-F(3,2)*dy(f)/p
        if mutation=='scalar_compensator_actual_sign':zeta=-zeta
        deltaR=-3*dy(dy(f))-12*a*dy(f)+2*(lam-3)*f/rho**2
        deltaScalar=-2*p*dy(zeta)-2*phy['vp']*zeta
        identity('Hamiltonian_scalar_sector_'+str((0,2,4)[j]),deltaR+deltaScalar,-sector)
    radialG=phy['Gmom']
    if mutation=='spacetime_affine_radial_momentum':radialG+=6*uy*ut/rho
    # A physical K2=q gamma0 in four spatial dimensions yields -3 dq.
    tracefactor=4 if mutation=='wrong_spatial_trace_factor' else 3
    radial=-tracefactor*dy(impl['K'])-p*impl['Pi']+radialG
    for ell,c in zip((0,2,4),(F(1,6),F(1,2),F(1))):
        identity('second_physical_radial_momentum_'+str(ell),c*radial)
    for ell,c in zip((2,4),(F(1,2),F(1))):
        identity('second_physical_angular_momentum_'+str(ell),c*(-tracefactor*impl['K']+phy['Fmom']))
    # Finite differentiation keeps all expressions in the packet-jet ideal.
    ideal_entries={name:phy[name] for name in ('u','v','ut','vt','A','D','Fmom','Gmom','K','Pi')}
    for j in range(3):
        ideal_entries['f_'+str(j)]=phy['fs'][j]
        ideal_entries['fy_'+str(j)]=phy['fy'][j]
        ideal_entries['zeta_'+str(j)]=-F(3,2)*phy['fy'][j]/p
    for name,expr in ideal_entries.items():
        minimum_degree=1 if name in ('u','v','ut','vt') else 2
        need(all(all(power[alg.INDEX['Q'+str(j)]]>=0 for j in range(9)) and
                 sum(power[alg.INDEX['Q'+str(j)]] for j in range(9))>=minimum_degree
                 for power in expr.terms),'PACKET_IDEAL:'+name)
    need(authenticate(root)==pins,'SOURCES_CHANGED_AFTER_AUDIT')
    need(sha(Path(__file__).resolve().with_name('laurent.py').read_bytes())==ENGINE_PIN,'ENGINE_CHANGED_AFTER_AUDIT')
    return {'status':'PASS_EXACT_ACTUAL_PHYSICAL_INITIAL_BINDING','input_sha256':pins,
      'arithmetic':{'coefficient_field':'Q','formal_generators':list(alg.NAMES),
        'division':'nonzero rho and p only, apart from rational constants',
        'zero_test':'every collected Laurent coefficient exactly zero',
        'scope':'identities valid at every point satisfying retained exact background/packet derivation'},
      'identities':rows,'identity_count':len(rows),
      'physical_constraint_modal_rows':{'Hamiltonian':3,'radial_momentum':3,'angular_gradient_momentum':2},
      'packet_ideal_checked':True,
      'packet_ideal_entries':sorted(ideal_entries),
      'actual_potential_polynomial_and_four_derivatives_bound':True,
      'scope':{'registered_actual_prefix_bound_to_physical_affine_sources':True,
        'all_point_in_actual_regular_branch':True,'initial_ADM_full23_assembly_audited_here':False,
        'full_evolution_or_boundary_subsidiary_proved_here':False,
        'finite_amplitude_nonlinear_family':False,'cosmological_origin_verified':False}}

MUTATIONS=('wrong_hamiltonian_denominator','scalar_compensator_actual_sign',
 'omit_Pi_cross_term','wrong_scalar_velocity_sign','wrong_hamiltonian_gradient',
 'spacetime_affine_radial_momentum','wrong_spatial_trace_factor',
 'wrong_y_z_derivative_scale')
EXPECTED_FAILURES={
 'wrong_hamiltonian_denominator':'Hamiltonian_scalar_sector_2',
 'scalar_compensator_actual_sign':'Hamiltonian_scalar_sector_0',
 'omit_Pi_cross_term':'physical_radial_momentum_G_minus_Fy_equals_pPi',
 'wrong_scalar_velocity_sign':'actual_profile_vt',
 'wrong_hamiltonian_gradient':'Hamiltonian_D_vs_registered',
 'spacetime_affine_radial_momentum':'second_physical_radial_momentum_0',
 'wrong_spatial_trace_factor':'second_physical_radial_momentum_0',
 'wrong_y_z_derivative_scale':'Hamiltonian_A_unsimplified_vs_registered',
}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--mutation',choices=MUTATIONS);ap.add_argument('--run-controls',action='store_true')
    args=ap.parse_args();report=audit(args.repo_root.resolve(),args.mutation)
    if args.run_controls:
        controls=[]
        with tempfile.TemporaryDirectory(prefix='hdblast-physical-controls-') as tmp:
            for mutation in MUTATIONS:
                output=Path(tmp)/(mutation+'.json')
                command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])
                command += [str(Path(__file__).resolve()),'--repo-root',str(args.repo_root.resolve()),'--output',str(output),'--mutation',mutation]
                result=subprocess.run(command,capture_output=True,text=True,timeout=60)
                reason=next((line.split('ValueError: ',1)[1] for line in reversed(result.stderr.splitlines()) if line.startswith('ValueError: ')),'')
                expected='NONZERO_IDENTITY:'+EXPECTED_FAILURES[mutation]+':terms='
                need(result.returncode!=0 and not output.exists() and reason.startswith(expected),'MUTATION_ESCAPED_OR_UNEXPECTED_FAILURE:'+mutation)
                controls.append({'mutation':mutation,'status':'REJECTED','reason':reason})
        report['wrong_formula_controls']=controls
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'identities':report['identity_count'],'controls':len(report.get('wrong_formula_controls',[]))}))
if __name__=='__main__':main()
