#!/usr/bin/env python3
"""Independent exact arbitrary-eight-mode covector wave/gradient binding."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
PINS={'covector_wave_rows.py':'cb025648deda6eec4cfd77b3bb96c6cf6e0713f2b55987b3fdfff22dbb066bbf',
      'codazzi_laurent.py':'178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2',
      'angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9'}
CONTROLS=('omit_radial_p2','wrong_temporal_Ricci_sign','unnormalized_symgrad',
          'omit_wave_half','omit_derivative_slot','positive_time_contraction',
          'omit_physical_time_rho','frozen_rho','time_connection_sign',
          'omit_radial_warp_connection','omit_FRW_connection',
          'scalar_instead_of_gradient_laplacian','wrong_gradient_time_coupling')


def need(ok,message):
    if not ok:
        raise ValueError(message)


def authenticate():
    records={}
    for name,pin in PINS.items():
        path=HERE/name
        need(path.is_file() and not path.is_symlink(),'HELPER_TYPE:'+name)
        raw=path.read_bytes()
        need(hashlib.sha256(raw).hexdigest()==pin,'HELPER_PIN:'+name)
        records[name]={'sha256':pin,'bytes':len(raw)}
    for path in (Path(__file__).resolve(),*(HERE/name for name in PINS)):
        tree=ast.parse(path.read_text())
        need(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'ASSERT:'+path.name)
        allowed={'argparse','ast','fractions','hashlib','importlib','itertools','json',
                 'math','pathlib','subprocess','sys','tempfile','functools'}
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):
                need(all(alias.name.split('.')[0] in allowed for alias in node.names),'IMPORT_CLOSURE:'+path.name)
            if isinstance(node,ast.ImportFrom):
                need(node.module and node.module.split('.')[0] in allowed,'IMPORT_CLOSURE:'+path.name)
            if isinstance(node,ast.Call) and isinstance(node.func,ast.Name):
                need(node.func.id not in {'eval','exec','compile','__import__'},'DYNAMIC_CODE:'+path.name)
    return records


def load(name):
    spec=importlib.util.spec_from_file_location('wave_'+name,HERE/(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calculate(mutation=None):
    before=authenticate()
    wave_rows=load('covector_wave_rows').wave_rows
    alg=load('codazzi_laurent');sphere=load('angular_generic');P=alg.P
    I=range(6);J=range(4);sign=(1,-1,1,1,1,1)
    rho,c,a,theta,p=[P.var(n) for n in ('rho','c','a','theta','p')]
    ay=-rho**-2-p*p/3
    potential=6/rho**2+p*p/2-6*a*a
    hubble=theta/rho
    names=[prefix+str(l) for prefix in ('f','g') for l in (0,2,4)]+['h2','h4']
    coefficients={name:P.var(name) for name in names}
    def derivative(x,axis):
        x=P.cv(x)
        rules={'rho':a*rho,'a':ay} if axis==0 else {'c':c*theta,'theta':1-theta*theta}
        if mutation=='frozen_rho' and axis==0:rules['rho']=P()
        result=P()
        for name in x.names():
            if name.split('_')[0] in coefficients:
                base,_,suffix=name.partition('_')
                new=suffix+('y' if axis==0 else 't')
                need(len(new)<=2,'JET_ORDER:'+name)
                value=P.var(base+'_'+'y'*new.count('y')+'t'*new.count('t'))
            else:value=rules.get(name,P())
            result+=x.diff(name)*value
        return result
    def quotient(x):
        out=P()
        for mon,co in P.cv(x).t.items():
            powers=dict(mon);n=powers.pop('theta',0)
            need(n>=0,'NEGATIVE_THETA')
            for k in range(n//2+1):
                z=dict(powers);z['theta']=n%2;z['c']=z.get('c',0)-2*k
                out+=P({tuple(sorted((key,e) for key,e in z.items() if e)):co*comb(n//2,k)*(-1)**k})
        return out
    add=sphere.add;scale=sphere.scale;mul=sphere.mul
    def ordinary(poly,axis):
        return sphere.clean({e:derivative(v,axis) for e,v in poly.items()})
    gamma={(u,d,v):{} for u,d,v in product(I,repeat=3)}
    gamma[1,1,0]={sphere.ZERO:a}
    gamma[0,1,1]={sphere.ZERO:a}
    for i,k in product(J,repeat=2):
        projector=sphere.P[i,k]
        gamma[k+2,i+2,0]=scale(a,projector)
        gamma[0,i+2,k+2]=scale(-a,projector)
        gamma[k+2,i+2,1]=scale(hubble,projector)
        gamma[1,i+2,k+2]=scale(hubble,projector)
    if mutation=='time_connection_sign':gamma[0,1,1]={sphere.ZERO:-a}
    if mutation=='omit_radial_warp_connection':
        for key in gamma:
            if key==(1,1,0) or key==(0,1,1) or (key[1]>=2 and (key[0]==0 or key[2]==0)):
                gamma[key]={}
    if mutation=='omit_FRW_connection':
        for u,d,v in gamma:
            if d>=2 and (u==1 or v==1):gamma[u,d,v]={}
    def cov(t,indices,d):
        if d==0:out=ordinary(t[indices],0)
        elif d==1:
            factor=1 if mutation=='omit_physical_time_rho' else 1/rho
            out=scale(factor,ordinary(t[indices],1))
        else:
            k=d-2
            out=add(*(mul(sphere.P[k,v],sphere.diff(t[indices],v)) for v in J))
            for slot,index in enumerate(indices):
                if index>=2:
                    replacement=list(indices);replacement[slot]=d
                    out=add(out,mul(sphere.X[index-2],t[tuple(replacement)]))
            out=scale(1/(rho*c),out)
        for slot,index in enumerate(indices):
            if mutation=='omit_derivative_slot' and len(indices)==2 and slot==0:
                continue
            for upper in I:
                replacement=list(indices);replacement[slot]=upper
                out=add(out,scale(-1,mul(gamma[upper,d,index],t[tuple(replacement)])))
        return out
    z={i:{} for i in I}
    for j,l in enumerate((0,2,4)):
        z[0]=add(z[0],scale(coefficients['f'+str(l)],sphere.SCALARS[j][()]))
        z[1]=add(z[1],scale(coefficients['g'+str(l)],sphere.SCALARS[j][()]))
    for j,l in enumerate((2,4)):
        for i in J:z[i+2]=add(z[i+2],scale(coefficients['h'+str(l)],sphere.GRADS[j][i,]))
    tensor_z={(b,):x for b,x in z.items()}
    dz={(d,b):cov(tensor_z,(b,),d) for d,b in product(I,repeat=2)}
    dz2={(d,e,b):cov(dz,(e,b),d) for d,e,b in product(I,repeat=3)}
    signature=(1,1,1,1,1,1) if mutation=='positive_time_contraction' else sign
    rough={b:add(*(scale(signature[d],dz2[d,d,b]) for d in I)) for b in I}
    ricci={b:scale(F(2,3)*potential,z[b]) for b in I}
    if mutation!='omit_radial_p2':ricci[0]=add(ricci[0],scale(p*p,z[0]))
    if mutation=='wrong_temporal_Ricci_sign':ricci[1]=scale(-F(2,3)*potential,z[1])
    wave={b:add(rough[b],ricci[b]) for b in I}
    symfactor=1 if mutation=='unnormalized_symgrad' else F(1,2)
    sym={(i,j):scale(symfactor,add(dz[i,j],dz[j,i])) for i,j in product(I,repeat=2)}
    b_sym={b:add(*(scale(sign[d],add(cov(sym,(d,b),d),scale(-F(1,2),cov(sym,(d,d),b)))) for d in I)) for b in I}
    checks=[]
    def equal_polys(label,left,right):
        residual=add(left,scale(-1,right))
        for exp,x in residual.items():
            rem=quotient(x)
            need(not rem,'WAVE_IDENTITY_NONZERO:'+label+':'+str(exp)+':'+json.dumps(rem.encode()[:3],separators=(',',':')))
        checks.append({'name':label,'left_angular_terms':len(left),'right_angular_terms':len(right),'residual_terms':0})
    targetfactor=1 if mutation=='omit_wave_half' else F(1,2)
    for b in I:equal_polys('B0_SYMGRAD_VS_DIRECT_WAVE:'+str(b),b_sym[b],scale(targetfactor,wave[b]))
    def project(t,basis,label):
        coeff=[]
        for base in basis:
            numerator=sum((sphere.integrate(mul(t[k],v)) for k,v in base.items()),P())
            coeff.append(quotient(numerator/sphere.inner(base,base)))
        for k,value in t.items():
            reconstruct=add(*(scale(cf,base[k]) for cf,base in zip(coeff,basis)))
            equal_polys(label+':'+str(k),value,reconstruct)
        return coeff
    rows={'y':project({():wave[0]},sphere.SCALARS,'WAVE_y_RECONSTRUCTION'),
          '0':project({():wave[1]},sphere.SCALARS,'WAVE_0_RECONSTRUCTION'),
          'a':project({(i,):wave[i+2] for i in J},sphere.GRADS,'WAVE_a_RECONSTRUCTION')}
    inputs={'y':[coefficients['f'+str(l)] for l in (0,2,4)],
            '0':[coefficients['g'+str(l)] for l in (0,2,4)],
            'a':[coefficients['h'+str(l)] for l in (2,4)]}
    formula=wave_rows(inputs,derivative,rho,a,p,c,theta)
    if mutation=='scalar_instead_of_gradient_laplacian':
        formula['a'][0]-=2*coefficients['h2']/(rho**2*c**2)
    if mutation=='wrong_gradient_time_coupling':
        formula['a'][0]+=4*theta*coefficients['g2']/(rho**2*c)
    for field in ('y','0','a'):
        for index,(left,right) in enumerate(zip(rows[field],formula[field])):
            rem=quotient(left-right)
            need(not rem,'WAVE_MODAL_NONZERO:'+field+':'+str(index)+':'+json.dumps(rem.encode()[:3],separators=(',',':')))
            checks.append({'name':'DIRECT_WAVE_VS_EXPORTED_MODAL:'+field+':'+str(index),'residual_terms':0})
    # Modal inner products distinguish the eight coefficients uniquely.
    for label,basis in (('scalar',sphere.SCALARS),('gradient',sphere.GRADS)):
        for i,j in product(range(len(basis)),repeat=2):
            inner=sphere.inner(basis[i],basis[j])
            need(bool(inner) if i==j else not inner,'GRAM:'+label+':'+str((i,j)))
    need(authenticate()==before,'SOURCE_CLOSURE_CHANGED')
    return {'status':'PASS_ARBITRARY_EIGHT_MODE_COVECTOR_WAVE',
        'helper_authentication':before,'checks':checks,'modal_rows':{k:[x.encode() for x in v] for k,v in rows.items()},
        'ordinary_covector_profiles':8,'ordinary_covector_jets_through_order_two':48,
        'identity':'B0(sym nabla Z)=(1/2)(Box0+Ric0)Z',
        'wave_derivation':'Direct sum of second covariant derivatives with every tensor/derivative slot; wave not defined by the Bianchi identity.',
        'physical_dimension':5,'ambient_sphere_slots':4,
        'Ricci_endomorphism':'(2/3)V Z_b+p² delta_b,y Z_y',
        'background_relations':'rho_y=a rho; a_y=-rho^-2-p²/3; V=6rho^-2+p²/2-6a²; c_t=c theta; theta_t=1-theta²; theta²=1-c^-2',
        'scope':'Finite exact operator binding in retained scalar/gradient modes. Existing analytic subsidiary proof and actual smooth field assumptions remain; no new PDE existence, numerical accuracy or cosmological claim.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true')
    args=ap.parse_args()
    need(not args.output.exists(),'REFUSE_OVERWRITE')
    report=calculate(args.mutation)
    if args.run_controls:
        rejected=[]
        with tempfile.TemporaryDirectory(prefix='covector-wave-control-') as temp:
            for mutation in CONTROLS:
                output=Path(temp)/(mutation+'.json')
                command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])
                command += [str(Path(__file__).resolve()),'--output',str(output),'--mutation',mutation]
                process=subprocess.run(command,capture_output=True,text=True,timeout=120)
                lines=process.stderr.strip().splitlines()
                need(process.returncode!=0 and not output.exists() and lines and
                     (lines[-1].startswith('ValueError: WAVE_IDENTITY_NONZERO:') or
                      lines[-1].startswith('ValueError: WAVE_MODAL_NONZERO:')),
                     'WRONG_FORMULA_NOT_MATHEMATICALLY_REJECTED:'+mutation+':'+process.stderr[-3000:])
                rejected.append({'mutation':mutation,'returncode':process.returncode,
                    'successful_output_absent':True,'exact_nonzero_witness':lines[-1].split('ValueError: ',1)[1]})
        report['rejected_corrupted_executions']=rejected
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'modal_rows':8,'checks':len(report['checks'])}))


if __name__=='__main__':main()
