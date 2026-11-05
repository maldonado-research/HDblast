"""Independent exact static-collar and negative-time coefficient majorants."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PINS_SHA = '2446c1f664df99bb8a418ac2c0faccc208ff43963aad3514ff022b8753c3b550'


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def add(*xs):
    return sum(x[0] for x in xs),sum(x[1] for x in xs)


def scale(c,x):
    return (c*x[0],c*x[1]) if c>=0 else (c*x[1],c*x[0])


def mul(x,y):
    v=[a*b for a in x for b in y]
    return min(v),max(v)


def fixed(v):
    return F(v),F(v)


def authenticate():
    raw=(HERE/'SOURCE_PINS.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==PINS_SHA,'SOURCE_PIN_FILE')
    pins=json.loads(raw)['files']
    data={}
    for rel,pin in pins.items():
        path=HERE/rel
        need(not path.is_symlink(),'SYMLINK_INPUT')
        raw=path.read_bytes()
        need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'PIN:'+rel)
        data[Path(rel).name]=json.loads(raw) if rel.endswith('.json') else raw.decode()
    return data,pins


def run():
    data,pins=authenticate()
    data['BACKGROUND_DERIVATIVE_BOUNDS.json']=data['WIDE_BACKGROUND.json']
    data['PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json']=data['SHELL_BOX.json']
    bg=data['BACKGROUND_DERIVATIVE_BOUNDS.json']
    shell=data['PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json']
    ell,Z=F(1,10**6),F(8,10**6)
    need(bg['candidate_state_exact']['Z']==str(Z)
         and bg['candidate_state_exact']['ell']==str(ell),'WIDE_COLLAR_8ELL')
    outer={k:tuple(map(F,v['outer']['exact']))
           for k,v in bg['first_exit']['component_margins'].items()}
    expected={'rho':(F(129),F(131)),'eta':(-F(9,100000),-F(8,100000)),
              'p':(-F(7,50000),-F(3,25000)),
              'a':(F(111,1000),F(14,125))}
    need(outer==expected,'STATIC_OUTER_RECTANGLE')
    rho,eta,p,a=(outer[k] for k in ('rho','eta','p','a'))
    e=(-eta[1],-eta[0])
    # Independent monomial summation on e>0, rather than importing the
    # retained producer or copying its Horner evaluation/flow table.
    u1=add(*(scale(c,(e[0]**n,e[1]**n)) for n,c in
             ((1,-F(28,9)),(2,F(50,9)),(3,F(2,3)),
              (4,-F(20,9)),(5,F(4,9)))))
    rhs={
        'rho':scale(-1,mul((rho[0]**2,rho[1]**2),a)),
        'eta':scale(-1,mul(rho,p)),
        'p':scale(-1,mul(rho,add(u1,scale(-4,mul(a,p))))),
        'a':add((1/rho[1],1/rho[0]),
                scale(F(1,3),mul(rho,(p[1]**2,p[0]**2)))),
    }
    first_exit={}
    for k in expected:
        initial=tuple(map(F,shell['boxes'][k+'_b']['exact']))
        tube=add(initial,mul((F(0),Z),rhs[k]))
        need(outer[k][0]<tube[0]<=tube[1]<outer[k][1],
             'STRICT_WIDE_FIRST_EXIT:'+k)
        first_exit[k]={'shell':list(map(str,initial)),
                       'independent_monomial_rhs':list(map(str,rhs[k])),
                       'flow_tube':list(map(str,tube)),
                       'strict_lower_margin':str(tube[0]-outer[k][0]),
                       'strict_upper_margin':str(outer[k][1]-tube[1])}
    need(131*Z<9,'PROPER_RADIAL_ENDPOINT_MARGIN')

    delta,emin=e[1],e[0]
    pmin,pmax=-p[1],-p[0]
    hlo=(F(28,9)*emin-F(50,9)*delta**2-F(2,3)*delta**3-F(4,9)*delta**5)/pmax
    hhi=F(28,9)*delta/pmin
    need(hlo>F(7,4) and hhi==F(7,3),'H_BOUNDS')
    # U2=28/9+(100/9)eta-2eta²-(80/9)eta³-(20/9)eta⁴.
    u2err=F(100,9)*delta+2*delta**2+F(80,9)*delta**3+F(20,9)*delta**4
    need(F(28,9)-u2err>F(31,10) and F(28,9)+u2err<F(25,8),'U2_BOUNDS')
    qlo=2*F(7,4)**2-10*a[1]*F(7,3)-F(25,8)-pmax**2/6
    qhi=2*F(7,3)**2-10*a[0]*F(7,4)+F(35,4)*a[1]**2-F(31,10)
    need(qlo>F(1,3) and qhi<6,'Q_BOUNDS')
    past_time=ell
    wlo=F(129)**2/3-F(35,4)*past_time**2
    whi=F(131)**2*6+F(35,4)*past_time**2
    need(wlo>5500 and whi<103000,'FRESH_NEGATIVE_TIME_POTENTIAL_BOUND')
    need(ell**2*103000==F(103,10**9),'SCALED_POTENTIAL')
    mu=F(103,10**9)
    tau=F(1)
    beta=mu*tau*tau/2
    B=(1+tau*ell)/(1-beta)
    kappa=ell/2+mu*tau*B/2
    need(0<beta<1 and B<F(1000002,1000000),'UNIFORM_VOLterra_C0_L2_BOUND')
    need(0<kappa<F(1,10**6),'TWO_NONZERO_BRANCHES')
    share_lo=(1-kappa)**2/((1-kappa)**2+(1+kappa)**2)
    share_hi=(1+kappa)**2/((1-kappa)**2+(1+kappa)**2)
    need(share_lo>F(499999,10**6) and share_hi<F(500001,10**6),'BOTH_PRINCIPAL_WEIGHT_SHARES')
    need(1-kappa>F(999999,10**6) and 1+kappa<F(1000001,10**6),'BRANCH_NORM_BOUNDS')
    quarter=F(1,4)
    need(1-quarter==F(3,4) and 2+quarter==F(9,4),'QUARTER_SUPPORT')
    need(3*ell<Z,'WHOLE_OPEN_PRESHELL_SUPPORT_INSIDE_STATIC_COLLAR')
    # Exact Taylor2 plus a geometric tail proves exp(4/3)<4.
    arg=F(4,3)
    exp_upper=1+arg+arg**2/2+(arg**3/6)/(1-arg/4)
    need(exp_upper==F(103,27) and exp_upper<4,'DIAMOND_BUMP_LOWER_BOUND')
    beta_half=mu/8
    B_half=(1+ell/2)/(1-beta_half)
    face_lower=F(1,2)-ell/4-beta_half*B_half
    apex_upper=-ell*(1-2*ell)/16+beta_half*B_half
    need(face_lower>F(499999,10**6),'TWO_PAST_FACE_MIDPOINTS_POSITIVE')
    need(apex_upper<-F(1,25000000),'COMMON_PAST_APEX_NEGATIVE')
    goursat_q=mu/16
    need(goursat_q<1,'CHARACTERISTIC_DIAMOND_VOLterra_CONTRACTION')
    # Check the constant and all radial coefficients in the conjugated
    # potential after inserting the explicitly derived C_z formula.
    cz={'h2':F(2),'ah':F(-10),'a2':F(5),'U2':F(-2),'p2':-F(5,3)}
    c_square_quarter={'h2':F(1),'ah':F(-5),'a2':F(25,4),'U2':F(0),'p2':F(0)}
    combined={k:cz[k]/2+c_square_quarter[k]+(F(2,3) if k=='p2' else F(0)) for k in cz}
    need(combined=={'h2':F(2),'ah':F(-10),'a2':F(35,4),'U2':F(-1),'p2':-F(1,6)},'MASTER_CONJUGATION_RADIAL_COEFFICIENTS')
    need(F(35,4)-F(25,4)-F(5,2)==0,'MASTER_CONJUGATION_CONSTANT_CANCELLATION')
    need(authenticate()[1]==pins,'INPUT_AUTHENTICATION_AFTER')
    result={'status':'PASS_ACTUAL_REGISTERED_PAST_MASTER_PRECONTACT',
            'source_pins_sha256':PINS_SHA,'selected_local_inputs_authenticated':len(pins),'static_collar_length':str(Z),
            'first_exit':first_exit,'negative_time_interval':[str(-past_time),'0'],
            'radial_target_scaled':['0','3'],
            'concrete_quarter_slab_radial_target_scaled':['3/4','9/4'],
            'derived_h_interval':['7/4','7/3'],'derived_U2_interval':['31/10','25/8'],
            'derived_Q_interval':['1/3','6'],
            'independent_Q_majorants':[str(qlo),str(qhi)],
            'fresh_negative_time_W_interval':['5500','103000'],
            'fresh_scaled_potential_upper':'103/1000000000',
            'temporal_input':'All-real |tanh(t)|<=|t| and W=rho²Q-(35/4)tanh²(t), derived in the accompanying proof; no positive-time potential cap is transferred.',
            'analytic_premises_retained':['Actual smooth radial branch and authenticated shell enclosure.','Standard continuation/first-exit principle and direct negative-time background extension.'],
            'uniform_past_scaled_time':{'interval':'-1<x<=0','closed_slabs':'Every -tau<=x<=0 with 0<tau<1; tau=1 is used only as an upper bound in the constants.','support':'[1-|x|,2+|x|], contained in (0,3) for -1<x<=0.','quarter_support':['3/4','9/4'],'shell_contact':'None on the open precontact interval; no endpoint interaction assertion.'},
            'uniform_estimate':{'scaled_potential_upper':str(mu),'volterra_contraction_upper':str(beta),'C0_and_relative_L2_upper':str(B),'characteristic_relative_error_upper':str(mu*B/2),'kappa':str(kappa),'branch_norm_lower_factor':str(1-kappa),'branch_norm_upper_factor':str(1+kappa),'principal_share_lower':str(share_lo),'principal_share_upper':str(share_hi),'safe_strict_share_interval':['499999/1000000','500001/1000000'],'D_definition':'D=||Q_r||_L2; D^2>=4 for the unit bump.','norm_convention':'Jplus=U_x+U_r travels toward the shell in future time; Jminus=U_x-U_r travels away.'},
            'constructive_past_diamond':{'past_apex':['-1/2','3/2'],'null_coordinates':'u=x+r, v=x-r','domain':'1<=u<=2,-2<=v<=-1,u+v<=0','data_faces':['u=1','v=-2'],'common_corner_condition':'f(-2)=g(1)','exp_4_over_3_upper':str(exp_upper),'B_half':str(B_half),'C0_free_error_upper':str(beta_half*B_half),'two_face_midpoint_lower':str(face_lower),'common_apex_upper':str(apex_upper),'goursat_contraction_upper':str(goursat_q),'conclusion':'Both actual compatible characteristic traces are nonconstant and change sign; this reconstructs the selected linear Cauchy solution, not a physical preparation law.'},
            'physical_scope':['The master reconstructs the coupled first linearized Einstein-scalar perturbation by retained differential identities; it is not an arbitrary test scalar.','Harmonic and exact-family identification retain full one-jet and local development/uniqueness premises; no effective nonlinear parameter interval is supplied.','The normalized-master principal weights are not physical energy or radiation flux, and their purity obstruction is limited to this finite past interval.','No second-order constrained state, past incoming selection, physical source, absolute proper observable, SI scale or cosmological origin is obtained.'],
            'activity':{'project_imports':0,'legacy_checker_invocations':0,'background_oracle_calls':0,'source_oracle_calls':0,'bump_queries':0,'numerical_PDE_solves':0,'exact_background_ODE_image_recomputations':4}}
    return result



def main():
    class Parser(argparse.ArgumentParser):
        def error(self,message):
            raise ValueError('CLI_ARGUMENTS')
    parser=Parser(add_help=False)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    need(args.output.is_absolute(),'ABSOLUTE_OUTPUT_REQUIRED')
    need(not any(p.is_symlink() for p in (args.output,*args.output.parents)),'SYMLINK_OUTPUT')
    output=args.output.resolve()
    need(not output.exists(),'OUTPUT_EXISTS')
    need(not output.is_relative_to(HERE),'OUTPUT_INSIDE_SOURCE')
    need(output.parent.is_dir(),'OUTPUT_PARENT_MISSING')
    result=run()
    raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    with output.open('xb') as stream:
        stream.write(raw)
    print(json.dumps({'status':result['status'],'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':
    try:
        main()
    except (ValueError,OSError,ArithmeticError,KeyError,TypeError,IndexError,AttributeError,UnicodeError) as exc:
        witness=str(exc) if isinstance(exc,ValueError) else type(exc).__name__
        print(json.dumps({'status':'REJECTED','witness':witness},sort_keys=True))
        raise SystemExit(2)
