import json

def run(repo_root,pins,mutation=None):
    import sympy as S
    y=S.symbols('y')
    rho,a,p,f,xi,Qh,Qm,H=[S.Function(n)(y) for n in ('rho','a','p','f','xi','Qh','Qm','H')]
    V1=S.symbols('V1')
    D=lambda x:S.diff(x,y)
    V=p*p/2-6*a*a+6/rho**2
    rules={D(rho):a*rho,D(a):-rho**-2-p*p/3,D(p):V1-4*a*p}
    
    def reduce(x):
        for k in range(5):x=S.expand(x.doit().subs(rules,simultaneous=True))
        return S.factor(x)
    
    checks=[];controls=[]
    def eq(name,x):
        r=reduce(x)
        if r!=0:raise ValueError(name+': '+str(r))
        checks.append(name)
    def bad(name,x):
        if reduce(x)==0:raise ValueError('Wrong formula accepted: '+name)
        controls.append(name)
    def Pg(A,B):
        return -3*D(D(B))/2+3*a*D(A)/2-6*a*D(B)+(3*D(a)+6*a*a+p*p/2)*A-3*B/rho**2
    def Pm(C,T):return -2*(D(C)+3*a*(C-T))
    C=-3*(D(D(f))+3*a*D(f))/2-V*f
    eq('conformal_particular_operator',Pg(f,f)-C)
    eq('scalar_fixed_gauge_operator',Pg(2*D(xi),2*a*xi)-D(rho**4*p*p*xi)/rho**4)
    Clag=rho**4*C
    Euler=S.diff(Clag,f)-D(S.diff(Clag,D(f)))+D(D(S.diff(Clag,D(D(f)))))
    eq('compact_conformal_integral',Euler+S.Rational(9,2)*rho*rho)
    I,J=S.symbols('I J',nonzero=True)
    eq('bump_coefficient_zero_total_primitive',I-S.Rational(9,2)*(2*I/(9*J))*J)
    xiH=H/(rho**4*p*p)
    eq('full_metric_particular_constraint',
       reduce(Pg(f+2*D(xiH),f+2*a*xiH)+Qh).subs(D(H),-rho**4*(Qh+C)))
    momentum=-Qm/(6*a) if mutation!='wrong_momentum' else Qm/(3*a)
    eq('canonical_momentum_factor',Pm(S.Integer(0),momentum)+Qm)
    bad('suggested_momentum_sign_and_factor',Pm(S.Integer(0),Qm/(3*a))+Qm)
    bad('wrong_bump_coefficient',I-S.Rational(9,2)*(I/(9*J))*J)
    bad('wrong_primitive_sign',C+(Qh+C)+Qh)
    s,u,j=S.symbols('s u j')
    T=S.Integer(0);U=-Qm/(6*a)
    eq('physical_mean_K_radial',4*T/3-2*U+6*u*s-(Qm/(3*a)+6*u*s))
    eq('physical_mean_K_angular',-2*T/3-u*s+u*s)
    eq('coordinate_mean_velocity_radial',2*rho*(4*T/3-2*U+16*u*s/3)-2*rho*(Qm/(3*a)+16*u*s/3))
    eq('coordinate_mean_velocity_angular',2*rho*(-2*T/3-2*u*s/3)+4*rho*u*s/3)
    e,Y=S.symbols('e Y')
    volume=S.sqrt((1-4*e*u*Y)*(1+2*e*u*Y)**3)
    eq('zero_second_coordinate_scalar_velocity',S.diff(rho*(1+e*u*Y)*e*j*Y/volume,e,2).subs(e,0))
    bgrow=next(row for row in pins['sources'] if row['key']=='SB')
    bg=json.loads((repo_root/bgrow['path']).read_text())['background_intervals']
    alo,ahi=map(S.Rational,bg['a']['exact']);plo,phi=map(S.Rational,bg['p']['exact'])
    if not(alo>0 and ahi>alo and plo<phi<0):raise ValueError('Unsafe a or p divisor')
    checks.append('retained_nonzero_a_and_p')
    return {'status':'PASS','checks':checks,'check_count':len(checks),'in_process_controls':controls,'in_process_control_count':len(controls),'real_subprocess_control_count':0,'sympy':S.__version__,'scope':'Explicit preparation algebra only; retained analytic inverse, IFT and evolution remain conditional.'}
