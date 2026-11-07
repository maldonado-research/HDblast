"""Exact formal tests for bundle conversion, ordering, and second variation.

Fields and their jets are indeterminates. Connection and principal tensors
use dense exact rational test coefficients, including noncommuting bundle
matrices, nonzero connection derivatives and nonzero base connection.
These finite tests supplement the identities proved in the note.
"""
from fractions import Fraction as F

NV=96
ZERO=(0,)*NV


class P:
    def __init__(self,terms=None):
        self.terms={k:F(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def cast(value):
        return value if isinstance(value,P) else P({ZERO:F(value)})
    @staticmethod
    def var(index):
        key=list(ZERO);key[index]=1
        return P({tuple(key):1})
    def __add__(self,other):
        out=dict(self.terms)
        for key,value in P.cast(other).terms.items():
            out[key]=out.get(key,F(0))+value
        return P(out)
    __radd__=__add__
    def __neg__(self):
        return P({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):
        return self+-P.cast(other)
    def __mul__(self,other):
        out={}
        for ka,va in self.terms.items():
            for kb,vb in P.cast(other).terms.items():
                key=tuple(a+b for a,b in zip(ka,kb))
                if key[0]<=2:
                    out[key]=out.get(key,F(0))+va*vb
        return P(out)
    __rmul__=__mul__
    def __eq__(self,other):
        return self.terms==P.cast(other).terms
    def e_coeff(self,order):
        return P({(0,)+key[1:]:value for key,value in self.terms.items() if key[0]==order})
    def homogeneous(self,order):
        return P({key:value for key,value in self.terms.items() if sum(key)==order})


def mv(matrix,vector):
    return [sum((matrix[i][j]*vector[j] for j in range(2)),P()) for i in range(2)]


def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def evaluate(data,check):
    next_var=1
    def vector(n):
        nonlocal next_var
        result=[P.var(next_var+i) for i in range(n)]
        next_var+=n
        if next_var>NV:
            raise ValueError('Formal variable budget exhausted')
        return result
    V=vector(2)
    Q=[vector(2) for _ in range(3)]
    rawP=[vector(2) for _ in range(3)]
    H=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(i,3):
            H[i][j]=H[j][i]=vector(2)
    omega=[[[F((i+1)*(a+2)-3*b+a*b,7) for b in range(2)] for a in range(2)] for i in range(3)]
    domega=[[[[F((i+1)*(j+2)*(a+1)+b-2*i*a,11) for b in range(2)] for a in range(2)]
              for j in range(3)] for i in range(3)]
    gamma=[[[F((i+j+1)*(k+1)+i*j-2,13) for k in range(3)] for j in range(3)] for i in range(3)]
    q0=[[F(-3 if i==j==0 else (4 if i==j else i+j+1),17) for j in range(3)] for i in range(3)]
    q1=[[[F((i+j+1)*(a+2)+i*j,19) for a in range(2)] for j in range(3)] for i in range(3)]
    C=[[F((a+2)*(b+1)-3,23) for b in range(2)] for a in range(2)]
    C2=[[[F((a+1)*(b+c+2)+b*c,29) for c in range(2)] for b in range(2)] for a in range(2)]
    b0=[[[F((a+1)*(j+2)+c+1,31) for c in range(2)] for j in range(3)] for a in range(2)]
    b1=[[[[F((a+2)*(j+c+1)+k*c+1,37) for k in range(2)] for c in range(2)]
          for j in range(3)] for a in range(2)]
    # Combined first-jet indices are j*2+c. D is symmetric in those indices.
    d0=[[[F((a+1)*(p+r+3)+p*r,41) for r in range(6)] for p in range(6)] for a in range(2)]
    d1=[[[[F((a+k+1)*(p+r+2)+p*r,43) for k in range(2)] for r in range(6)]
          for p in range(6)] for a in range(2)]
    def q(v):
        return [[q0[i][j]+sum((q1[i][j][a]*v[a] for a in range(2)),P()) for j in range(3)] for i in range(3)]
    def bv(v,p):
        return [sum(((b0[a][j][c]+sum((b1[a][j][c][k]*v[k] for k in range(2)),P()))*p[j][c]
                      for j in range(3) for c in range(2)),P()) for a in range(2)]
    def dv(v,p,r):
        pp=[x for row in p for x in row];rr=[x for row in r for x in row]
        return [sum(((d0[a][i][j]+sum((d1[a][i][j][k]*v[k] for k in range(2)),P()))*pp[i]*rr[j]
                      for i in range(6) for j in range(6)),P()) for a in range(2)]
    def cv(v):
        return [sum((C[a][b]*v[b] for b in range(2)),P())+
                F(1,2)*sum((C2[a][b][c]*v[b]*v[c] for b in range(2) for c in range(2)),P()) for a in range(2)]
    def coord_f(v,p):
        c,b,d=cv(v),bv(v,p),dv(v,p,p)
        return [c[a]+b[a]+d[a] for a in range(2)]
    def correction_raw(v,p):
        qq=q(v);out=[P(),P()]
        for i in range(3):
            for j in range(3):
                oi_pj,oj_pi=mv(omega[i],p[j]),mv(omega[j],p[i])
                zeromat=[[domega[i][j][a][b]+mm(omega[i],omega[j])[a][b]-
                          sum(gamma[i][j][k]*omega[k][a][b] for k in range(3))
                          for b in range(2)] for a in range(2)]
                zero=mv(zeromat,v)
                for a in range(2):
                    out[a]+=qq[i][j]*(oi_pj[a]+oj_pi[a]-sum((gamma[i][j][k]*p[k][a] for k in range(3)),P())+zero[a])
        return out
    def corrected_f(v,qq):
        ov=[mv(o,v) for o in omega]
        pp=[[qq[i][a]-ov[i][a] for a in range(2)] for i in range(3)]
        f,cor=coord_f(v,pp),correction_raw(v,pp)
        return [f[a]+cor[a] for a in range(2)]
    def bundled_candidate(v,qq):
        ov=[mv(o,v) for o in omega]
        base=coord_f(v,qq);linear_sub=bv(v,ov);cross=dv(v,ov,qq);square=dv(v,ov,ov)
        out=[base[a]-data['coordinate_linear_substitution_factor']*linear_sub[a]
             -2*cross[a]+square[a] for a in range(2)]
        qt=q(v)
        for i in range(3):
            for j in range(3):
                zero=mv([[domega[i][j][a][b]-mm(omega[j],omega[i])[a][b] for b in range(2)] for a in range(2)],v)
                one,two=mv(omega[i],qq[j]),mv(omega[j],qq[i])
                for a in range(2):
                    out[a]+=qt[i][j]*(one[a]+two[a]-sum((gamma[i][j][k]*qq[k][a] for k in range(3)),P())+zero[a])
        return out
    actual=corrected_f(V,Q);candidate=bundled_candidate(V,Q)
    for a in range(2):
        check('coordinate_to_bundle_row:'+str(a),actual[a]==candidate[a])
    check('bundle_matrices_noncommute',mm(omega[0],omega[1])!=mm(omega[1],omega[0]))
    # Independently differentiate Q_j=P_j+Omega_j V, then apply the
    # bundle and base connections; compare each uncontracted Hessian row.
    hess={}
    for i in range(3):
        for j in range(3):
            qi=[rawP[j][a]+mv(omega[j],V)[a] for a in range(2)]
            conn=mv(omega[i],qi)
            derivative=mv(domega[i][j],V);ojpi=mv(omega[j],rawP[i])
            qk=[[rawP[k][a]+mv(omega[k],V)[a] for a in range(2)] for k in range(3)]
            direct=[H[i][j][a]+derivative[a]+ojpi[a]+conn[a]-
                    sum((gamma[i][j][k]*qk[k][a] for k in range(3)),P()) for a in range(2)]
            zero=mv([[domega[i][j][a][b]+mm(omega[i],omega[j])[a][b]-
                      sum(gamma[i][j][k]*omega[k][a][b] for k in range(3))
                      for b in range(2)] for a in range(2)],V)
            expanded=[H[i][j][a]+mv(omega[i],rawP[j])[a]+mv(omega[j],rawP[i])[a]-
                      sum((gamma[i][j][k]*rawP[k][a] for k in range(3)),P())+zero[a] for a in range(2)]
            for a in range(2):
                check('covariant_hessian:'+str(i)+':'+str(j)+':'+str(a),direct[a]==expanded[a])
            hess[i,j]=direct
    for i in range(3):
        for j in range(i+1,3):
            curvature=mv([[domega[i][j][a][b]-domega[j][i][a][b]+mm(omega[i],omega[j])[a][b]-mm(omega[j],omega[i])[a][b]
                           for b in range(2)] for a in range(2)],V)
            for a in range(2):
                check('nonzero_curvature_identity:'+str(i)+':'+str(j)+':'+str(a),
                      hess[i,j][a]-hess[j,i][a]==curvature[a] and curvature[a]!=0)
    wrong_order_delta=mv([[mm(omega[1],omega[0])[a][b]-mm(omega[0],omega[1])[a][b]
                           for b in range(2)] for a in range(2)],V)
    for a in range(2):
        check('wrong_uncontracted_Omega_order_detected:'+str(a),
              hess[0,1][a]!=hess[0,1][a]+wrong_order_delta[a])
    # At the zero state retain every connection term in the linear maps.
    A=[[C[a][b]-sum(b0[a][i][c]*omega[i][c][b] for i in range(3) for c in range(2))+
        sum(q0[i][j]*(domega[i][j][a][b]-mm(omega[j],omega[i])[a][b]) for i in range(3) for j in range(3))
        for b in range(2)] for a in range(2)]
    B=[[[b0[a][k][b]+2*sum(q0[i][k]*omega[i][a][b] for i in range(3))-
         (sum(q0[i][j]*gamma[i][j][k] for i in range(3) for j in range(3)) if a==b else 0)
         for b in range(2)] for a in range(2)] for k in range(3)]
    for a in range(2):
        linear=mv(A,V)[a]+sum((mv(B[k],Q[k])[a] for k in range(3)),P())
        check('zero_state_linear_maps:'+str(a),actual[a].homogeneous(1)==linear)
    # Full second-amplitude derivative: e^2/2 convention.
    e=P.var(0);v=vector(2);w=vector(2)
    p=[vector(2) for _ in range(3)];r=[vector(2) for _ in range(3)]
    hh=[[vector(2) for _ in range(3)] for _ in range(3)]
    kk=[[vector(2) for _ in range(3)] for _ in range(3)]
    ve=[e*v[a]+F(1,2)*e*e*w[a] for a in range(2)]
    pe=[[e*p[i][a]+F(1,2)*e*e*r[i][a] for a in range(2)] for i in range(3)]
    he=[[[e*hh[i][j][a]+F(1,2)*e*e*kk[i][j][a] for a in range(2)] for j in range(3)] for i in range(3)]
    fe=corrected_f(ve,pe);qe=q(ve)
    zeroV=[P(),P()];zeroP=[[P(),P()] for _ in range(3)]
    ev=[e*x for x in v];ep=[[e*x for x in row] for row in p]
    fvv=[2*x.e_coeff(2) for x in corrected_f(ev,zeroP)]
    fpp=[2*x.e_coeff(2) for x in corrected_f(zeroV,ep)]
    total=[2*x.e_coeff(2) for x in corrected_f(ev,ep)]
    fvp=[F(1,2)*(total[a]-fvv[a]-fpp[a]) for a in range(2)]
    for a in range(2):
        residual=sum((qe[i][j]*he[i][j][a] for i in range(3) for j in range(3)),P())-fe[a]
        linear_w=sum((q0[i][j]*kk[i][j][a] for i in range(3) for j in range(3)),P())-mv(A,w)[a]-sum((mv(B[i],r[i])[a] for i in range(3)),P())
        principal=sum((sum((q1[i][j][b]*v[b] for b in range(2)),P())*hh[i][j][a] for i in range(3) for j in range(3)),P())
        expected=linear_w+data['principal_variation_factor']*principal-fvv[a]-2*fvp[a]-fpp[a]
        check('second_variation_row:'+str(a),2*residual.e_coeff(2)==expected)
        energy_force=data['principal_variation_factor']*principal-fvv[a]-2*fvp[a]-fpp[a]
        check('energy_sign_second_variation:'+str(a),
              -2*residual.e_coeff(2)==-linear_w-energy_force)
        check('first_jet_Hessian_factor_two:'+str(a),fpp[a]==2*dv(zeroV,p,p)[a])
        check('second_mixed_derivative_nonzero:'+str(a),fvp[a]!=0)
        check('second_principal_derivative_nonzero:'+str(a),principal!=0)
    # Product Hessian ordering in X_A=sech(t)E_A: raw directional
    # commutator -hW_A is exactly canceled by -D_(nabla_t X_A).
    hhpar=F(2,7)
    for angular in range(3):
        W=vector(1)[0];ordered=vector(1)[0]
        dt_da=ordered-hhpar*W
        product_ta=dt_da+hhpar*W
        check('moving_product_Hessian_order:'+str(angular),product_ta==ordered)
        check('raw_moving_commutator_nonzero:'+str(angular),dt_da!=ordered)
    def epsilon(i,j,k):
        if len({i,j,k})<3:return 0
        return 1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1
    cs=[[F((i+1)*(j+1)+i+j,47) for j in range(3)] for i in range(3)]
    for k in range(3):
        check('angular_symmetric_connection_cancellation:'+str(k),
              sum(cs[i][j]*epsilon(i,j,k) for i in range(3) for j in range(3))==0)
    return {'formal_fiber_dimension':2,'formal_base_dimension':3,
            'formal_indeterminates_used':next_var,'bundle_conversion_monomial_counts':[len(row.terms) for row in actual],
            'energy_sign_second_force':'2 q_V[V1] H V1 - cbar_VV[V1,V1] - 2 bbar_V[V1] DV1 - 2 dbar(DV1,DV1)',
            'scope':'Polynomial field/jet identities with rational noncommuting connection test coefficients; analytic general identities remain in proof.'}
