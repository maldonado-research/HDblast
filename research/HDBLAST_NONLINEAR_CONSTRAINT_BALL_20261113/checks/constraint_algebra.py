#!/usr/bin/env python3
"""Finite exact nonlinear constraint Hessian witnesses, not a ball estimate.

The independent oracle differentiates rational two-parameter metric/canonical
families and uses chart Christoffel/Ricci and density divergence definitions.
Production uses background-covariant curvature and full canonical Hessians.
No analytic norm, positivity radius, K, or contraction certificate is inferred.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import sympy as S
import mpmath

e,f=S.symbols("e f")
FACTORS=("inverse_mixed_factor","det_trace_product_factor","det_trace_cross_factor",
 "connection_mixed_factor","connection_state_inverse_mixed_factor","Ricci_quadratic_cross_factor",
 "Ricci_state_connection_cross_factor","scalar_curvature_inverse_cross_factor",
 "scalar_curvature_inverse_mixed_factor","Hamiltonian_curvature_factor","Hamiltonian_momentum_factor",
 "Hamiltonian_trace_squared_factor","Hamiltonian_scalar_density_factor","Hamiltonian_scalar_gradient_factor",
 "Hamiltonian_potential_factor","momentum_h_dPi_factor","momentum_Pi_dh_factor",
 "momentum_Pi_trace_dh_factor","momentum_matter_factor")
F={k:S.Symbol("claim_"+k) for k in FACTORS}


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def clean(v):return S.cancel(S.expand(v))
def trunc(v):
    P=S.Poly(S.expand(v),e,f)
    return S.Add(*(a*e**i*f**j for (i,j),a in P.terms() if i<=1 and j<=1))
def cross(v):return clean(S.diff(v,e,f).subs({e:0,f:0}))
def scalar(M):return M[0]
def contract(A,B):return sum(A[i,j]*B[i,j] for i in range(4) for j in range(4))


def build_checks():
    coords=S.symbols("x0:4");t,x,y,z=coords;n=4;at={a:0 for a in coords}
    g0=S.diag(4,1+t+t*t,1,1)
    h=S.Matrix([[1+t+x*x,1+t*y,2+t+z,x*y],
       [1+t*y,2+x+y*y,y+t*z,1+t*z],
       [2+t+z,y+t*z,-1+x*z,1+y],
       [x*y,1+t*z,1+y,1+t+x*x]])
    k=S.Matrix([[2+t*y,t+x,1+y*z,1+z],
       [t+x,-1+z*z,1+t+y,x],
       [1+y*z,1+t+y,1+x+y*z,2+t],
       [1+z,x,2+t,-2+y]])
    A=g0.inv();G=g0.subs(at);H=h.subs(at);K=k.subs(at);A0=A.subs(at)
    ge=g0+e*h+f*k;gp=ge.subs(at)
    inverse_exact=gp.inv()
    inv=inverse_exact.applyfunc(lambda a:clean(a.subs({e:0,f:0}))+clean(S.diff(a,e).subs({e:0,f:0}))*e+
          clean(S.diff(a,f).subs({e:0,f:0}))*f+cross(a)*e*f)
    dg=[ge.diff(a).subs(at) for a in coords]
    ddg=[[ge.diff(a).diff(b).subs(at) for b in coords] for a in coords]
    dinv=[(-inv*dg[a]*inv).applyfunc(trunc) for a in range(n)]
    def W(l,i,j):return dg[i][l,j]+dg[j][l,i]-dg[l][i,j]
    gamma=[[[trunc(sum(inv[a,l]*W(l,i,j) for l in range(n))/2) for j in range(n)] for i in range(n)] for a in range(n)]
    def dgamma(c,a,i,j):
        return trunc(sum(dinv[c][a,l]*W(l,i,j)+inv[a,l]*(ddg[c][i][l,j]+ddg[c][j][l,i]-ddg[c][l][i,j]) for l in range(n))/2)
    ric=S.Matrix(n,n,lambda i,j:trunc(sum(dgamma(a,a,i,j)-dgamma(j,a,i,a) for a in range(n))+
       sum(gamma[a][a][b]*gamma[b][i][j]-gamma[a][j][b]*gamma[b][i][a] for a in range(n) for b in range(n))))
    R_oracle=trunc(contract(inv,ric))
    Gamma0=[[[clean(sum(A[a,l]*(S.diff(g0[l,j],coords[i])+S.diff(g0[l,i],coords[j])-S.diff(g0[i,j],coords[l])) for l in range(n))/2)
              for j in range(n)] for i in range(n)] for a in range(n)]
    gam0=[[[v.subs(at) for v in row] for row in section] for section in Gamma0]
    def C(direction):
        def cov(a,i,j):return S.diff(direction[i,j],coords[a])-sum(Gamma0[l][a][i]*direction[l,j]+Gamma0[l][a][j]*direction[i,l] for l in range(n))
        return [[[clean(sum(A[a,l]*(cov(i,l,j)+cov(j,l,i)-cov(l,i,j)) for l in range(n))/2)
                  for j in range(n)] for i in range(n)] for a in range(n)]
    Ch=C(h);Ck=C(k)
    def tensor_point(T):return [[[v.subs(at) for v in row] for row in section] for section in T]
    ch=tensor_point(Ch);ck=tensor_point(Ck)
    dh=[tensor_point([[[S.diff(v,q) for v in row] for row in section] for section in Ch]) for q in coords]
    dk=[tensor_point([[[S.diff(v,q) for v in row] for row in section] for section in Ck]) for q in coords]
    Ah=-A0*H*A0;Ak=-A0*K*A0
    Aef=F["inverse_mixed_factor"]*(A0*H*A0*K*A0+A0*K*A0*H*A0)
    # C_hk=−A h C_k−A k C_h, with consistent fixed-background connection.
    chk=[[[F["connection_mixed_factor"]*sum((A0*H)[a,l]*ck[l][i][j]+(A0*K)[a,l]*ch[l][i][j] for l in range(n))
          for j in range(n)] for i in range(n)] for a in range(n)]
    dchk=[]
    for r in range(n):
        dAH=(A.diff(coords[r])*h+A*h.diff(coords[r])).subs(at)
        dAK=(A.diff(coords[r])*k+A*k.diff(coords[r])).subs(at)
        dchk.append([[[F["connection_mixed_factor"]*sum(dAH[a,l]*ck[l][i][j]+(A0*H)[a,l]*dk[r][l][i][j]+
                         dAK[a,l]*ch[l][i][j]+(A0*K)[a,l]*dh[r][l][i][j] for l in range(n))
                      for j in range(n)] for i in range(n)] for a in range(n)])
    def covC(T,dT,r,a,i,j):return dT[r][a][i][j]+sum(gam0[a][r][l]*T[l][i][j]-gam0[l][r][i]*T[a][l][j]-gam0[l][r][j]*T[a][i][l] for l in range(n))
    def L(T,dT):return S.Matrix(n,n,lambda i,j:sum(covC(T,dT,a,a,i,j)-covC(T,dT,j,a,i,a) for a in range(n)))
    def Q(T,U):return S.Matrix(n,n,lambda i,j:sum(T[a][a][b]*U[b][i][j]-T[a][j][b]*U[b][i][a] for a in range(n) for b in range(n)))
    Rich=L(ch,dh);Rick=L(ck,dk)
    Ricmixed=L(chk,dchk)+F["Ricci_quadratic_cross_factor"]*(Q(ch,ck)+Q(ck,ch))
    Ric0=ric.subs({e:0,f:0})
    Rprod=(contract(A0,Ricmixed)+F["scalar_curvature_inverse_cross_factor"]*(contract(Ah,Rick)+contract(Ak,Rich))+
        F["scalar_curvature_inverse_mixed_factor"]*contract(Aef,Ric0))
    # Exact inverse and determinant Hessians use noncommuting symmetric H,K.
    det_exact=gp.det();D0=G.det();th=S.trace(A0*H);tk=S.trace(A0*K);trhk=S.trace(A0*H*A0*K)
    det_prod=D0*(F["det_trace_product_factor"]*th*tk+F["det_trace_cross_factor"]*trhk)
    J=1/D0;Jh=-J*th;Jk=-J*tk;Jhk=J*(th*tk+trhk)
    checks=[("matrix_inverse_mixed_Hessian",list(Aef-inverse_exact.applyfunc(cross))),
       ("matrix_determinant_mixed_Hessian",[det_prod-cross(det_exact)]),
       ("reciprocal_determinant_mixed_Hessian",[Jhk-cross(1/det_exact)]),
       ("connection_mixed_Hessian",[chk[a][i][j]-cross(gamma[a][i][j]) for a in range(n) for i in range(n) for j in range(n)]),
       ("Ricci_mixed_Hessian",list(Ricmixed-ric.applyfunc(cross))),
       ("scalar_curvature_mixed_Hessian",[Rprod-cross(R_oracle)])]
    # A second production route uses a fixed flat chart connection. Its state
    # C0=Gamma(g0) is nonzero. Thus the nonlinear A_hk*S(g0) term and the
    # Ricci products Q(C_hk,C0)+Q(C0,C_hk) cannot disappear by choosing the
    # comparison connection to equal the state connection.
    def ww(dfield,l,i,j):return dfield[i][l,j]+dfield[j][l,i]-dfield[l][i,j]
    dhg=[h.diff(q).subs(at) for q in coords];dkg=[k.diff(q).subs(at) for q in coords]
    d0g=[g0.diff(q).subs(at) for q in coords]
    dd0=[[g0.diff(q).diff(r).subs(at) for r in coords] for q in coords]
    ddh=[[h.diff(q).diff(r).subs(at) for r in coords] for q in coords]
    ddk=[[k.diff(q).diff(r).subs(at) for r in coords] for q in coords]
    def ddw(ddfield,r,l,i,j):return ddfield[r][i][l,j]+ddfield[r][j][l,i]-ddfield[r][l][i,j]
    statefactor=F["connection_state_inverse_mixed_factor"]
    cf=[[[sum(statefactor*Aef[a,l]*ww(d0g,l,i,j)+Ah[a,l]*ww(dkg,l,i,j)+Ak[a,l]*ww(dhg,l,i,j) for l in range(n))/2
           for j in range(n)] for i in range(n)] for a in range(n)]
    dcf=[]
    def prod(matrices):
        answer=S.eye(n)
        for matrix in matrices:answer=answer*matrix
        return answer
    def product_derivative(matrices,derivatives):
        return sum((prod(matrices[:i]+[derivatives[i]]+matrices[i+1:]) for i in range(len(matrices))),S.zeros(n))
    for r in range(n):
        da=A.diff(coords[r]).subs(at)
        dah=-(da*H*A0+A0*dhg[r]*A0+A0*H*da)
        dak=-(da*K*A0+A0*dkg[r]*A0+A0*K*da)
        daef=F["inverse_mixed_factor"]*(product_derivative([A0,H,A0,K,A0],[da,dhg[r],da,dkg[r],da])+
                product_derivative([A0,K,A0,H,A0],[da,dkg[r],da,dhg[r],da]))
        dcf.append([[[sum(statefactor*(daef[a,l]*ww(d0g,l,i,j)+Aef[a,l]*ddw(dd0,r,l,i,j))+
                  dah[a,l]*ww(dkg,l,i,j)+Ah[a,l]*ddw(ddk,r,l,i,j)+dak[a,l]*ww(dhg,l,i,j)+Ak[a,l]*ddw(ddh,r,l,i,j)
                  for l in range(n))/2 for j in range(n)] for i in range(n)] for a in range(n)])
    Ricref=(S.Matrix(n,n,lambda i,j:sum(dcf[a][a][i][j]-dcf[j][a][i][a] for a in range(n)))+
        F["Ricci_state_connection_cross_factor"]*(Q(cf,gam0)+Q(gam0,cf))+Q(ch,ck)+Q(ck,ch))
    checks.extend([("fixed_chart_nonlinear_connection_state_term",[cf[a][i][j]-cross(gamma[a][i][j]) for a in range(n) for i in range(n) for j in range(n)]),
       ("fixed_chart_Ricci_state_connection_products",list(Ricref-ric.applyfunc(cross)))])
    # Every canonical sector is nonzero at the reference germ.
    Pi0=S.Matrix([[2+t,1+x,1+y,z],[1+x,-1+y,2+z,1+t],[1+y,2+z,3+x,1+y],[z,1+t,1+y,1+z]])
    Ph=S.Matrix([[1+x,2+t,y,1+z],[2+t,-1+y,1+x,z],[y,1+x,2+z,1+t],[1+z,z,1+t,-2+x]])
    Pk=S.Matrix([[2+y,1+z,1+t,x],[1+z,1+x,z,2+y],[1+t,z,-1+y,1+x],[x,2+y,1+x,3+t]])
    phi_background_offset=S.Symbol("phi_background_offset",real=True)
    phi0=1+phi_background_offset+t+2*x+3*y+4*z+t*y
    phih=2+2*t-x+y+3*z+x*z
    phik=-1-t+x+2*y-2*z+t*x
    phi=phi0+e*phih+f*phik
    varpi=3+t+e*(2+x)+f*(-1+z)
    Pi=Pi0+e*Ph+f*Pk
    Pip=Pi.subs(at);dPi=[Pi.diff(q).subs(at) for q in coords]
    v=S.Matrix([S.diff(phi0,q).subs(at) for q in coords]);vh=S.Matrix([S.diff(phih,q).subs(at) for q in coords]);vk=S.Matrix([S.diff(phik,q).subs(at) for q in coords])
    phi_p=phi.subs(at);varpi_p=varpi.subs(at)
    grad=S.Matrix([phi.diff(q).subs(at) for q in coords])
    B=G*Pi0.subs(at);Bh=H*Pi0.subs(at)+G*Ph.subs(at);Bk=K*Pi0.subs(at)+G*Pk.subs(at)
    Bhk=H*Pk.subs(at)+K*Ph.subs(at)
    T=S.trace(B);Th=S.trace(Bh);Tk=S.trace(Bk);Thk=S.trace(Bhk)
    tracefactor=F["Hamiltonian_trace_squared_factor"]
    V=tracefactor*T*T-S.trace(B*B)
    Vh=2*tracefactor*T*Th-2*S.trace(B*Bh)
    Vk=2*tracefactor*T*Tk-2*S.trace(B*Bk)
    Vhk=2*tracefactor*(Th*Tk+T*Thk)-2*S.trace(Bh*Bk+B*Bhk)
    Hmom_prod=F["Hamiltonian_momentum_factor"]*(Jhk*V+Jh*Vk+Jk*Vh+J*Vhk)
    bp=gp*Pip
    Hmom_exact=2*(S.trace(bp)**2/3-S.trace(bp*bp))/det_exact
    Hmom_oracle=cross(Hmom_exact)
    wh=phih.subs(at);wk=phik.subs(at);w0=phi0.subs(at)
    vp0=varpi_p.subs({e:0,f:0});vph=S.diff(varpi_p,e);vpk=S.diff(varpi_p,f)
    Hvp_prod=F["Hamiltonian_scalar_density_factor"]*(Jhk*vp0**2+2*Jh*vp0*vpk+2*Jk*vp0*vph+2*J*vph*vpk)
    Hvp_oracle=cross(-varpi_p**2/(2*det_exact))
    Hgrad_prod=F["Hamiltonian_scalar_gradient_factor"]*(scalar(v.T*Aef*v)+2*scalar(vk.T*Ah*v)+2*scalar(vh.T*Ak*v)+2*scalar(vh.T*A0*vk))
    Hgrad_oracle=cross(-scalar(grad.T*inv*grad)/2)
    eta=S.Symbol("eta");U=-S.Rational(2,27)+S.Rational(14,9)*eta**2+S.Rational(50,27)*eta**3-eta**4/6-S.Rational(4,9)*eta**5-S.Rational(2,27)*eta**6
    Hpot_prod=F["Hamiltonian_potential_factor"]*S.diff(U,eta,2).subs(eta,w0-1)*wh*wk
    Hpot_oracle=cross(-U.subs(eta,phi_p-1))
    Hfull_prod=F["Hamiltonian_curvature_factor"]*Rprod+Hmom_prod+Hvp_prod+Hgrad_prod+Hpot_prod
    Hfull_oracle=cross(R_oracle/2+Hmom_exact-varpi_p**2/(2*det_exact)-scalar(grad.T*inv*grad)/2-U.subs(eta,phi_p-1))
    checks.extend([("Hamiltonian_metric_momentum_Hessian",[Hmom_prod-Hmom_oracle]),
       ("Hamiltonian_scalar_density_Hessian",[Hvp_prod-Hvp_oracle]),
       ("Hamiltonian_scalar_gradient_Hessian",[Hgrad_prod-Hgrad_oracle]),
       ("Hamiltonian_potential_Hessian",[Hpot_prod-Hpot_oracle]),
       ("complete_nonzero_state_Hamiltonian_Hessian",[Hfull_prod-Hfull_oracle])])
    # Independent covariant density divergence: both contravariant slots and
    # the weight-one density connection are charged before cancellation.
    explicit=[];divergence=[];mixed=[];wrong_weight=[]
    for i in range(n):
        density_expression=trunc(-2*sum(gp[i,k]*dPi[j][j,k] for j in range(n) for k in range(n))-
            2*sum(Pip[j,k]*dg[j][k,i] for j in range(n) for k in range(n))+
            sum(Pip[j,k]*dg[i][j,k] for j in range(n) for k in range(n))+varpi_p*grad[i])
        density_div=[];weight0_div=[]
        for kidx in range(n):
            raw=sum(dPi[j][j,kidx]+sum(gamma[j][j][l]*Pip[l,kidx]+gamma[kidx][j][l]*Pip[j,l]-gamma[l][l][j]*Pip[j,kidx]
                                      for l in range(n)) for j in range(n))
            density_div.append(trunc(raw))
            weight0_div.append(trunc(raw+sum(gamma[l][l][j]*Pip[j,kidx] for j in range(n) for l in range(n))))
        divergence.append(trunc(-2*sum(gp[i,kidx]*density_div[kidx] for kidx in range(n))+varpi_p*grad[i]))
        wrong_weight.append(trunc(-2*sum(gp[i,kidx]*weight0_div[kidx] for kidx in range(n))+varpi_p*grad[i]))
        explicit.append(density_expression)
        term1=sum(H[i,kidx]*Pk.diff(coords[j])[j,kidx].subs(at)+K[i,kidx]*Ph.diff(coords[j])[j,kidx].subs(at) for j in range(n) for kidx in range(n))
        term2=sum(Ph[j,kidx].subs(at)*k.diff(coords[j])[kidx,i].subs(at)+Pk[j,kidx].subs(at)*h.diff(coords[j])[kidx,i].subs(at) for j in range(n) for kidx in range(n))
        term3=sum(Ph[j,kidx].subs(at)*k.diff(coords[i])[j,kidx].subs(at)+Pk[j,kidx].subs(at)*h.diff(coords[i])[j,kidx].subs(at) for j in range(n) for kidx in range(n))
        matter=vph*vk[i]+vpk*vh[i]
        mixed.append(F["momentum_h_dPi_factor"]*term1+F["momentum_Pi_dh_factor"]*term2+
                     F["momentum_Pi_trace_dh_factor"]*term3+F["momentum_matter_factor"]*matter)
    checks.extend([("full_density_divergence_vs_coordinate_formula",[trunc(a-b) for a,b in zip(divergence,explicit)]),
       ("complete_density_momentum_Hessian",[mixed[i]-cross(divergence[i]) for i in range(n)]),
       ("fixed_background_density_source_normalization",[mixed[i]/2-cross(divergence[i]/2) for i in range(n)])])
    # Normalize the canonical density by the fixed background rho0, including
    # its spatial derivative. Compare the complete family in the background
    # connection with the independent raw density divergence.
    tau=Pip/2;omega=varpi_p/2
    density_gradient=[sum(gam0[a][a][r] for a in range(n)) for r in range(n)]
    dtau=[dPi[r]/2-density_gradient[r]*tau for r in range(n)]
    divtau=[sum(dtau[j][j,kidx]+sum(gam0[j][j][l]*tau[l,kidx]+gam0[kidx][j][l]*tau[j,l] for l in range(n)) for j in range(n)) for kidx in range(n)]
    nabg=[S.Matrix(n,n,lambda i,j:dg[r][i,j]-sum(gam0[l][r][i]*gp[l,j]+gam0[l][r][j]*gp[i,l] for l in range(n))) for r in range(n)]
    normalized=[trunc(-2*sum(gp[i,kidx]*divtau[kidx] for kidx in range(n))-
       2*sum(tau[j,kidx]*nabg[j][kidx,i] for j in range(n) for kidx in range(n))+
       sum(tau[j,kidx]*nabg[i][j,kidx] for j in range(n) for kidx in range(n))+omega*grad[i]) for i in range(n)]
    checks.append(("full_fixed_background_tau_connection_formula",[trunc(normalized[i]-divergence[i]/2) for i in range(n)]))
    if not any(trunc(a-b)!=0 for a,b in zip(wrong_weight,divergence)):raise ValueError("density weight premise witness degenerate")
    if not any(cross(divergence[i]/S.sqrt(det_exact))-cross(divergence[i]/2)!=0 for i in range(n)):raise ValueError("fixed-density premise witness degenerate")
    if not any(v!=0 for v in Ric0):raise ValueError("curved background witness degenerate")
    metadata={"dimension":4,"background_metric_at_germ":[4,1,1,1],"fixed_background_density_at_germ":"2",
      "germ":"x0=x1=x2=x3=0","nonzero_background_connection_and_curvature":True,
      "density_weight_zero_rejected":True,"perturbed_volume_source_normalization_rejected":True,
      "nonzero_state_connection_terms_tested_in_fixed_chart":True,
      "scalar_background_offset":"generic symbolic parameter; all polynomial potential coefficients are retained",
      "matter_sectors":["metric momentum density","scalar gradient","scalar momentum density","actual polynomial potential","varpi dphi"],
      "scope":"Finite exact rational/polynomial two-direction Hessian witnesses at a nonzero curved 4D canonical state. These checks do not prove general analytic estimates, a positive product ball, a K constant, or a contraction radius."}
    return checks,metadata


def evaluate(claims,check):
    if S.__version__!="1.14.0" or mpmath.__version__!="1.3.0":raise ValueError("requires pinned exact engine")
    if claims.get("schema")!="hdblast-nonlinear-ball-constraint-algebra-claims-v1" or set(claims)!={"schema",*FACTORS}:raise ValueError("invalid claims schema/keys")
    checks,metadata=build_checks();sub={F[k]:S.Rational(claims[k]) for k in FACTORS}
    records=[];failures=[]
    for name,residuals in checks:
        bad=[clean(v.subs(sub)) for v in residuals];bad=[v for v in bad if v!=0]
        check("constraint_algebra/"+name,not bad)
        records.append({"name":name,"components":len(residuals),"status":"fail" if bad else "pass"})
        if bad:failures.append({"name":name,"nonzero_exact_residual":str(bad[0])})
    if failures:raise ValueError("constraint oracle rejects claims: "+json.dumps(failures,sort_keys=True))
    return {"schema":"hdblast-nonlinear-ball-constraint-algebra-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),
       "checks":records,"metadata":metadata,"engine":"sympy==1.14.0; mpmath==1.3.0"}


def controls(claims):
    records=[]
    for name in FACTORS:
        bad=copy.deepcopy(claims);bad[name]=str(S.Rational(bad[name])+1);encoded=canonical(bad)
        try:evaluate(json.loads(encoded),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name+"_plus_one","self_rehashed_sha256":hashlib.sha256(encoded).hexdigest(),
            "rejected_by":"fresh exact independent canonical-family oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving canonical factor fault: "+name)
    return {"schema":"hdblast-nonlinear-ball-constraint-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("constraint_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text());report=evaluate(claims,lambda name,ok:None);encoded=canonical(report)
        if args.report:args.report.write_bytes(encoded)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(encoded);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":raise SystemExit(main())
