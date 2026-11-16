#!/usr/bin/env python3
"""Finite exact witnesses for subsidiary-closure interfaces.

These algebraic examples test the five-dimensional trace reverse, the current
half-action densities, full covector Robin conversion, coupled wave changes
of frame and boundary flux. A separate small auxiliary chart/connection model
tests reverse source conversion. They do not prove the actual continuum
subsidiary equation, weak boundary passage, uniqueness or physical closure.
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

FACTORS=("symmetric_gradient_factor","Einstein_trace_factor","inverse_trace_factor",
 "initial_normal_projection_factor","initial_mixed_projection_factor",
 "canonical_K_factor","canonical_K_trace_factor","density_determinant_exponent",
 "canonical_scalar_velocity_factor","canonical_H_curvature_factor","canonical_H_metric_factor",
 "canonical_H_scalar_factor","momentum_divergence_factor","momentum_matter_factor",
 "momentum_projection_sign","coordinate_Robin_K_factor","y_to_z_sign",
 "conjugator_first_transport_factor","conjugator_second_factor","conjugator_time_factor",
 "conjugator_angular_factor","conjugator_lower_derivative_factor","coupled_lower_factor",
 "lower_matrix_conjugation_mode","principal_covector_mode","normal_time_cross",
 "normal_angular_cross","split_metric_cross","Dirichlet_time_trace","Neumann_normal_trace",
 "firstjet_connection_factor","Q_first_connection_factor","Q_connection_derivative_factor",
 "Q_connection_square_factor","Q_base_connection_factor","chart_base_second_factor",
 "chart_base_state_factor","chart_state_second_factor","chart_Jacobian_mode",
 "normalized_mixed_factor","normalized_spatial_factor","normalized_source_factor",
 "unnormalize_factor","reference_K_factor","reference_dK_factor",
 "metric_wave_row_factor","scalar_wave_row_factor","literal_row_order")
MUTATIONS=(("symmetric_gradient_factor","1"),("Einstein_trace_factor","1/3"),
 ("inverse_trace_factor","1/2"),("initial_normal_projection_factor","1"),
 ("initial_mixed_projection_factor","1"),("canonical_K_factor","1"),
 ("canonical_K_trace_factor","1/4"),("density_determinant_exponent","-1"),
 ("canonical_scalar_velocity_factor","2"),("canonical_H_curvature_factor","1"),
 ("canonical_H_metric_factor","1"),("canonical_H_scalar_factor","-1"),
 ("momentum_divergence_factor","-1"),("momentum_matter_factor","0"),
 ("momentum_projection_sign","1"),("coordinate_Robin_K_factor","1"),
 ("y_to_z_sign","1"),("conjugator_first_transport_factor","1"),
 ("conjugator_second_factor","0"),("conjugator_time_factor","0"),
 ("conjugator_angular_factor","0"),("conjugator_lower_derivative_factor","0"),
 ("coupled_lower_factor","0"),("lower_matrix_conjugation_mode","unconjugated"),
 ("principal_covector_mode","omit_time"),("normal_time_cross","1/7"),
 ("normal_angular_cross","1/7"),("split_metric_cross","1/7"),
 ("Dirichlet_time_trace",False),("Neumann_normal_trace",False),
 ("firstjet_connection_factor","1"),("Q_first_connection_factor","1"),
 ("Q_connection_derivative_factor","0"),("Q_connection_square_factor","1"),
 ("Q_base_connection_factor","0"),("chart_base_second_factor","0"),
 ("chart_base_state_factor","1"),("chart_state_second_factor","0"),
 ("chart_Jacobian_mode","identity"),("normalized_mixed_factor","1"),
 ("normalized_spatial_factor","0"),("normalized_source_factor","1"),
 ("unnormalize_factor","1"),("reference_K_factor","0"),("reference_dK_factor","0"),
 ("metric_wave_row_factor","-1"),("scalar_wave_row_factor","2"),
 ("literal_row_order","before_J"))


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()
def clean(value):return S.cancel(S.expand(value))
def zero(value):return all(clean(x)==0 for x in value) if isinstance(value,S.MatrixBase) else clean(value)==0
def symmatrix(prefix,n):
    out=S.zeros(n)
    for i in range(n):
        for j in range(i,n):out[i,j]=out[j,i]=S.Symbol(prefix+str(i)+str(j))
    return out


def evaluate(claims,check):
    records=[]
    def need(name,condition):
        passed=bool(condition);check("subsidiary_oracles/"+name,passed)
        records.append({"name":name,"status":"pass" if passed else "fail"})
        if not passed:raise ValueError("subsidiary oracle rejects "+name)
    need("schema",isinstance(claims,dict) and claims.get("schema")=="hdblast-subsidiary-oracle-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    f=lambda name:S.Rational(claims[name])

    # Nonzero lapse/shift, fully non-diagonal spatial metric. The selected
    # initial family has zero shift; this is an auxiliary chart identity.
    L=S.Matrix([[2,0,0,0],[S.Rational(1,3),3,0,0],[S.Rational(1,5),S.Rational(2,7),1,0],
                [S.Rational(2,9),S.Rational(1,11),S.Rational(1,13),1]])
    gamma=L*L.T;A=gamma.inv();sg=S.sqrt(gamma.det());rho0=S.Rational(5,7)
    alpha=S.Rational(7,5);beta=S.Matrix([S.Rational(1,7),-S.Rational(2,9),S.Rational(1,11),S.Rational(3,13)])
    g=S.zeros(5);g[0,0]=-alpha**2+(beta.T*gamma*beta)[0]
    for i in range(4):
        g[0,i+1]=g[i+1,0]=(gamma*beta)[i]
        for j in range(4):g[i+1,j+1]=gamma[i,j]
    inverse=g.inv();u=S.Matrix([1/alpha,*(-beta/alpha)])
    spatial=[S.eye(5)[:,i+1] for i in range(4)];X=S.Matrix(S.symbols("initial_X0:5"))
    normal_covector=g*u
    derivative=-normal_covector*X.T
    Eoracle=(derivative+derivative.T)/2
    E=f("symmetric_gradient_factor")*(derivative+derivative.T)
    trE=S.trace(inverse*Eoracle);Horacle=Eoracle-g*trE/2
    H=E-f("Einstein_trace_factor")*g*S.trace(inverse*E)
    need("ADM_unit_normal_and_tangential_split",zero((u.T*g*u)[0]+1) and all(zero((u.T*g*e)[0]) for e in spatial))
    need("five_component_symmetric_covector_gradient",zero(E-Eoracle))
    need("five_dimensional_Einstein_trace_reverse",zero(H-Horacle))
    need("inverse_five_dimensional_trace_reverse",zero(H-f("inverse_trace_factor")*g*S.trace(inverse*H)-Eoracle))
    Xu=(X.T*u)[0]
    need("initial_normal_constraint_projection",zero((u.T*Horacle*u)[0]-f("initial_normal_projection_factor")*Xu))
    need("initial_four_momentum_constraint_projections",all(zero((u.T*Horacle*spatial[i])[0]-f("initial_mixed_projection_factor")*X[i+1]) for i in range(4)))
    projections=S.Matrix([(u.T*Horacle*u)[0],*((u.T*Horacle*e)[0] for e in spatial)])
    need("all_five_normal_derivatives_determined_by_constraints",projections.jacobian(X).det()!=0)

    # The physical forward momentum map is the independent oracle for the
    # canonical inverse. No doubled-action normalization is used.
    tau=symmatrix("tau",4);T=S.trace(gamma*tau);omega=S.Symbol("omega")
    determinant_ratio=gamma.det()/rho0**2
    lam=determinant_ratio**f("density_determinant_exponent")
    K=f("canonical_K_factor")*lam*(gamma*tau*gamma-f("canonical_K_trace_factor")*T*gamma)
    traceK=S.trace(A*K)
    pi_forward=sg*(A*K*A-traceK*A)/2
    need("half_action_metric_density_inverse",zero(pi_forward-rho0*tau))
    Pi=f("canonical_scalar_velocity_factor")*lam*omega
    need("half_action_scalar_density_inverse",zero(sg*Pi-rho0*omega))
    R,U=S.symbols("spatial_R potential_U");p=S.Matrix(S.symbols("spatial_phi0:4"))
    physical_H=(R+traceK**2-S.trace(A*K*A*K))/2-Pi**2/2-(p.T*A*p)[0]/2-U
    canonical_H=(f("canonical_H_curvature_factor")*R+f("canonical_H_metric_factor")*lam**2*(T*T/3-S.trace(gamma*tau*gamma*tau))+
       f("canonical_H_scalar_factor")*lam**2*omega**2-(p.T*A*p)[0]/2-U)
    need("full_half_action_Hamiltonian_projection",zero(canonical_H-physical_H))
    xx=S.symbols("initial_space0:4");Kjets=[symmatrix("dK"+str(i)+"_",4) for i in range(4)]
    localK=K+sum((xx[j]*Kjets[j] for j in range(4)),S.zeros(4))
    localtrace=S.trace(A*localK);localpi=sg*(A*localK*A-localtrace*A)/2
    momentum=S.Matrix([f("momentum_divergence_factor")*sum(gamma[i,k]*S.diff(localpi[j,k],xx[j]) for j in range(4) for k in range(4))+
       f("momentum_matter_factor")*sg*Pi*p[i] for i in range(4)])
    codazzi=S.Matrix([sum(A[j,k]*S.diff(localK[k,i],xx[j]) for j in range(4) for k in range(4))-S.diff(localtrace,xx[i])-Pi*p[i] for i in range(4)])
    need("full_metric_scalar_momentum_density_projection",zero(momentum/sg-f("momentum_projection_sign")*codazzi))

    # Independent linearized Ricci/contracted-Christoffel principal symbol.
    amplitude=symmatrix("metric_amplitude",5)
    frequency=S.Matrix([1,-S.Rational(2,3),S.Rational(3,5),S.Rational(1,7),-S.Rational(1,11)])
    connection=[[[sum(inverse[a,b]*(frequency[m]*amplitude[b,n]+frequency[n]*amplitude[b,m]-frequency[b]*amplitude[m,n]) for b in range(5))/2
                   for n in range(5)] for m in range(5)] for a in range(5)]
    Ricci=S.Matrix(5,5,lambda m,n:sum(frequency[a]*connection[a][m][n]-frequency[n]*connection[a][m][a] for a in range(5)))
    gauge_up=S.Matrix([sum(inverse[m,n]*connection[a][m][n] for m in range(5) for n in range(5)) for a in range(5)])
    gauge_down=g*gauge_up
    gauge_derivative=(frequency*gauge_down.T+gauge_down*frequency.T)/2
    box=(frequency.T*inverse*frequency)[0]
    need("metric_wave_to_literal_Einstein_row",zero(Ricci-gauge_derivative-f("metric_wave_row_factor")*box*amplitude))
    scalar_amplitude=S.Symbol("scalar_amplitude")
    need("scalar_wave_to_literal_scalar_row",zero(box*scalar_amplitude-f("scalar_wave_row_factor")*box*scalar_amplitude))

    # The covector connection contributes one K, and E(n,a)=0 contributes
    # another K. In coordinates the resulting Robin coefficient is full2K.
    shell_h=S.Matrix([[-2,S.Rational(1,5),0,0],[S.Rational(1,5),3,S.Rational(1,7),0],
                      [0,S.Rational(1,7),4,S.Rational(1,9)],[0,0,S.Rational(1,9),5]])
    shape_low=symmatrix("shell_K",4);shape_cov=shape_low*shell_h.inv()
    Ctan=S.Matrix(S.symbols("shell_C0:4"));Cyjet=S.Matrix(S.symbols("shell_dy_C0:4"))
    radial_root=S.Rational(6,5);rhob=S.Rational(129,1)
    covariant_normal=Cyjet/radial_root-shape_cov*Ctan
    geometric_robin=covariant_normal-shape_cov*Ctan
    coordinate_robin=Cyjet/radial_root-f("coordinate_Robin_K_factor")*shape_cov*Ctan
    need("covector_Robin_is_full_two_K",zero(coordinate_robin-geometric_robin))
    z=S.Symbol("z");y=S.Symbol("y");map_y=S.Rational(4,9)-rhob*z
    Cchart=S.Matrix([y+y*y,1+2*y+y**3,y*y+y**3,2-y+y**4])
    need("literal_inward_y_to_z_chain_sign",zero(Cchart.subs(y,map_y).diff(z)-f("y_to_z_sign")*rhob*Cchart.diff(y).subs(y,map_y)))
    five_covector=S.Matrix(S.symbols("coordinate_C0:5"))
    jacobian=S.diag(1,-rhob,1,1,1)
    covector_pullback=jacobian.T*five_covector
    need("normal_covector_component_coordinate_conversion",zero(covector_pullback[1]-f("y_to_z_sign")*rhob*five_covector[1]))
    sigma=S.Symbol("sigma")
    need("umbilic_full_Robin_is_Nsigma_over_three",zero(f("coordinate_Robin_K_factor")*radial_root*(sigma*S.eye(4)/6)-radial_root*sigma*S.eye(4)/3))

    # Block-preserving matrix conjugator with genuine time/angular dependence;
    # fully coupled five-component lower matrices are retained throughout.
    t,z,x1,x2,x3=S.symbols("t z x1 x2 x3");coords=(t,z,x1,x2,x3)
    Rmat=S.zeros(5)
    # Bundle order is (normal,time,three angular); base order is
    # (t,z,x1,x2,x3). This is an invertible POLYNOMIAL conjugator witness,
    # not a calculation of the exponential's Frechet derivative formula.
    Rmat[1,2]=1+t+x1+t*x1;Rmat[2,3]=2+t*t+x2;Rmat[3,4]=3+x3+x1*x2
    transform=S.eye(5)-z*Rmat
    transform_inv=S.eye(5)+z*Rmat+z*z*Rmat**2+z**3*Rmat**3
    need("block_matrix_conjugator_inverse",zero(transform*transform_inv-S.eye(5)))
    need("Robin_conjugator_boundary_jet",zero((transform.diff(z)+Rmat*transform).subs(z,0)))
    need("nonzero_mixed_time_angular_conjugator_jet",not zero(transform.diff(t,x1)))
    WW=S.Matrix([1+t+z*z+x1,2+t*z+x2*x2,3+t*t+x1*x3,4+z*x2+x3*x3,5+t*x3+z*x1])
    q=S.Matrix([[-2,S.Rational(1,7),S.Rational(1,11),-S.Rational(1,13),S.Rational(1,17)],
       [S.Rational(1,7),3,S.Rational(1,5),S.Rational(1,9),-S.Rational(1,15)],
       [S.Rational(1,11),S.Rational(1,5),4,S.Rational(1,8),S.Rational(1,10)],
       [-S.Rational(1,13),S.Rational(1,9),S.Rational(1,8),5,S.Rational(1,12)],
       [S.Rational(1,17),-S.Rational(1,15),S.Rational(1,10),S.Rational(1,12),6]])
    lower=[S.Matrix(5,5,lambda i,j:S.Rational((i+1)*(j+2)+mu+1,19+mu)) for mu in range(5)]
    zeroth=S.Matrix(5,5,lambda i,j:S.Rational(i+2*j+1,23))
    physical=transform*WW
    literal=(sum((q[i,j]*physical.diff(coords[i],coords[j]) for i in range(5) for j in range(5)),S.zeros(5,1))+
       sum((lower[i]*physical.diff(coords[i]) for i in range(5)),S.zeros(5,1))+zeroth*physical)
    def selected(matrix):return S.Matrix(5,5,lambda i,j:matrix[i,j]*(1 if i==j else f("coupled_lower_factor")))
    adjusted_lower=[selected(B) for B in lower];adjusted_zeroth=selected(zeroth)
    def jetweight(indices):
        value=S.Integer(1)
        if 0 in indices:value*=f("conjugator_time_factor")
        if any(index>=2 for index in indices):value*=f("conjugator_angular_factor")
        return value
    principal=sum((q[i,j]*WW.diff(coords[i],coords[j]) for i in range(5) for j in range(5)),S.zeros(5,1))
    first=sum((f("conjugator_first_transport_factor")*q[i,j]*jetweight((i,))*transform_inv*transform.diff(coords[i])*WW.diff(coords[j]) for i in range(5) for j in range(5)),S.zeros(5,1))
    second=sum((f("conjugator_second_factor")*q[i,j]*jetweight((i,j))*transform_inv*transform.diff(coords[i],coords[j])*WW for i in range(5) for j in range(5)),S.zeros(5,1))
    if claims["lower_matrix_conjugation_mode"]=="conjugated":
        first+=sum((transform_inv*adjusted_lower[i]*transform*WW.diff(coords[i]) for i in range(5)),S.zeros(5,1))
        second+=transform_inv*adjusted_zeroth*transform*WW
    elif claims["lower_matrix_conjugation_mode"]=="unconjugated":
        first+=sum((adjusted_lower[i]*WW.diff(coords[i]) for i in range(5)),S.zeros(5,1));second+=adjusted_zeroth*WW
    else:raise ValueError("unsupported lower matrix conjugation mode")
    second+=sum((f("conjugator_lower_derivative_factor")*jetweight((i,))*transform_inv*adjusted_lower[i]*transform.diff(coords[i])*WW for i in range(5)),S.zeros(5,1))
    need("coupled_five_covector_wave_matrix_conjugation",zero(transform_inv*literal-principal-first-second))
    principal_symbol=S.Symbol("five_component_principal_symbol")
    production_symbol=principal_symbol*S.eye(5)
    if claims["principal_covector_mode"]=="omit_time":production_symbol[1,1]=0
    elif claims["principal_covector_mode"]!="all_five":raise ValueError("unsupported covector component mode")
    need("principal_wave_retains_all_five_covector_components",production_symbol==principal_symbol*S.eye(5))

    # Fixed positive covector energy metric, rank1/4 orthogonal shell split.
    velocity=S.Matrix(S.symbols("boundary_velocity0:5"));normal=S.Matrix(S.symbols("boundary_normal0:5"))
    tangents=[S.Matrix(S.symbols("boundary_tangent"+str(i)+"_0:5")) for i in range(3)]
    positive_metric=S.diag(2,3,4,5,6)
    positive_metric[0,1]=positive_metric[1,0]=f("split_metric_cross")
    flux=-(f("normal_time_cross")*(velocity.T*positive_metric*velocity)[0]+
       3*(velocity.T*positive_metric*normal)[0]+f("normal_angular_cross")*sum((velocity.T*positive_metric*v)[0] for v in tangents))
    trace_rules={}
    if claims["Dirichlet_time_trace"] is True:trace_rules[velocity[0]]=0
    if claims["Neumann_normal_trace"] is True:trace_rules.update({normal[j]:0 for j in range(1,5)})
    need("block_zero_shell_flux_with_free_complementary_jets",zero(flux.subs(trace_rules,simultaneous=True)))

    # Reverse the fixed product connection in a SMALL independent chart model.
    # Actual evolution has16 fields/5 base directions; this auxiliary model
    # has3 fields/4 directions and does not evaluate its physical coefficients.
    tt,zz,xx,yy=S.symbols("aux_t aux_z aux_x aux_y");cc=(tt,zz,xx,yy)
    V=S.Matrix([S.Rational(1,9)+(tt+zz*xx+yy*yy)/20,S.Rational(1,11)+(zz+tt*xx+xx*yy)/23,
                -S.Rational(1,13)+(yy+tt*zz+xx*xx)/29])
    O=[S.zeros(3) for _ in range(4)]
    O[2]=S.Matrix([[0,-(1+yy)/7,0],[(1+yy)/7,0,0],[0,0,0]])
    O[3]=S.Matrix([[0,0,0],[0,0,-(1+xx)/11],[0,(1+xx)/11,0]])
    need("auxiliary_angular_connections_do_not_commute",not zero(O[2]*O[3]-O[3]*O[2]))
    QQ=[V.diff(cc[i])+O[i]*V for i in range(4)]
    PP=[QQ[i]+f("firstjet_connection_factor")*O[i]*V for i in range(4)]
    need("reverse_covariant_firstjet_P_equals_Q_minus_OmegaV",all(zero(PP[i]-V.diff(cc[i])) for i in range(4)))
    aq=S.Matrix([[-2,S.Rational(1,7),-S.Rational(1,11),S.Rational(1,13)],
       [S.Rational(1,7),3,S.Rational(1,5),-S.Rational(1,9)],
       [-S.Rational(1,11),S.Rational(1,5),4,S.Rational(1,8)],
       [S.Rational(1,13),-S.Rational(1,9),S.Rational(1,8),5]])
    Gamma={ (2,3,3):S.Rational(1,9),(3,2,3):S.Rational(1,10),(3,3,2):S.Rational(1,10)}
    HH=[[QQ[j].diff(cc[i])+O[i]*QQ[j]-sum((Gamma.get((k,i,j),0)*QQ[k] for k in range(4)),S.zeros(3,1)) for j in range(4)] for i in range(4)]
    correction=sum((aq[i,j]*(f("Q_first_connection_factor")*O[i]*QQ[j]+
       (f("Q_connection_derivative_factor")*O[j].diff(cc[i])+f("Q_connection_square_factor")*O[i]*O[j])*V+
       f("Q_base_connection_factor")*sum((Gamma.get((k,i,j),0)*QQ[k] for k in range(4)),S.zeros(3,1))) for i in range(4) for j in range(4)),S.zeros(3,1))
    second_coordinate=sum((aq[i,j]*V.diff(cc[i],cc[j]) for i in range(4) for j in range(4)),S.zeros(3,1))
    second_covariant=sum((aq[i,j]*HH[i][j] for i in range(4) for j in range(4)),S.zeros(3,1))
    need("reverse_connection_square_is_minus_in_Q_variables",zero(second_covariant-second_coordinate-correction))

    vv=S.Matrix(S.symbols("chart_V0:3"))
    Psi=S.Matrix([vv[0]+vv[0]*vv[1]/5+tt*vv[2]/7+zz**2/11,
                  vv[1]+vv[1]*vv[2]/6+xx*vv[0]/13+tt*yy/17,
                  vv[2]+vv[0]**2/8+yy*vv[1]/19+xx*zz/23])
    substitution=dict(zip(vv,V));Xfield=Psi.subs(substitution,simultaneous=True)
    Jexpr=Psi.jacobian(vv).subs(substitution,simultaneous=True)
    point={tt:S.Rational(1,7),zz:S.Rational(1,5),xx:-S.Rational(1,11),yy:S.Rational(2,13)}
    at=lambda expr:expr.subs(point,simultaneous=True).applyfunc(clean) if isinstance(expr,S.MatrixBase) else clean(expr.subs(point,simultaneous=True))
    actualJ=at(Jexpr)
    if claims["chart_Jacobian_mode"]=="actual":chosenJ=actualJ
    elif claims["chart_Jacobian_mode"]=="identity":chosenJ=S.eye(3)
    else:raise ValueError("unsupported chart Jacobian mode")
    chain=S.zeros(3,1)
    for i in range(4):
        for j in range(4):
            state_hessian=S.Matrix([sum(Psi[o].diff(vv[a],vv[b])*PP[i][a]*PP[j][b] for a in range(3) for b in range(3)) for o in range(3)])
            cross=Psi.diff(cc[i]).jacobian(vv)*PP[j]
            chain+=aq[i,j]*(f("chart_base_second_factor")*Psi.diff(cc[i],cc[j])+f("chart_base_state_factor")*cross+f("chart_state_second_factor")*state_hessian).subs(substitution,simultaneous=True)
    chainpoint=at(chain)
    literal_second=at(sum((aq[i,j]*Xfield.diff(cc[i],cc[j]) for i in range(4) for j in range(4)),S.zeros(3,1)))
    need("full_base_state_chart_second_chain",zero(literal_second-chosenJ*at(second_coordinate)-chainpoint))
    produced_first=[at(Psi.diff(cc[i]).subs(substitution,simultaneous=True)+Jexpr*PP[i]) for i in range(4)]
    actual_first=[at(Xfield.diff(c)) for c in cc]
    need("full_base_state_chart_first_chain",all(zero(produced_first[i]-actual_first[i]) for i in range(4)))
    reference=S.Matrix(3,3,lambda i,j:S.Rational(i+2*j+1,17)+(i+1)*tt/19+(j+1)*zz/23+xx*yy/(i+j+3))
    Xpoint=at(Xfield)
    nonlinear=S.Matrix([Xpoint[0]*Xpoint[1],actual_first[0][0]*actual_first[2][1],Xpoint[2]**2+actual_first[1][2]])
    actual_reference=at(reference)*Xpoint+sum((at(reference.diff(cc[i]))*actual_first[i] for i in range(4)),S.zeros(3,1))
    physical_source=actual_reference+nonlinear
    produced_source=(f("reference_K_factor")*at(reference)*Xpoint+
       f("reference_dK_factor")*sum((at(reference.diff(cc[i]))*produced_first[i] for i in range(4)),S.zeros(3,1))+nonlinear)
    covariant_source=chosenJ.inv()*(produced_source-chainpoint)+at(correction)
    lapse=-aq[0,0]
    linear_normalized=at(V.diff(tt,2)-sum((f("normalized_mixed_factor")*aq[0,i]/lapse*QQ[i].diff(tt) for i in range(1,4)),S.zeros(3,1))-
       sum((f("normalized_spatial_factor")*aq[i,j]/lapse*HH[i][j] for i in range(1,4) for j in range(1,4)),S.zeros(3,1)))
    normalized_source=f("normalized_source_factor")*covariant_source/lapse
    normalized_residual=linear_normalized-normalized_source
    physical_wave_residual=literal_second-physical_source
    reconstructed=f("unnormalize_factor")*lapse*chosenJ*normalized_residual
    need("minus_a_J_reconstructs_full_reference_wave_residual",zero(reconstructed-physical_wave_residual))
    row=S.diag(f("metric_wave_row_factor"),f("metric_wave_row_factor"),f("scalar_wave_row_factor"))
    literal_rows=S.diag(-S.Rational(1,2),-S.Rational(1,2),1)*physical_wave_residual
    if claims["literal_row_order"]=="after_J":produced_rows=row*reconstructed
    elif claims["literal_row_order"]=="before_J":produced_rows=f("unnormalize_factor")*lapse*chosenJ*row*normalized_residual
    else:raise ValueError("unsupported literal row order")
    need("metric_half_and_scalar_full_rows_follow_Jacobian",zero(produced_rows-literal_rows))
    return {"schema":"hdblast-subsidiary-oracle-report-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "finite_exact":{"physical_dimension":5,"spatial_dimension":4,"subsidiary_covector_components":5,"shell_D_rank":1,"shell_N_rank":4,
          "wave_bundle_order":["normal","time","angular1","angular2","angular3"],"wave_base_order":["t","z","x1","x2","x3"],
          "nonzero_shift_initial_chart":True,"auxiliary_reverse_chain_fields":3,"auxiliary_reverse_chain_directions":4},
       "scope":"Finite exact trace-reverse/half-action constraint projections, full2K covector Robin and inward coordinate signs, coupled five-component polynomial conjugation and block flux, and a small auxiliary noncommuting connection/chart/reference reverse chain. The polynomial witness does not check exponential Frechet formulas. No actual continuum subsidiary PDE or weak boundary closure, coefficient ceiling, uniqueness, amplitude change or physical confirmation is proved by these witnesses."}


def controls(claims):
    records=[]
    for field,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[field]=value;raw=canonical(bad);failures=[]
        def check(name,ok):
            if not bool(ok):failures.append(name)
        try:evaluate(json.loads(raw),check)
        except ValueError as exc:
            if len(failures)!=1 or failures[0].endswith(("/schema","/exact_engine")):
                raise ValueError("mutation did not reach a named mathematical oracle: "+field) from exc
            records.append({"fault":field,"self_rehashed_sha256":hashlib.sha256(raw).hexdigest(),"failed_check":failures[0],"diagnostic":str(exc),"rejected_by":"fresh named exact subsidiary oracle"})
        else:raise ValueError("surviving subsidiary fault: "+field)
    return {"schema":"hdblast-subsidiary-oracle-controls-v1","status":"pass","controls":records}
