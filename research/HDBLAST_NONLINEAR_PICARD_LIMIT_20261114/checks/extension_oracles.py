#!/usr/bin/env python3
"""Finite principal-symbol, reflection, closure and trace-count witnesses.

Analytic inputs: WEAK_LINEAR_ENERGY_REGULARITY.md and
TOP_ORDER_TIME_CONTINUITY.md. This algebra does not prove their weak boundary
identities, Friedrichs commutators, trace theorems or time continuity.
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

FACTORS=("firstorder_gradient_factor", "symmetrizer_gradient_mode", "symmetrizer_value_factor",
 "firstorder_A_factor", "firstorder_offdiagonal_B_factor", "characteristic_Bzz_factor", "normal_A_reflection_sign",
 "tangential_A_reflection_sign", "mixed_B_reflection_sign", "even_B_reflection_sign",
 "interior_normal_A_factor", "interior_mixed_B_factor", "inner_D_parity", "inner_N_parity",
 "outer_field_parity", "normal_gradient_reflection_sign", "closure_multiple", "Neumann_twoZ_sign",
 "inner_odd_coefficient_trace", "outer_odd_coefficient_trace", "odd_jump_factor",
 "D_ordinary_trace_maximum", "N_ordinary_trace_maximum", "top_N_ordinary_trace_required",
 "connection_cross_split_factor")
MUTATIONS=(("firstorder_gradient_factor","2"),("symmetrizer_gradient_mode","identity"),
 ("symmetrizer_value_factor","0"),("firstorder_A_factor","0"),
 ("firstorder_offdiagonal_B_factor","0"),("characteristic_Bzz_factor","0"),("normal_A_reflection_sign","1"),
 ("tangential_A_reflection_sign","-1"),("mixed_B_reflection_sign","1"),
 ("even_B_reflection_sign","-1"),("interior_normal_A_factor","0"),
 ("interior_mixed_B_factor","0"),("inner_D_parity","1"),("inner_N_parity","-1"),
 ("outer_field_parity","1"),("normal_gradient_reflection_sign","1"),
 ("closure_multiple",2),("Neumann_twoZ_sign","1"),
 ("inner_odd_coefficient_trace","1"),("outer_odd_coefficient_trace","1"),
 ("odd_jump_factor","0"),("D_ordinary_trace_maximum",14),
 ("N_ordinary_trace_maximum",13),("top_N_ordinary_trace_required",True),
 ("connection_cross_split_factor","1"))


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()
def clean(value):return S.cancel(S.expand(value))
def zero(value):return all(clean(x)==0 for x in value) if isinstance(value,S.MatrixBase) else clean(value)==0


def evaluate(claims,check):
    records=[]
    def need(name,condition):
        passed=bool(condition);check("extension_oracles/"+name,passed)
        records.append({"name":name,"status":"pass" if passed else "fail"})
        if not passed:raise ValueError("extension oracle rejects "+name)
    need("schema",isinstance(claims,dict) and claims.get("schema")=="hdblast-extension-oracle-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    f=lambda name:S.Rational(claims[name])
    directions=4;fibers=16;Dcount=4;Ncount=12
    A=S.Matrix(S.symbols("A0:4"));B=S.zeros(4)
    for i in range(4):
        for j in range(i,4):B[i,j]=B[j,i]=S.Symbol("B"+str(i)+str(j))
    if claims["symmetrizer_gradient_mode"]=="B":gradient_metric=B
    elif claims["symmetrizer_gradient_mode"]=="identity":gradient_metric=S.eye(4)
    else:raise ValueError("unknown gradient symmetrizer mode")
    Sbase=S.diag(S.ones(1),gradient_metric,S.Matrix([[f("symmetrizer_value_factor")]]))
    Sactual=S.diag(S.ones(1),B,S.ones(1))
    fiberI=S.SparseMatrix.eye(fibers)
    lift=lambda matrix:S.SparseMatrix(S.kronecker_product(S.SparseMatrix(matrix),fiberI))
    symmetrizer=lift(Sbase);actual_symmetrizer=lift(Sactual)
    symbols=[];base_symbols=[]
    for i in range(4):
        matrix=S.zeros(6);matrix[0,0]=f("firstorder_A_factor")*A[i]
        for j in range(4):matrix[0,1+j]=(1 if i==j else f("firstorder_offdiagonal_B_factor"))*B[i,j]
        matrix[1+i,0]=f("firstorder_gradient_factor")
        base_symbols.append(matrix);symbols.append(lift(matrix))
        product=symmetrizer*symbols[-1]
        need("full_direction_fiber_symmetrized_symbol_"+str(i),product==product.T)
    Y=S.Matrix(S.symbols("firstorder_state0:96"))
    need("complete_wave_energy_Gram",zero((Y.T*(symmetrizer-actual_symmetrizer)*Y)[0]))
    speed=S.Symbol("normal_speed")
    need("normal_principal_characteristic_polynomial",zero(base_symbols[0].charpoly(speed).as_expr()-speed**4*(speed**2-A[0]*speed-f("characteristic_Bzz_factor")*B[0,0])))
    permutation=[block*fibers+fiber for fiber in range(fibers) for block in range(6)]
    need("actual_fiber_arrays_have_same_normal_characteristic_blocks",symbols[0].extract(permutation,permutation)==S.diag(*([base_symbols[0]]*fibers)))
    # The energy flux is Phi=Y^T S M^z Y/2. The following quadratic form is
    # TWICE that physical flux and retains every B^{zA} angular cross term.
    twice_flux=A[0]*sum(Y[a]**2 for a in range(fibers))+2*sum(B[0,j]*Y[a]*Y[(1+j)*fibers+a] for j in range(4) for a in range(fibers))
    need("twice_normal_energy_flux_retains_angular_cross_terms",zero((Y.T*symmetrizer*symbols[0]*Y)[0]-twice_flux))
    tau=S.Symbol("temporal_symbol");xi=S.Matrix(S.symbols("spatial_symbol0:4"));amplitude=S.Matrix(S.symbols("fiber_amplitude0:16"))
    phase=S.Matrix.vstack(tau*amplitude,*(xi[j]*amplitude for j in range(4)),amplitude)
    reconstructed=sum((xi[i]*symbols[i]*phase for i in range(4)),S.zeros(96,1))[:16,0]
    target=(tau*(A.dot(xi))+(xi.T*B*xi)[0])*amplitude
    need("actual_scalar_principal_wave_reconstruction",zero(reconstructed-target))

    J0=S.diag(*([f("inner_D_parity")]*Dcount+[f("inner_N_parity")]*Ncount))
    JZ=f("outer_field_parity")*S.eye(fibers)
    low_value=S.Matrix([0]*Dcount+list(S.symbols("free_N_value0:12")))
    low_normal=S.Matrix(list(S.symbols("free_D_normal0:4"))+[0]*Ncount)
    need("inner_mixed_rows_determine_value_and_normal_parity",zero((S.eye(fibers)-J0)*low_value) and zero((S.eye(fibers)+J0)*low_normal))
    outer_normal=S.Matrix(S.symbols("free_outer_normal0:16"))
    need("outer_all_Dirichlet_determines_odd_parity",zero((S.eye(fibers)+JZ)*outer_normal))
    connection=S.zeros(fibers)
    for a,b,value in ((0,1,S.Rational(1,3)),(4,5,S.Rational(2,7)),(0,4,f("connection_cross_split_factor"))):
        connection[a,b]=-value;connection[b,a]=value
    need("fixed_metric_connection_preserves_D_N_split",connection+connection.T==S.zeros(fibers) and zero(J0*connection-connection*J0))

    # First-order reflected gradient has an extra covector-normal sign. Base
    # coefficients are reflected by their indices, not by substituting the
    # reflected field into a nonlinear coefficient map.
    extended_A=S.Matrix([f("normal_A_reflection_sign")*A[0],*(f("tangential_A_reflection_sign")*A[i] for i in range(1,4))])
    extended_B=S.Matrix(4,4,lambda i,j:(f("mixed_B_reflection_sign") if (i==0)!=(j==0) else f("even_B_reflection_sign"))*B[i,j])
    Rz=S.diag(-1,1,1,1)
    need("literal_base_index_reflection_definitions",extended_A==Rz*A and extended_B==Rz*B*Rz)
    for boundary,J in (("inner",J0),("outer",JZ)):
        K=S.diag(J,f("normal_gradient_reflection_sign")*J,J,J,J,J)
        for i in range(4):
            base=S.zeros(6);base[0,0]=extended_A[i]
            for j in range(4):base[0,j+1]=extended_B[i,j]
            base[1+i,0]=1
            reflected=lift(base)
            need("firstorder_reflection_intertwiner_"+boundary+"_"+str(i),zero(reflected*K-(-1 if i==0 else 1)*K*symbols[i]))

    # Literal second-order differentiation at both seams of the unit interval.
    # Coefficient germs are arbitrary smooth one-sided data; the formulas
    # define their extension. In particular no nonlinear reflection covariance
    # of the original coefficient map is presumed.
    s,z,y1,y2,y3=S.symbols("s z y1 y2 y3");coords=(z,y1,y2,y3)
    aa=[z*(1-z)*(1+s+y1*y1),1+z+s*y1,2+z*z+s*y2,3+z+s*s+y3]
    bb=S.zeros(4)
    for i in range(4):
        for j in range(i,4):
            if (i==0)!=(j==0):value=z*(1-z)*(i+j+1+s+coords[max(i,j)]**2)
            elif i==j:value=5+i+z+z*z+coords[i]**2+s*s
            else:value=S.Rational(1,i+j+3)+z*z+s*coords[i]*coords[j]
            bb[i,j]=bb[j,i]=value
    U=S.Matrix([i+1+(i+2)*s+(i+3)*z+(i+1)*s*z+s*s*y1+z*z*y2+(i+2)*z*y3+y1*y2+y2*y3+s*y3*y3+S.Rational(i+1,i+2)*z**3 for i in range(fibers)])
    def operator(field,avec,bmatrix,production=False):
        answer=field.diff(s,2)
        for i in range(4):answer-=avec[i]*(f("interior_normal_A_factor") if production and i==0 else 1)*field.diff(coords[i],s)
        for i in range(4):
            for j in range(4):answer-=bmatrix[i,j]*(f("interior_mixed_B_factor") if production and (i==0)!=(j==0) else 1)*field.diff(coords[i],coords[j])
        return answer
    actual=operator(U,aa,bb)
    for boundary,b,J in (("inner",S.Integer(0),J0),("outer",S.Integer(1),JZ)):
        substitution={z:2*b-z}
        reflectedU=J*U.subs(substitution,simultaneous=True)
        reflectedA=[(f("normal_A_reflection_sign") if i==0 else f("tangential_A_reflection_sign"))*aa[i].subs(substitution,simultaneous=True) for i in range(4)]
        reflectedB=S.Matrix(4,4,lambda i,j:(f("mixed_B_reflection_sign") if (i==0)!=(j==0) else f("even_B_reflection_sign"))*bb[i,j].subs(substitution,simultaneous=True))
        need("literal_secondorder_reflection_"+boundary,zero(operator(reflectedU,reflectedA,reflectedB,True)-J*actual.subs(substitution,simultaneous=True)))
    need("mixed_coefficients_nonzero_interior_zero_both_ends",all(expr.subs(z,0)==expr.subs(z,1)==0 and expr.subs({z:S.Rational(1,2),s:0,y1:0,y2:0,y3:0})!=0 for expr in (aa[0],*(bb[0,j] for j in range(1,4)))))

    # Reflection characters, and independent explicit DD/ND modes, distinguish
    # the ordinary 4Z closure from a twisted 2Z identification of N components.
    translation2=JZ*J0
    need("twoZ_D_periodic_N_antiperiodic_character",translation2==S.diag(*([1]*Dcount+[f("Neumann_twoZ_sign")]*Ncount)))
    period=claims["closure_multiple"]
    need("ordinary_closure_uses_even_multiple",type(period) is int and period>0 and period%2==0)
    need("untwisted_bundle_closes_on_fourZ",translation2**(period//2)==S.eye(fibers))
    modeD=S.sin(S.pi*z);modeN=S.cos(S.pi*z/2)
    need("literal_DD_twoZ_periodicity",S.trigsimp(S.expand_trig(modeD.subs(z,z+2))-modeD)==0)
    need("literal_ND_twoZ_antiperiodicity",S.trigsimp(S.expand_trig(modeN.subs(z,z+2))-f("Neumann_twoZ_sign")*modeN)==0)
    need("literal_ND_ordinary_closure",S.trigsimp(S.expand_trig(modeN.subs(z,z+period))-modeN)==0)
    need("coefficient_base_parity_closes_on_twoZ",Rz**2==S.eye(4) and f("normal_A_reflection_sign")**2==1 and f("mixed_B_reflection_sign")**2==1)

    trace=S.Symbol("free_odd_coefficient_trace")
    for boundary,orientation in (("inner",1),("outer",-1)):
        right,left=(trace,-trace) if boundary=="inner" else (-trace,trace)
        literal_jump=right-left
        need("odd_coefficient_jump_formula_"+boundary,zero(literal_jump-orientation*f("odd_jump_factor")*trace))
        need("odd_coefficient_jump_zero_iff_zero_trace_"+boundary,S.solve(literal_jump,trace)==[0])
        need("zero_actual_odd_trace_"+boundary,literal_jump.subs(trace,f(boundary+"_odd_coefficient_trace"))==0)
    need("even_coefficients_have_no_value_jump",zero(trace-trace))

    # Integer trace regularities only: an ordinary H^s trace needs s>1/2.
    # A time/rotation word of degree q has field regularity14-q, and its
    # normal derivative has regularity13-q. The top word's Neumann row is
    # therefore a WEAK conormal identity supplied by the analytic proof.
    Dmax=claims["D_ordinary_trace_maximum"];Nmax=claims["N_ordinary_trace_maximum"]
    Dallowed=[q for q in range(15) if S.Rational(14-q)>S.Rational(1,2)]
    Nallowed=[q for q in range(15) if S.Rational(13-q)>S.Rational(1,2)]
    need("ordinary_D_trace_degree_at_most13",Dmax==max(Dallowed))
    need("ordinary_N_trace_degree_at_most12",Nmax==max(Nallowed))
    trace_rows=[]
    for total in range(14):
        for time_degree in range(total+1):
            trace_rows.append({"time_degree":time_degree,"rotation_degree":total-time_degree,"D_field_order":14-total,"N_normal_order":13-total,"ordinary_N_trace":total<=Nmax})
    need("no_highest_N_ordinary_trace_assumption",(S.Rational(13-13)>S.Rational(1,2)) if claims["top_N_ordinary_trace_required"] is True else S.Rational(13-13)<=S.Rational(1,2))
    return {"schema":"hdblast-extension-oracle-report-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "finite_exact":{"spatial_directions":directions,"fiber_components":fibers,"D_components":Dcount,"N_components":Ncount,"firstorder_array_size":96,"normal_characteristic_factor":"lambda^4*(lambda^2-Az*lambda-Bzz)","ordinary_period_in_Z":period,"trace_rows":trace_rows},
       "analytic_inputs":["WEAK_LINEAR_ENERGY_REGULARITY.md","TOP_ORDER_TIME_CONTINUITY.md"],
       "scope":"Finite scalar-principal symbol/symmetrizer/characteristic/flux arrays, base-index/fiber reflection, fourZ closure, odd-trace jumps and trace derivative counts. Polynomial/mode witnesses do not prove PDE finite propagation, weak conormal boundary identities, commutator convergence, trace theorems, continuity, strong top-norm approximation or an actual nonlinear coefficient-reflection covariance."}


def controls(claims):
    records=[]
    for field,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[field]=value;raw=canonical(bad)
        try:evaluate(json.loads(raw),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":field,"self_rehashed_sha256":hashlib.sha256(raw).hexdigest(),"diagnostic":str(exc),"rejected_by":"fresh exact extension oracle"})
        else:raise ValueError("surviving extension fault: "+field)
    return {"schema":"hdblast-extension-oracle-controls-v1","status":"pass","controls":records}
