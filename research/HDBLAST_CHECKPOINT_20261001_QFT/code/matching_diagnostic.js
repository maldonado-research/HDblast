"use strict";
function matchingDiagnostic(){
const rb=7.838285538073243,phi=1.0361021913225144,x=phi-.5,G=100,M=G*x,cut=53.61022411468734,b=cut**-3,delta=.1,c=.5975949350280132,d=-3.106933495673783,W=1-phi+phi**3/3;
const sigma=2*W+delta*(1+c*phi+d*phi*phi/2),sigmap=2*(phi*phi-1)+delta*(c+d*phi),sh=rb*sigma,sph=rb*sigmap;
const rows=[M,1].map(mu=>{const L=Math.log(M*M/(mu*mu)),V=M**4/(64*Math.PI**2)*(L-1.5),J=G**4*x**3/(16*Math.PI**2)*(L-1),ds=b*V,dsp=b*J;
return {mu_hat:mu,V1_hat:V,J1_hat:J,Delta_sigma_hat:ds,Delta_sigma_phi_hat:dsp,fixed_phi_vacuum_decomposition_shift:sh*ds/18-sph*dsp/24+ds*ds/36-dsp*dsp/48};});
return {status:"exploratory scheme-dependent matching illustration",endpoint:"older A1 endpoint H0tau=7.353355110666192; not the new mode endpoint",rb,phi,G,M_hat:M,characteristic_cutoff_hat:cut,b_mass_screen:b,sigma_hat:sh,sigma_phi_hat:sph,rows,limitations:["unretuned flat-space MS-bar local terms, not physical scale-dependent predictions","not produced radiation; no absolute renormalized stress inferred","fixed-phi decomposition excludes induced curvature operators and bulk/scalar readjustment","characteristic mass screen is not full-band UV certification","large-mass local approximation is inapplicable at the mass-zero crossing"]};
}
console.log(JSON.stringify(matchingDiagnostic(),null,2));
