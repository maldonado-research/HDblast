"use strict";
// Copyright (c) 2026 Ricardo Maldonado. MIT License.
// Post hoc diagnostic; older registered verdicts are unchanged.
const plateaus = [{"label":"Y=0.5 dc=0.01","tag":"main_Y0.5_dc1e-2_dzf5e-4","classification":"FAIL-Weyl","plateau":{"H0tau":10.529955723450618,"window_start_H0tau":4.063044072959733,"r":0.22365100578666058,"Omega_r":0.9433102843484779,"Wa4":84.72064416279474,"Ra4":166.91040650510894,"R_over_sigma":0.0013782979510745732,"H_over_H0":0.04694018819872893,"phi_b":1.0361021908827632,"Omega_vac":-0.15428249624684884,"ln_a":3.028287481492495,"efolds_to_turnaround":null,"Omega_W":0.21097229386343785,"efolds_to_turnaround_estimate":0.5031122430936011,"I_v2a4":2616.5781857840793,"Wa4_over_I_v2a4":0.03237841109548481,"Ra4_over_Y_I_v2a4":1.0000017829170131}},{"label":"Y=0.7 dc=0.01","tag":"main_Y0.7_dc1e-2_dzf5e-4","classification":"FAIL-Weyl","plateau":{"H0tau":11.47604028476125,"window_start_H0tau":4.6498049119568154,"r":0.1325787616239562,"Omega_r":1.035868251747044,"Wa4":197.67588304515343,"Ra4":656.9786746203879,"R_over_sigma":0.0013482417585224712,"H_over_H0":0.044302570029942916,"phi_b":1.0361021889708826,"Omega_vac":-0.17320227015647424,"ln_a":3.3763480827690704,"efolds_to_turnaround":null,"Omega_W":0.1373341300221956,"efolds_to_turnaround_estimate":0.47825806601686643,"I_v2a4":7356.558542371278,"Wa4_over_I_v2a4":0.02687070073684692,"Ra4_over_Y_I_v2a4":0.9999991192387125}},{"label":"Y=1.5 dc=0.01","tag":"main_Y1.5_dc1e-2_dzf5e-4","classification":"PASS-conservative","plateau":{"H0tau":15.442739568636451,"window_start_H0tau":6.954608345312029,"r":0.03602188025961443,"Omega_r":1.2359184721023289,"Wa4":11152.105875742182,"Ra4":136439.13774229953,"R_over_sigma":0.000993809990889772,"H_over_H0":0.03481894473279056,"phi_b":1.0361021677427096,"Omega_vac":-0.2804380654097001,"ln_a":4.78659631853216,"efolds_to_turnaround":0.38344228391613555,"Omega_W":0.044520107212715716,"efolds_to_turnaround_estimate":0.3796512596028925,"I_v2a4":712971.1363872605,"Wa4_over_I_v2a4":0.015641735417581834,"Ra4_over_Y_I_v2a4":0.9999927211749594}},{"label":"Y=2 dc=0.0001","tag":"fine_Y2_dc1e-4_dzf2.5e-4","classification":"PASS-combined","plateau":{"H0tau":21.65883094452056,"window_start_H0tau":11.90017112149363,"r":0.021299443638161504,"Omega_r":1.3557093051324265,"Wa4":216123022779.3105,"Ra4":4472232962087.1,"R_over_sigma":0.0007948613704113396,"H_over_H0":0.029730324616866127,"phi_b":1.0361021963123889,"Omega_vac":-0.38458513247056386,"ln_a":9.168759653973671,"efolds_to_turnaround":null,"Omega_W":0.028875853934399214,"efolds_to_turnaround_estimate":0.3202476685827485,"I_v2a4":17527529648066.729,"Wa4_over_I_v2a4":0.01233048964222683,"Ra4_over_Y_I_v2a4":0.9999880089630393}},{"label":"Y=2 dc=0.01","tag":"fine_Y2_dc1e-2_dzf2.5e-4","classification":"PASS-combined","plateau":{"H0tau":18.20157110470227,"window_start_H0tau":8.468273748775829,"r":0.021299555782589328,"Omega_r":1.3574158351536791,"Wa4":191113.94262868533,"Ra4":3954704.448797571,"R_over_sigma":0.0007922721837079378,"H_over_H0":0.029663180305346936,"phi_b":1.0361021963130264,"Omega_vac":-0.38632816250071517,"ln_a":5.684952205852737,"efolds_to_turnaround":0.3195150947587777,"Omega_W":0.02891235430102587,"efolds_to_turnaround_estimate":0.31943169261182525,"I_v2a4":15499239.202276252,"Wa4_over_I_v2a4":0.012330537011172647,"Ra4_over_Y_I_v2a4":0.9999878795280308}},{"label":"Y=3 dc=0.01","tag":"fine_Y3_dc1e-2_dzf1.25e-4","classification":"PASS-combined","plateau":{"H0tau":23.54814560408794,"window_start_H0tau":11.34094062408419,"r":0.010115968133639859,"Omega_r":1.6585125391359155,"Wa4":83495237.85767728,"Ra4":3638297572.0889716,"R_over_sigma":0.0005538457783922844,"H_over_H0":0.022436036086624606,"phi_b":1.036102199672502,"Omega_vac":-0.6752903595135945,"ln_a":7.480548208686116,"efolds_to_turnaround":0.22710084256926777,"Omega_W":0.01677745999514105,"efolds_to_turnaround_estimate":0.2271497006152965,"I_v2a4":9506197275.631224,"Wa4_over_I_v2a4":0.008783242703337765,"Ra4_over_Y_I_v2a4":0.9999797821580676}},{"label":"Y=5 dc=0.01","tag":"fine_Y5_dc1e-2_dzf1.25e-4_cfl0.25","classification":"UNRELIABLE","plateau":{"H0tau":35.08091613831256,"window_start_H0tau":17.004741189730353,"r":0.004408284500412442,"Omega_r":3.2169010916082157,"Wa4":37448451875640.56,"Ra4":3745052010119921.5,"R_over_sigma":0.00032520566947923904,"H_over_H0":0.012343731506890523,"phi_b":1.0361021901977028,"Omega_vac":-2.2310803346581762,"ln_a":11.074762870462434,"efolds_to_turnaround":0.09258636967137157,"Omega_W":0.014181015221496361,"efolds_to_turnaround_estimate":0.09258279388944618,"I_v2a4":5871593542678431,"Wa4_over_I_v2a4":0.006377902626168124,"Ra4_over_Y_I_v2a4":0.9998916579250956}}];
function requireCheck(ok,msg){if(!ok)throw new Error(msg);}
function frozen(L,A,B) {
  requireCheck(L>0 && B>=0 && A+B-L>0,"expanding initial branch required");
  const disc=Math.sqrt(A*A+4*L*B), x=(A+disc)/(2*L);
  return {x_turn:x,N_turn:Math.log(x)/4,tau_turn:Math.acos((2*L-A)/disc)/(4*Math.sqrt(L))};
}
function rk4Turn(L,A,B,steps) {
  const exact=frozen(L,A,B),dt=exact.tau_turn/steps;
  let x=1,p=4*Math.sqrt(A+B-L),maxInvariant=0;
  const accel=z=>8*A-16*L*z;
  for(let n=0;n<steps;n++){
    const k1x=p,k1p=accel(x);
    const k2x=p+dt*k1p/2,k2p=accel(x+dt*k1x/2);
    const k3x=p+dt*k2p/2,k3p=accel(x+dt*k2x/2);
    const k4x=p+dt*k3p,k4p=accel(x+dt*k3x);
    x+=dt*(k1x+2*k2x+2*k3x+k4x)/6;
    p+=dt*(k1p+2*k2p+2*k3p+k4p)/6;
    const inv=p*p-16*(-L*x*x+A*x+B);
    maxInvariant=Math.max(maxInvariant,Math.abs(inv)/(16*(L*x*x+Math.abs(A*x)+B)));
  }
  return {steps,x_at_analytic_turn:x,p_at_analytic_turn:p,N_turn_numeric:Math.log(x)/4,
    absolute_N_error:Math.abs(Math.log(x)/4-exact.N_turn),max_relative_invariant:maxInvariant};
}
function gateWindow(q,w,ell,threshold=.9){
  const c=1/threshold-1-Math.abs(w),discriminant=c*c-4*ell*q;
  if(c<0||discriminant<0)return null;
  const root=Math.sqrt(discriminant);
  const xlo=q===0?0:2*q/(c+root);
  const xhi=ell===0?Infinity:(c+root)/(2*ell);
  return {N_entry:xlo===0?null:Math.log(xlo)/4,N_exit:xhi===Infinity?null:Math.log(xhi)/4,
    interpretation:"source-free scalar-free gate; null endpoint means unbounded"};
}
const rows=plateaus.map(({label,tag,classification,plateau:p})=>{
  if(!p)return {label,tag,missing_plateau:true,old_classification:classification};
  const q=p.R_over_sigma/2,ell_total=Math.abs(p.Omega_vac)/p.Omega_r,r=p.r;
  const ell_linear=ell_total*(1+q),w_linear=r*(1+q);
  const F0=1/((1+q)*(1+r+ell_total));
  const fixed=frozen(ell_total,1/(1+q)+r,q/(1+q));
  const linear_one_efold_injection_over_Rs=ell_linear*Math.exp(4)/.1-1;
  return {label,tag,old_classification:classification,
    q_over_linear:q,vac_over_total:ell_total,vac_over_linear:ell_linear,W_over_linear:w_linear,
    F0_upper_full_absolute_fraction:F0,instant_gate_pass:F0>=.9,
    local_vac_frozen_N:fixed.N_turn,local_vac_grouped_a4_N:Math.log((1+r)/ell_total)/4,
    archived_measured_N_to_turnaround:p.efolds_to_turnaround,
    conditional_linear_one_efold_injection_over_Rs:linear_one_efold_injection_over_Rs,
    frozen_gate_window:gateWindow(q,w_linear,ell_linear)};
});
const negative={Omega_r:1/(1+.02-.99),F_abs:1/(1+.02+.99)};
const positive={Omega_r:1/(1+.02+.99),F_abs:1/(1+.02+.99)};
const highQ=frozen(.2,1,10),highQRK=[500,1000,2000].map(n=>rk4Turn(.2,1,10,n));
for(const r of highQRK)requireCheck(r.absolute_N_error<1e-8 && r.max_relative_invariant<1e-8,"quadratic-density ODE control failed");
requireCheck(Math.abs(negative.F_abs-positive.F_abs)<1e-15 && negative.Omega_r>1,"cancellation control");
const y2=rows.find(r=>r.label==="Y=2 dc=0.01");
requireCheck(Math.abs(y2.F0_upper_full_absolute_fraction-.7654491662636435)<1e-12,"linear radiation normalization");
function endpointInjection(Rs,a4s,a4f,Vf,alpha,eps,beta=1/36){
  const target=Math.abs(Vf)/eps;
  const Rreq=beta===0?target/alpha:2*target/(alpha+Math.sqrt(alpha*alpha+4*beta*target));
  return {R_required:Rreq,weighted_injection_required:Math.max(0,a4f*Rreq-a4s*Rs)};
}
const report={status:"PASS",scope:"post hoc frozen and saved-point diagnostics, not a new 5D evolution",
  definition:"linear radiation over radiation+quadratic+absolute vacuum/Weyl/scalar/source terms",
  plateaus:rows,controls:{negative_vacuum:negative,positive_vacuum:positive,high_quadratic_exact:highQ,
    high_quadratic_RK4:highQRK,high_quadratic_gate:gateWindow(1,.02,1e-5)},
  endpoint_inverse_example:{linear:endpointInjection(1,1,Math.exp(4),.2847183530796413,1,.1,0),
    quadratic:endpointInjection(1,1,Math.exp(4),.2847183530796413,1,.1)},
  caveats:["F0 is an upper bound on full F_abs, not a whole-trajectory verdict",
    "frozen continuation assumes no transfer or scalar force and fixed vacuum/tension/Weyl",
    "source requirement is an endpoint target with endpoint quantities, not sufficient era duration"]};
console.log(JSON.stringify(report,null,2));
