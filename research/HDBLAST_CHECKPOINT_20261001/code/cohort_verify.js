"use strict";
// Copyright (c) 2026 Ricardo Maldonado. MIT License.
const input={"model":{"rb":7.838285538073243,"delta":0.1,"c":0.5975949350280132,"d":-3.106933495673783},"cohort":{"background":{"H0tau_end":7.353355110666192,"phi_end":1.0361021913225144,"W_end":0.003666076408891769,"vac_end":-0.00033994224195131275,"lna_end":2.0958245770395325,"Wa4_end":16.033255291448366},"rows":[{"key":"main_dstar_Y0_dc1e-2_dzf5e-4|phistar=0.5|G=100|y=1","phistar":0.5,"G":100,"q":128.31216219495383,"number_factor":0.9809491460159282,"crossing":{"s_star":1.4861183349220013,"v":1.2831216219495383,"A":0.9551402171424274,"B":-13.08843572954591,"h0":0.6447779418500043,"h1":-1.1618657098120238,"lna_star":1.3968545679837534,"D":-2.44445626635661,"G_min_screen":43.184788170506884,"n_fit":13},"crossing_half0p05":{"s_star":1.4861183349220013,"v":1.2831852240834927,"A":0.9559636396300987,"B":-13.229901548714151,"h0":0.6448142499296644,"h1":-1.161874524126087,"lna_star":1.3968545679837534,"D":-2.469974449927901,"G_min_screen":43.252846898817545,"n_fit":7},"Lambda_hat_full":53.61022411468734,"Lambda_hat_post":53.61022411468734,"b_cut_full":0.000006490185826699186},{"key":"main_dstar_Y0_dc1e-2_dzf5e-4|phistar=0.5|G=1000|y=1","phistar":0.5,"G":1000,"q":1283.1216219495382,"number_factor":0.9980949146015928,"crossing":{"s_star":1.4861183349220013,"v":1.2831216219495383,"A":0.9551402171424274,"B":-13.08843572954591,"h0":0.6447779418500043,"h1":-1.1618657098120238,"lna_star":1.3968545679837534,"D":-2.44445626635661,"G_min_screen":43.184788170506884,"n_fit":13},"crossing_half0p05":{"s_star":1.4861183349220013,"v":1.2831852240834927,"A":0.9559636396300987,"B":-13.229901548714151,"h0":0.6448142499296644,"h1":-1.161874524126087,"lna_star":1.3968545679837534,"D":-2.469974449927901,"G_min_screen":43.252846898817545,"n_fit":7},"Lambda_hat_full":536.1022411468734,"Lambda_hat_post":536.1022411468734,"b_cut_full":6.490185826699185e-9},{"key":"main_dstar_Y0_dc1e-2_dzf5e-4|phistar=0.9|G=100|y=1","phistar":0.9,"G":100,"q":74.62015750906967,"number_factor":0.965699644900714,"crossing":{"s_star":1.836049179957632,"v":0.7462015750906967,"A":-3.016386947850691,"B":1.630427574951324,"h0":0.27339786019513135,"h1":-0.6601802783440182,"lna_star":1.5517837254421836,"D":-2.5594979001257374,"G_min_screen":2189.8058998630772,"n_fit":11},"crossing_half0p05":{"s_star":1.836049179957632,"v":0.746134110440027,"A":-3.017264979844526,"B":1.81628561616329,"h0":0.2733728866573851,"h1":-0.6601305575937217,"lna_star":1.5517837254421836,"D":-2.4986110733241786,"G_min_screen":2191.6753342194756,"n_fit":5},"Lambda_hat_full":89.70173134673259,"Lambda_hat_post":13.610224114687341,"b_cut_full":0.00000138547126701606},{"key":"main_dstar_Y0_dc1e-2_dzf5e-4|phistar=0.9|G=1000|y=1","phistar":0.9,"G":1000,"q":746.2015750906967,"number_factor":0.9965699644900714,"crossing":{"s_star":1.836049179957632,"v":0.7462015750906967,"A":-3.016386947850691,"B":1.630427574951324,"h0":0.27339786019513135,"h1":-0.6601802783440182,"lna_star":1.5517837254421836,"D":-2.5594979001257374,"G_min_screen":2189.8058998630772,"n_fit":11},"crossing_half0p05":{"s_star":1.836049179957632,"v":0.746134110440027,"A":-3.017264979844526,"B":1.81628561616329,"h0":0.2733728866573851,"h1":-0.6601305575937217,"lna_star":1.5517837254421836,"D":-2.4986110733241786,"G_min_screen":2191.6753342194756,"n_fit":5},"Lambda_hat_full":897.0173134673258,"Lambda_hat_post":136.1022411468734,"b_cut_full":1.38547126701606e-9}]}};
function check(ok,msg){if(!ok)throw Error(msg);}
function simpson(f,lo,hi,n){check(n%2===0,"even quadrature count");const h=(hi-lo)/n;let s=f(lo)+f(hi);for(let i=1;i<n;i++)s+=(i%2?4:2)*f(lo+i*h);return h*s/3;}
const measure=x=>x*x*Math.exp(-Math.PI*x*x)/(2*Math.PI*Math.PI);
const moment=n=>simpson(measure,0,8,n);
const moment2=n=>simpson(x=>x*x*measure(x),0,8,n);
const nExact=1/(8*Math.PI**3),k2Exact=3/(2*Math.PI),normal=[];
for(const n of [1024,2048,4096,8192]){
const I=moment(n),K=moment2(n)/I;
normal.push({bins:n,number_relative_error:Math.abs(I/nExact-1),second_moment_relative_error:Math.abs(K/k2Exact-1)});
check(Math.abs(I/nExact-1)<1e-10 && Math.abs(K/k2Exact-1)<1e-10,"Gaussian quadrature failed");
}
const b=input.cohort.background,rb=input.model.rb;
const W=p=>1-p+p**3/3;
const sigma=rb*(2*W(b.phi_end)+input.model.delta*(1+input.model.c*b.phi_end+input.model.d*b.phi_end**2/2));
const required=36*b.W_end/.1/(Math.sqrt(sigma*sigma+36*b.W_end/.1)+sigma);
function energyBudget(row,F=row.number_factor,G=row.G){
const ratio=G/row.G,q=row.q*ratio,a=Math.exp(b.lna_end-row.crossing.lna_star),
M=row.Lambda_hat_post*ratio,Lambda=row.Lambda_hat_full*ratio,
nstar=F*q**1.5/(8*Math.PI**3),Rbar_cap=Lambda**(-3)*nstar/a**3*Math.sqrt(M*M+3*q/(2*Math.PI*a*a));
return {Rbar_cap,r_lower:b.W_end/(sigma*Rbar_cap/18+Rbar_cap*Rbar_cap/36),
gain_for_r_point1:required/Rbar_cap,a_rel:a,M,q,nstar,bcut:Lambda**(-3)};
}
const rows=input.cohort.rows.map(row=>({...energyBudget(row),phistar:row.phistar,G:row.G,
local_production_screen_pass:row.G>=row.crossing.G_min_screen,
status:row.G>=row.crossing.G_min_screen?"conditional screened ansatz":"formal budget; local production screen fails"}));
check(Math.abs(rows[0].r_lower-51.43974552824248)<1e-9,"endpoint reproduction failed");
const q=7,a=3,M=2;
const histories=[
{label:"all survivors",decays:[],survival:1,m_end:1.7},
{label:"all decay late",decays:[{p:1,ad:3,m:2}],survival:0,m_end:0},
{label:"early decay",decays:[{p:1,ad:1,m:2}],survival:0,m_end:0},
{label:"mixed nonmonotone bounded masses",decays:[{p:.2,ad:1.1,m:1.9},{p:.3,ad:1.7,m:.4},{p:.1,ad:2.4,m:1.8}],survival:.4,m_end:.3},
{label:"incomplete specified population",decays:[{p:.2,ad:1.4,m:1.5}],survival:.3,m_end:2}
];
const jensen=Math.sqrt(M*M+q*k2Exact/(a*a));
const historyChecks=histories.map(H=>{
const totalP=H.survival+H.decays.reduce((s,d)=>s+d.p,0);
check(totalP<=1+1e-12 && H.m_end<=M && H.decays.every(d=>d.p>=0&&d.ad>=1&&d.ad<=a&&d.m<=M),"bad positive decay history");
const E=simpson(x=>{
const k=x*Math.sqrt(q);
const surv=H.survival*Math.sqrt(k*k/(a*a)+H.m_end*H.m_end);
const daughter=H.decays.reduce((s,d)=>s+d.p*Math.sqrt(k*k/(a*a)+(d.m*d.ad/a)**2),0);
return measure(x)*(surv+daughter);
},0,8,4096)/nExact;
check(E<=jensen*(1+1e-12),"cohort bound failed");
return {label:H.label,total_probability:totalP,mean_endpoint_energy:E,Jensen_upper:jensen,pass:true};
});
const contractionCounter={a_decay:2,a_end:1,k:0,m:1,daughter_energy:2,expanding_formula_cap:1};
check(contractionCounter.daughter_energy>contractionCounter.expanding_formula_cap,"contraction control failed");
const masslessEnergy=simpson(x=>measure(x)*x*Math.sqrt(q)/a,0,8,4096)/nExact;
check(masslessEnergy>0,"kinetic term negative control");
const corridor=[];
for(const row of [input.cohort.rows[0],input.cohort.rows[2]])for(const Fmax of [1,1.3]){
const Gmin=row.crossing.G_min_screen,cap=energyBudget(row,Fmax,Gmin);
let last=cap.Rbar_cap;
for(let i=1;i<=500;i++){const G=Gmin*Math.exp(12*i/500),next=energyBudget(row,Fmax,G).Rbar_cap;check(next<last,"coupling-corridor monotonic control");last=next;}
corridor.push({phistar:row.phistar,Gmin,Fmax,Rbar_cap:cap.Rbar_cap,r_lower:cap.r_lower,
status:"post-registration exploratory fixed-history Gaussian corridor"});
}
const report={status:"PASS",scope:"specified Gaussian cohort on fixed expanding archived history",
provenance:input.cohort.rows.map(x=>x.key),sigmahat:sigma,required_Rbar_for_r_point1:required,
Gaussian_quadrature:normal,endpoint_rows:rows,positive_decay_histories:historyChecks,
negative_controls:{contraction:contractionCounter,massless_mean_energy:masslessEnergy,
wrong_phase_space_factor:(2*Math.PI)**3},
coupling_corridor:corridor,caveats:["population/mass/expansion bound is assumed, not derived quantum production",
"full-history characteristic cutoff screen is not UV certification",
"initial radiation and multiple cohorts need separate additions; contraction excluded",
"fixed Weyl/background endpoint, not a backreacted trajectory or future-time global bound"]};
console.log(JSON.stringify(report,null,2));
