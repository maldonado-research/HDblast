function auditInterference(jumps,af){
 const eta=jumps.map(x=>x.eta),d=jumps.map(x=>x.generic_first_jump);
 const spacing=Math.min(...eta.slice(1).map((x,i)=>x-eta[i]));
 if(!(spacing>0))throw Error("Nonpositive knot spacing");
 const lower=1e3/spacing,upper=1e6/spacing;
 function ci(x){
   let f=1/x,g=1/(x*x),tf=f,tg=g;
   for(let n=1;n<=5;n++){tf*=-(2*n)*(2*n-1)/(x*x);tg*=-(2*n+1)*(2*n)/(x*x);f+=tf;g+=tg;}
   return Math.sin(x)*f-Math.cos(x)*g;
 }
 const sumsq=d.reduce((s,x)=>s+x*x,0),diagonal=sumsq*Math.log(upper/lower);
 let cross=0,absoluteWeightedPairs=0;
 for(let i=0;i<d.length;i++)for(let j=i+1;j<d.length;j++){
   const gap=eta[j]-eta[i];absoluteWeightedPairs+=Math.abs(d[i]*d[j])/gap;cross+=2*d[i]*d[j]*(ci(2*gap*upper)-ci(2*gap*lower));
 }
 const ratio=(diagonal+cross)/diagonal;
 return {method:"Executed JavaScript: six-term large-argument Ci expansion; scipy cross-check supplied",knots:d.length,min_conformal_knot_spacing:spacing,K:[lower,upper],minimum_Ci_argument:2*spacing*lower,diagonal_integral:diagonal,interference_integral:cross,coefficient_ratio:ratio,phase_independent_cross_bound:2*(1/lower+1/upper)*absoluteWeightedPairs,phase_independent_relative_bound:2*(1/lower+1/upper)*absoluteWeightedPairs/diagonal,registered_relative_tolerance:1e-3,passed:Math.abs(ratio-1)<=1e-3,energy_integral_leading:(diagonal+cross)/(32*Math.PI**2*af**4),scope:"Integral of leading archive jump amplitudes only; no exact archive high-k mode evolution or physical energy",roundoff_caution:"Floating-point asymptotic diagnostic, not an interval certificate."};
}
// Reproduce with Node.js from the checkpoint directory.
const fs=require("node:fs");
const jumps=JSON.parse(fs.readFileSync("outputs/primary_hermite_jumps.json","utf8"));
const geom=JSON.parse(fs.readFileSync("outputs/geometry_audit.json","utf8"));
console.log(JSON.stringify(auditInterference(jumps,geom.reconstructions[0].a_final),null,2));
