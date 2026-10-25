#!/usr/bin/env node
"use strict";
/**
 * Independent exact audit of the actual finite B1 samples.
 * Usage: node independent_actual_b1.js NATIVE_OUTPUT_DIRECTORY
 * Reads REPORT.json, PROFILES.json, EVALUATED_CANDIDATE.json.
 * Uses built-in Node modules only; no project-kernel imports.
 */
const fs = require("fs");
const path = require("path");
const crypto = require("crypto");
const digest = value => crypto.createHash("sha256").update(value).digest("hex");
const audit = function audit(data,evaluated,report,inputHashes,checkerSha256){
const {R,r,z,vzero,va,vs,scale,mv,ma,ms,mn,mx}=(() => {
const gcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){const t=a%b;a=b;b=t;}return a;};
class R{
constructor(a=0n,b=1n){if(a instanceof R){this.n=a.n;this.d=a.d;return;}if(typeof a==="string"){const p=a.split("/");a=BigInt(p[0]);b=p.length>1?BigInt(p[1]):b;}else a=BigInt(a);b=BigInt(b);if(b===0n)throw Error("zero denominator");if(b<0n){a=-a;b=-b;}const g=gcd(a,b);this.n=a/g;this.d=b/g;}
add(x){x=new R(x);const g=gcd(this.d,x.d);return new R(this.n*(x.d/g)+x.n*(this.d/g),(this.d/g)*x.d);}
neg(){return new R(-this.n,this.d);}
sub(x){return this.add(new R(x).neg());}
mul(x){x=new R(x);const g=gcd(this.n,x.d),h=gcd(x.n,this.d);return new R((this.n/g)*(x.n/h),(this.d/h)*(x.d/g));}
div(x){x=new R(x);return this.mul(new R(x.d,x.n));}
pow(i){if(i<0)return new R(this.d,this.n).pow(-i);return new R(this.n**BigInt(i),this.d**BigInt(i));}
cmp(x){x=new R(x);const t=this.n*x.d-x.n*this.d;return t<0n?-1:t>0n?1:0;}
eq(x){return this.cmp(x)===0;}
toString(){return String(this.n)+(this.d===1n?"":"/"+this.d);}
}
const z=new R(0),r=x=>new R(x),vzero=()=>Array.from({length:5},()=>r(0)),va=(a,b)=>a.map((v,i)=>v.add(b[i])),vs=(a,b)=>a.map((v,i)=>v.sub(b[i])),scale=(c,a)=>a.map(v=>v.mul(c)),mv=(m,v)=>m.map(row=>row.reduce((a,c,j)=>a.add(c.mul(v[j])),r(0))),ma=(a,b)=>a.map((row,i)=>va(row,b[i])),ms=(a,b)=>a.map((row,i)=>vs(row,b[i]));
const mn=(a,b)=>a.cmp(b)<0?a:b,mx=(a,b)=>a.cmp(b)>0?a:b;
return {R,r,z,vzero,va,vs,scale,mv,ma,ms,mn,mx};
})();
let joins=0;
const need=(condition,label)=>{joins++;if(!condition)throw Error(label);};
need(data.mean_order.join(",")==="K0,L0,N0,J0,w0","mean order");
need(data.gram.join(",")==="1,2,1,3,1","Gram");
need(r(data.ell).eq("1/1000000"),"ell");
const nodes=data.nodes.map(r),par=[1,-1,1,1,1];
const mats=Object.fromEntries(Object.entries(data.matrices).map(([k,m])=>[k,m.map(row=>row.map(r))]));
const profiles=Object.fromEntries(Object.entries(data.profiles).map(([k,m])=>[k,m.map(row=>row.map(r))]));
const plus=ma(mats.A,mats.B),minus=ms(mats.B,mats.A);
const segments=[];let pbase=vzero();
for(let i=0;i<nodes.length-1;i++){
 const a=nodes[i],b=nodes[i+1],width=b.sub(a);
 const affine=name=>{const slope=scale(r(1).div(width),vs(profiles[name][i+1],profiles[name][i]));return [vs(profiles[name][i],scale(a,slope)),slope,vzero()];};
 const u=affine("u"),du=affine("du"),v=affine("v");
 const pv=[va(vs(pbase,scale(a,profiles.v[i])),scale(a.pow(2).div(2),v[1])),v[0],scale(r(1).div(2),v[1])];
 const sources=(u,du,v,pv)=>({
 F:[0,1,2].map(j=>scale(r(1).div(2),va(mv(plus,va(du[j],v[j])),mv(mats.D,va(u[j],pv[j]))))),
 G:[0,1,2].map(j=>scale(r(1).div(2),va(mv(minus,vs(du[j],v[j])),mv(mats.D,vs(u[j],pv[j])))))
 });
 segments.push({a,b,...sources(u,du,v,pv)});
 const reflect=(poly,baseSign)=>poly.map((row,j)=>row.map((x,c)=>x.mul(baseSign*par[c]*(j%2?-1:1))));
 segments.push({a:b.neg(),b:a.neg(),...sources(reflect(u,1),reflect(du,-1),reflect(v,1),reflect(pv,-1))});
 pbase=va(pbase,scale(width.div(2),va(profiles.v[i],profiles.v[i+1])));
}
function clip(poly,alpha,beta,bound,upper){
 const value=p=>p[0].mul(alpha).add(p[1].mul(beta)).sub(bound);
 const inside=v=>upper?v.cmp(0)<=0:v.cmp(0)>=0;
 const out=[];
 for(let i=0;i<poly.length;i++){
  const a=poly[i],b=poly[(i+1)%poly.length],fa=value(a),fb=value(b),ia=inside(fa),ib=inside(fb);
  if(ia)out.push(a);
  if(ia!==ib){const t=fa.div(fa.sub(fb));out.push([a[0].add(b[0].sub(a[0]).mul(t)),a[1].add(b[1].sub(a[1]).mul(t))]);}
 }
 return out;
}
function moments(poly,alpha,beta){
 let m0=r(0),m1=r(0),m2=r(0);
 for(let i=0;i<poly.length;i++){
  const a=poly[i],b=poly[(i+1)%poly.length],cross=a[0].mul(b[1]).sub(a[1].mul(b[0]));
  const ya=a[0].mul(alpha).add(a[1].mul(beta)),yb=b[0].mul(alpha).add(b[1].mul(beta));
  m0=m0.add(cross.div(2));m1=m1.add(cross.mul(ya.add(yb)).div(6));m2=m2.add(cross.mul(ya.pow(2).add(ya.mul(yb)).add(yb.pow(2))).div(12));
 }
 return [m0,m1,m2];
}
let polygonPieces=0,linePieces=0;
function integrateArea(poly,alpha,beta,name){
 let out=vzero();
 if(poly.length<3)return out;
 const vals=poly.map(p=>p[0].mul(alpha).add(p[1].mul(beta)));
 const lo=vals.reduce(mn),hi=vals.reduce(mx);
 for(const seg of segments){
  if(seg.b.cmp(lo)<=0||seg.a.cmp(hi)>=0)continue;
  if(seg[name].every(row=>row.every(c=>c.eq(0))))continue;
  const cut=clip(clip(poly,alpha,beta,seg.a,false),alpha,beta,seg.b,true);
  if(cut.length<3)continue;
  const mom=moments(cut,alpha,beta);
  for(let j=0;j<3;j++)out=va(out,scale(mom[j],seg[name][j]));
  polygonPieces++;
 }
 return out;
}
function integrateLine(lo,hi,slope,intercept,name){
 let out=vzero();if(lo.cmp(hi)>=0)return out;
 for(const seg of segments){
  let left=lo,right=hi;
  if(slope.eq(0)){if(intercept.cmp(seg.a)<0||intercept.cmp(seg.b)>0)continue;if(intercept.eq(seg.b)&&seg.b.cmp(nodes.at(-1))!==0)continue;}
  else{
   const aa=seg.a.sub(intercept).div(slope),bb=seg.b.sub(intercept).div(slope);
   left=mx(left,mn(aa,bb));right=mn(right,mx(aa,bb));
   if(left.cmp(right)>=0)continue;
  }
  const d1=right.sub(left),d2=right.pow(2).sub(left.pow(2)),d3=right.pow(3).sub(left.pow(3));
  const moms=[d1,slope.mul(d2).div(2).add(intercept.mul(d1)),slope.pow(2).mul(d3).div(3).add(slope.mul(intercept).mul(d2)).add(intercept.pow(2).mul(d1))];
  for(let j=0;j<3;j++)out=va(out,scale(moms[j],seg[name][j]));
  linePieces++;
 }
 return out;
}
function primitive(name,t,order){
 const negative=t.cmp(0)<0;const tt=negative?t.neg():t;
 if(order===0){
  for(let i=0;i<nodes.length-1;i++)if(tt.cmp(nodes[i])>=0&&tt.cmp(nodes[i+1])<=0){
   const w=tt.sub(nodes[i]).div(nodes[i+1].sub(nodes[i]));
   const v=va(scale(r(1).sub(w),profiles[name][i]),scale(w,profiles[name][i+1]));
   return negative?v.map((c,j)=>c.mul(par[j])):v;
  }
  throw Error("profile footprint");
 }
 let out=vzero();
 for(let i=0;i<nodes.length-1;i++){
  const a=nodes[i],b=mn(tt,nodes[i+1]);if(b.cmp(a)<=0)continue;
  const width=b.sub(a),left=profiles[name][i],slope=scale(r(1).div(nodes[i+1].sub(a)),vs(profiles[name][i+1],left));
  for(let j=0;j<order;j++){
   const bin=order===3&&j===1?2:1;
   const fact=order===3?2:1;
   const coeff=tt.sub(a).pow(order-1-j).mul(bin*(j%2?-1:1)).div(fact);
   out=va(out,scale(coeff,va(scale(width.pow(j+1).div(j+1),left),scale(width.pow(j+2).div(j+2),slope))));
  }
 }
 return negative?out.map((c,j)=>c.mul((order%2?-1:1)*par[j])):out;
}
const matches=[];let equalities=0;
for(const sample of evaluated.samples){
 const x=r(sample.x),rr=r(sample.r),a=rr.add(x),q=rr.sub(x),t=x.sub(rr);
 const cone=[[r(0),rr.sub(x)],[x,rr],[r(0),rr.add(x)]];
 const pos=clip(cone,0,1,r(0),false),neg=clip(cone,0,1,r(0),true);
 let kval=va(integrateArea(pos,1,1,"F"),integrateArea(pos,-1,1,"G"));
 const negIntegral=va(integrateArea(neg,1,-1,"F"),integrateArea(neg,-1,-1,"G"));
 kval=scale(r(1).div(2),va(kval,negIntegral.map((c,j)=>c.mul(par[j]))));
 const p=va(integrateLine(r(0),x,r(0),a,"F"),integrateLine(r(0),x,r(-2),a,"G"));
 const split=mx(r(0),t);
 let m=va(integrateLine(split,x,r(2),q,"F"),integrateLine(split,x,r(0),q,"G"));
 if(t.cmp(0)>0){const negm=va(integrateLine(r(0),t,r(0),t,"F"),integrateLine(r(0),t,r(-2),t,"G"));m=va(m,negm.map((c,j)=>c.mul(par[j])));}
 let out={value:kval,dx:scale(r(1).div(2),va(p,m)),dr:scale(r(1).div(2),vs(p,m))};
 if(t.cmp(0)>0){
  const PN=v=>v.map((c,j)=>j===1?r(0):c);
  const g=mv(mats.M,PN(va(primitive("u",t,0),primitive("v",t,1))));
  const pg=mv(mats.M,PN(va(primitive("u",t,1),primitive("v",t,2))));
  out={value:vs(out.value,pg),dx:vs(out.dx,g),dr:va(out.dr,g)};
 }
 for(const name of ["value","dx","dr"])for(let j=0;j<5;j++){
  equalities++;
  if(!out[name][j].eq(sample.feedback_exact[name][j]))throw Error("independent polygon/line feedback mismatch "+[sample.x,sample.r,name,j].join(","));
 }
 matches.push({x:sample.x,r:sample.r,exact_feedback_match:true});
}
let sampleJoins=0,noncollapsedDirichletBoxes=0;
const requireSample=(v,label)=>{sampleJoins++;if(!v)throw Error(label);};
const E=r(report.levels.at(-1).ledger.total_X_error_upper),ell=r(data.ell);
requireSample(evaluated.total_X_error_upper===E.toString(),"total E");
requireSample(evaluated.sample_count===6&&evaluated.samples.length===6,"six samples");
requireSample(evaluated.native_report_sha256===inputHashes["REPORT.json"],"native report hash join");
requireSample(evaluated.profile_data_sha256===inputHashes["PROFILES.json"],"profile hash join");
for(let s=0;s<evaluated.samples.length;s++){
 const sample=evaluated.samples[s],free=data.free_wave_samples[s];
 requireSample(sample.x===free.x&&sample.r===free.r,"sample coordinates");
 for(const name of ["value","dx","dr"])for(let j=0;j<5;j++)for(let side=0;side<2;side++){
  const expected=r(free.B0[name][j][side]).add(sample.feedback_exact[name][j]);
  requireSample(expected.eq(sample.complete_candidate_interval_evaluation[name][j][side]),"candidate sum");
  const widen=E.mul(name==="value"?3:1).mul(side===0?-1:1);
  requireSample(expected.add(widen).eq(sample.conditional_exact_response_component_enclosure[name][j][side]),"response widening");
 }
 for(let j=0;j<5;j++)for(let side=0;side<2;side++){
  requireSample(r(sample.complete_candidate_interval_evaluation.dx[j][side]).div(ell).eq(sample.complete_candidate_interval_evaluation.ordinary_t_velocity[j][side]),"candidate ordinary t");
  requireSample(r(sample.conditional_exact_response_component_enclosure.dx[j][side]).div(ell).eq(sample.conditional_exact_response_component_enclosure.ordinary_t_velocity[j][side]),"response ordinary t");
 }
 if(r(sample.r).eq(0))for(const name of ["value","dx"]){
  const pair=sample.complete_candidate_interval_evaluation[name][1].map(r);
  requireSample(pair[0].cmp(0)<=0,"Dirichlet lower enclosure");
  requireSample(pair[1].cmp(0)>=0,"Dirichlet upper enclosure");
  if(!pair[0].eq(pair[1]))noncollapsedDirichletBoxes++;
 }
 if(r(sample.x).add(sample.r).cmp(1)<=0)for(const name of ["value","dx","dr"])for(let j=0;j<5;j++)for(let side=0;side<2;side++)requireSample(r(sample.complete_candidate_interval_evaluation[name][j][side]).eq(0),"zero wedge evaluation");
}
need(E.cmp(10n**29n)<0,"target");
need(E.cmp("16090323000000000000000000000")<0,"outward displayed E");
return {
 schema:"HDBLAST_INDEPENDENT_ACTUAL_B1_POLYGON_AUDIT_V1",
 status:"PASS_INDEPENDENT_EXACT_POLYGON_CANDIDATE_AUDIT",
 checker_sha256:checkerSha256,
 input_sha256:inputHashes,
 method:"Independent exact polygon clipping and signed triangle moments for the complete physical ridge source, full J reflection on negative-r source, exact characteristic endpoint line integrals for dx/dr, and direct P1 segment-moment shell lifts.",
 imports_project_kernel:false,
 profile_nodes:nodes.length,
 polygon_pieces:polygonPieces,
 endpoint_line_pieces:linePieces,
 exact_feedback_component_equalities:equalities,
 candidate_enclosure_and_time_joins:sampleJoins,
 representation_and_target_checks:joins,
 dirichlet_interval_boxes_noncollapsed:noncollapsedDirichletBoxes,
 dirichlet_interval_note:"The defined candidate satisfies exact Dirichlet parity. Some reported interval boxes retain conservative same-expression cancellation width; every boundary L0 value/dx box contains zero.",
 total_X_error_upper:E.toString(),
 strict_upper_display:"16090323000000000000000000000",
 target_X_error:"100000000000000000000000000000",
 exact_error_below_display_and_target:true,
 samples:matches,
 scope:[
 "Recomputes all six frozen B1 samples and associated interval, whole-cone-error and ordinary-t arithmetic joins from the stored profile/matrix data.",
 "Does not independently regenerate the initial rational profile curvature/quantization bounds or reauthenticate the retained PDE/source closures.",
 "This checker uses JavaScript BigInt exact arithmetic. It does not run the native Python producer or verify GitHub Actions status.",
 "Five conditional C1 means only; no full angular waveform, pointwise Einstein residual, invariant observable value/sign, nonlinear family, finite-amplitude remainder or physical-origin claim."
 ]
};
};
if (process.argv.length !== 3) throw Error("usage: node independent_actual_b1.js NATIVE_OUTPUT_DIRECTORY");
const names = ["REPORT.json","PROFILES.json","EVALUATED_CANDIDATE.json"];
const expected = {
 "REPORT.json":"6a7ee36f95d2227405dd301ff5dba5e52ebf1b45b7944fd864310e04fa4749be",
 "PROFILES.json":"6d6f1c8b6fc11589255aa7031bc90b2612fbc022ab42046ea1114ee44ccda0e2",
 "EVALUATED_CANDIDATE.json":"1c0d6f62f3858b595335be8a6cd71966a76ce15fa2b921e36f951d1452d54cb6"
};
const raw = Object.fromEntries(names.map(name => [name,fs.readFileSync(path.join(process.argv[2],name))]));
const hashes = Object.fromEntries(names.map(name => [name,digest(raw[name])]));
for (const name of names) if (hashes[name] !== expected[name]) throw Error("frozen input hash "+name);
const result = audit(JSON.parse(raw["PROFILES.json"]),JSON.parse(raw["EVALUATED_CANDIDATE.json"]),JSON.parse(raw["REPORT.json"]),hashes,digest(fs.readFileSync(__filename)));
process.stdout.write(JSON.stringify(result,null,2)+"\n");
