"use strict";
const fs=require("node:fs"),path=require("node:path"),crypto=require("node:crypto");
const RATIONAL=(function(){
 function R(n,d=1n){n=BigInt(n);d=BigInt(d);if(d===0n)throw Error("zero denominator");if(d<0n){n=-n;d=-d;}let a=n<0n?-n:n,b=d;while(b){const t=a%b;a=b;b=t;}return[n/a,d/a];}
 const q=s=>{if(Array.isArray(s))return s;const a=String(s).split("/");return R(a[0],a[1]??"1");};
 const add=(a,b)=>{a=q(a);b=q(b);return R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);};
 const neg=a=>{a=q(a);return[-a[0],a[1]];},sub=(a,b)=>add(a,neg(b));
 const mul=(a,b)=>{a=q(a);b=q(b);return R(a[0]*b[0],a[1]*b[1]);};
 const div=(a,b)=>{a=q(a);b=q(b);return R(a[0]*b[1],a[1]*b[0]);};
 const eq=(a,b)=>{a=q(a);b=q(b);return a[0]===b[0]&&a[1]===b[1];};
 const lt=(a,b)=>{a=q(a);b=q(b);return a[0]*b[1]<b[0]*a[1];};
 const max=(a,b)=>lt(a,b)?q(b):q(a),abs=a=>{a=q(a);return[a[0]<0n?-a[0]:a[0],a[1]];};
 const sum=a=>a.reduce(add,R(0)),str=a=>{a=q(a);return a[1]===1n?String(a[0]):a[0]+"/"+a[1];};
 function isqrt(n){if(n<0n)throw Error("negative square");if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2)),y=(x+n/x)>>1n;while(y<x){x=y;y=(x+n/x)>>1n;}return x;}
 const sqrtUp=a=>{a=q(a);if(a[0]===0n)return R(0);const g=1n<<112n,v=(a[0]*g*g+a[1]-1n)/a[1];return R(isqrt(v)+1n,g);};
 return{R,q,add,sub,mul,div,eq,lt,max,abs,sum,str,sqrtUp};
})();

function audit(r,data,E,reportSHA) {
 const F=RATIONAL, {q,add,mul,div,eq,lt,max,sum,sqrtUp}=F;
 let checks=0;
 const need=(ok,s)=>{if(!ok)throw Error(s);checks++;};
 const Q=q("3409/8000000"),U=q("1030000000000000000000000000000000"),S=q("52162547172519/2"),one=q("1"),G=[1,2,1,3,1];
 const mNorm=m=>sqrtUp(sum(m.flatMap((row,i)=>row.map((v,j)=>mul(div(G[i],G[j]),mul(v,v))))));
 for(const name of ["A","B","D","M"])need(eq(mNorm(data.matrices[name]),r.constant_matrix_H_operator_upper[name]),"matrix norm "+name);
 for(const level of r.levels){
  const errors=[q(0),q(0),q(0)];let prior=q(0);
  for(const cell of level.cell_certificates){
   need(eq(cell.r[0],prior)&&lt(cell.r[0],cell.r[1]),"continuous cell cover");prior=q(cell.r[1]);
   for(let k=0;k<3;k++){
    const h=sqrtUp(sum(cell.component_error_caps[k].map((c,i)=>mul(G[i],mul(c,c)))));
    need(eq(h,cell.H_error_caps[k]),"weighted cell norm");
    errors[k]=max(errors[k],h);
   }
  }
  need(eq(prior,"7/4")&&level.cells===level.cell_certificates.length,"full profile cover");
  const L=level.ledger;
  for(let k=0;k<3;k++)need(eq(errors[k],L[["epsilon_u","epsilon_du","epsilon_v"][k]]),"max profile norm");
  const e0=sum([L.e_free_retained,L.epsilon_v,S]),x0=add(mul(add(one,Q),U),e0);
  const ce=r.coefficient_error_H_operator_upper,nm=r.constant_matrix_H_operator_upper;
  const dfc=mul(sum([ce.A,ce.B,mul("3/4",ce.D)]),x0),dgc=mul(mul("1/4",ce.M),x0);
  const dfp=add(mul(add(nm.A,nm.B),L.epsilon_du),mul(nm.D,L.epsilon_u)),dgp=mul(nm.M,L.epsilon_u);
  const df=add(dfc,dfp),dg=add(dgc,dgp),e1=add(mul("3/2",df),dg),tail=mul(mul(Q,Q),U),total=sum([tail,mul(add(one,Q),e0),e1]);
  for(const[name,v]of Object.entries({e0,B0_X_upper:x0,delta_f_coefficient:dfc,delta_g_coefficient:dgc,delta_f_profile:dfp,delta_g_profile:dgp,delta_f:df,delta_g:dg,e1,two_term_tail:tail,total_X_error_upper:total}))need(eq(v,L[name]),"complete ledger "+name);
  need(lt(total,r.target_X_error)===level.passes_target,"actual target");
 }
 const last=r.levels.at(-1);
 need(last.passes_target===r.passes_target&&eq(last.ledger.total_X_error_upper,E.total_X_error_upper),"evaluated response bound");
 need(last.nodes===data.nodes.length&&last.cells+1===data.nodes.length,"selected footprint");
 const edge=[last.cell_certificates[0].r[0],...last.cell_certificates.map(c=>c.r[1])];
 for(let i=0;i<edge.length;i++)need(eq(edge[i],data.nodes[i]),"stored common knots");
 for(const name of ["u","du","v"]){
  need(data.profiles[name].length===data.nodes.length,"profile row coverage");
  for(const row of data.profiles[name]){need(row.length===5,"all five means");for(const c of row)q(c);}
 }
 return {status:"PASS_INDEPENDENT_BIGINT_NATIVE_LEDGER_AUDIT",exact_checks:checks,
  method:"Independent BigInt fractions and integer-square-root norms; every cell cover, norm, maximum, coefficient/profile/forcing/tail join, stored knot and target decision.",
  report_sha256:reportSHA,total_X_error_upper:last.ledger.total_X_error_upper,
  accuracy_target_passed:r.passes_target,failed_meshes:r.levels.filter(l=>!l.passes_target).map(l=>l.mesh_denominator),
  selected_mesh:last.mesh_denominator,
  scope:"Saved native certificate arithmetic and representation-shape audit; not an independent reexecution of interval profile or coefficient producers."};
}

const root=path.resolve(process.argv[2]||path.join(__dirname,"results"));
const raw=fs.readFileSync(path.join(root,"REPORT.json"),"utf8");
const result=audit(JSON.parse(raw),
 JSON.parse(fs.readFileSync(path.join(root,"PROFILES.json"),"utf8")),
 JSON.parse(fs.readFileSync(path.join(root,"EVALUATED_CANDIDATE.json"),"utf8")),
 crypto.createHash("sha256").update(raw).digest("hex"));
process.stdout.write(JSON.stringify(result,null,2)+"\n");
