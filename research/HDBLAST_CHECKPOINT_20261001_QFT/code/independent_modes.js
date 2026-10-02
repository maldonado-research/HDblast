"use strict";
function independentRK4Modes(rows){
function at(t){let lo=0,hi=rows.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(rows[m][0]<=t)lo=m;else hi=m;}if(t===rows[hi][0])lo=Math.max(0,hi-1);const r=rows[lo],s=rows[lo+1],h=s[0]-r[0],z=(t-r[0])/h,z2=z*z,z3=z2*z;function hval(i,d){return[(2*z3-3*z2+1)*r[i]+(-2*z3+3*z2)*s[i]+h*(z3-2*z2+z)*r[d]+h*(z3-z2)*s[d],((6*z2-6*z)*r[i]+(-6*z2+6*z)*s[i])/h+(3*z2-4*z+1)*r[d]+(3*z2-2*z)*s[d]];}const n=hval(1,2),p=hval(3,4);return {a:Math.exp(n[0]),H:n[1],phi:p[0],v:p[1]};}
let l=1,r=2;for(let j=0;j<60;j++){const m=(l+r)/2;if(at(m).phi<.5)l=m;else r=m;}const cross=(l+r)/2,g=at(cross),q=100*Math.abs(g.v),ks=g.a*Math.sqrt(q),nodes=[];
for(let i=0;i<32;i++){let x=Math.cos(Math.PI*(i+.75)/(32+.5));for(let j=0;j<20;j++){let p0=1,p1=x;for(let n=2;n<=32;n++){const pn=((2*n-1)*x*p1-(n-1)*p0)/n;p0=p1;p1=pn;}const dp=32*(x*p1-p0)/(x*x-1),dx=p1/dp;x-=dx;if(Math.abs(dx)<2e-16)break;}nodes.push(.05+(x+1)*2.95/2);}nodes.sort((a,b)=>a-b);
const ids=[0,10,21,31],kappas=ids.map(i=>nodes[i]),ks4=kappas.map(x=>x*ks),runs=[];
for(const steps of [100000,200000]){const h=6.9/steps,ini=at(0),ys=ks4.map(k=>{const w=Math.sqrt(k*k/ini.a**2+10000*(ini.phi-.5)**2),wd=(-ini.H*k*k/ini.a**2+10000*(ini.phi-.5)*ini.v)/w,c=1/Math.sqrt(2*ini.a**3*w);return[c,0,(-1.5*ini.H-wd/(2*w))*c,-w*c];});let maxWr=0;
function f(y,k,b){const w2=k*k/b.a**2+10000*(b.phi-.5)**2;return[y[2],y[3],-3*b.H*y[2]-w2*y[0],-3*b.H*y[3]-w2*y[1]];}
for(let j=0;j<steps;j++){const t=j*h,b0=at(t),bm=at(t+h/2),b1=at((j+1)*h);for(let m=0;m<4;m++){const y=ys[m],k=ks4[m],f1=f(y,k,b0),f2=f(y.map((v,i)=>v+h*f1[i]/2),k,bm),f3=f(y.map((v,i)=>v+h*f2[i]/2),k,bm),f4=f(y.map((v,i)=>v+h*f3[i]),k,b1);for(let i=0;i<4;i++)y[i]+=h*(f1[i]+2*f2[i]+2*f3[i]+f4[i])/6;}if(j%1000===0||j===steps-1)for(const y of ys)maxWr=Math.max(maxWr,Math.abs(2*b1.a**3*(y[1]*y[2]-y[0]*y[3])-1));}
const end=at(6.9),result=ys.map((y,i)=>{const w=Math.sqrt(ks4[i]**2/end.a**2+10000*(end.phi-.5)**2);return {index:ids[i],kappa:kappas[i],k:ks4[i],n:end.a**3/(2*w)*((w*y[0]+y[3])**2+(w*y[1]-y[2])**2),wronskian:2*end.a**3*(y[1]*y[2]-y[0]*y[3])};});runs.push({steps,h,maximum_sampled_Wronskian_error:maxWr,rows:result});}
const changes=runs[0].rows.map((v,i)=>Math.abs(v.n-runs[1].rows[i].n));return {status:runs.every(r=>r.maximum_sampled_Wronskian_error<1e-6)&&changes.every(x=>x<1e-5)?"PASS":"FAIL",scope:"exploratory independent RK4 implementation of four registered primary nodes; primary DOP853 results are separate",crossing:{t:cross,...g,q,k_scale:ks},runs,absolute_step_refinement_changes:changes};
}
const rows=require('../data/independent_geometry.json');
console.log(JSON.stringify(independentRK4Modes(rows),null,2));
