function run() {
"use strict";
// Exact rational controls for an ordinary reflected Duhamel integral identity.
// No PDE/source premise, physical observable, nonlinear family or novelty claim.
function R(n,d=1n){n=BigInt(n);d=BigInt(d);if(!d)throw Error("zero denominator");if(d<0n){n=-n;d=-d;}let a=n<0n?-n:n,b=d;while(b){let t=a%b;a=b;b=t;}return[n/a,d/a];}
const z=R(0),o=R(1),h=R(1,2),I=R;
const add=(a,b)=>R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const neg=a=>[-a[0],a[1]], sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>R(a[0]*b[0],a[1]*b[1]), div=(a,b)=>R(a[0]*b[1],a[1]*b[0]);
const eq=(a,b)=>a[0]===b[0]&&a[1]===b[1],sgn=a=>a[0]<0n?-1:a[0]>0n?1:0;
const pow=(a,n)=>{let b=o;while(n--)b=mul(b,a);return b;};
const sum=a=>a.reduce(add,z),str=a=>a[1]===1n?String(a[0]):a[0]+"/"+a[1];
const choose=(n,k)=>{let a=1;for(let j=1;j<=k;j++)a=a*(n-j+1)/j;return a;};
function need(ok,label){if(!ok)throw Error(label);}
function profile(p=[],hinges=[]){return{p:p.map(v=>Array.isArray(v)?v:I(v)),hinges:hinges.map(([c,e,k])=>[I(c),I(e),Array.isArray(k)?k:I(k)])};}
function val(P,t,order=0){
 let out=z;
 for(let n=0;n<P.p.length;n++){let den=1;for(let j=1;j<=order;j++)den*=n+j;out=add(out,div(mul(P.p[n],pow(t,n+order)),I(den)));}
 for(const[c,e,k]of P.hinges){let u=sub(mul(e,t),k);if(sgn(u)>0){let den=order===0?1:order===1?2:6;let term=mul(c,pow(u,order+1));if(order===1)term=mul(term,e);out=add(out,div(term,I(den)));}}
 return out;
}
const D=(v,dx=z,dr=z)=>[v,dx,dr], dc=a=>D(a),da=(a,b)=>a.map((v,i)=>add(v,b[i])),dn=a=>a.map(neg),ds=(a,b)=>da(a,dn(b));
const dm=(a,b)=>[mul(a[0],b[0]),add(mul(a[1],b[0]),mul(a[0],b[1])),add(mul(a[2],b[0]),mul(a[0],b[2]))];
const scale=(a,b)=>dm(dc(a),b);
function pv(P,t,order){return[val(P,t[0],order),mul(val(P,t[0],order-1),t[1]),mul(val(P,t[0],order-1),t[2])];}
function kernel(F,G,x0,r0,J,force=null){
 const x=D(x0,o,z),r=D(r0,z,o),a=da(r,x),q=ds(r,x),t=dn(q),PF=u=>pv(F,u,1),QF=u=>pv(F,u,2),PG=u=>pv(G,u,1),QG=u=>pv(G,u,2);
 if(force===true||(force===null&&sgn(q[0])>=0)){
  return scale(h,ds(da(ds(dm(x,PF(a)),scale(h,ds(QF(a),QF(q)))),scale(h,ds(QG(a),QG(q)))),dm(x,PG(q))));
 }
 let physical=ds(da(ds(ds(dm(x,PF(a)),scale(h,QF(a))),scale(h,QF(t))),scale(h,da(QG(a),QG(dn(t))))),dm(r,PG(dn(t))));
 let reflected=da(ds(dm(t,PF(t)),QF(t)),scale(h,da(QG(t),QG(dn(t)))));
 return scale(h,da(physical,scale(I(J),reflected)));
}
function charDerivs(F,G,x,r,J){
 const a=add(r,x),q=sub(r,x),t=neg(q);
 let p=add(mul(x,val(F,a)),mul(h,sub(val(G,a,1),val(G,q,1)))),m;
 if(sgn(q)>=0)m=add(mul(h,sub(val(F,a,1),val(F,q,1))),mul(x,val(G,q)));
 else m=add(add(mul(h,sub(val(F,a,1),val(F,t,1))),mul(r,val(G,neg(t)))),mul(I(J),add(mul(t,val(F,t)),mul(h,sub(val(G,t,1),val(G,neg(t),1))))));
 return[mul(h,add(p,m)),mul(h,sub(p,m))];
}
// Independent exact polygon integration. Coordinates are (s,v), CCW.
// Green: integral s^i v^j ds dv = contour s^(i+1)v^j dv/(i+1).
function conv(a,b){let c=Array(a.length+b.length-1).fill(z);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=add(c[i+j],mul(a[i],b[j]));return c;}
function polyPower(a,n){let p=[o];while(n--)p=conv(p,a);return p;}
function moment(P,i,j){
 if(P.length<3)return z;let ans=z;
 for(let k=0;k<P.length;k++){let a=P[k],b=P[(k+1)%P.length],ds0=sub(b[0],a[0]),dv=sub(b[1],a[1]);
 let pol=conv(polyPower([a[0],ds0],i+1),polyPower([a[1],dv],j));
 ans=add(ans,mul(dv,sum(pol.map((v,n)=>div(v,I(n+1))))));}
 return div(ans,I(i+1));
}
function clip(P,as,av,c){
 const v=p=>add(add(mul(as,p[0]),mul(av,p[1])),c),out=[];
 if(!P.length)return out;
 for(let i=0;i<P.length;i++){let a=P[i],b=P[(i+1)%P.length],u=v(a),w=v(b),ina=sgn(u)>=0,inb=sgn(w)>=0;
 if(ina)out.push(a);if(ina!==inb){let t=div(u,sub(u,w));out.push([add(a[0],mul(t,sub(b[0],a[0]))),add(a[1],mul(t,sub(b[1],a[1])))]);}}
 return out;
}
function ridgeIntegral(T,as,av,P){
 let ans=z;
 for(let n=0;n<T.p.length;n++)for(let i=0;i<=n;i++)ans=add(ans,mul(mul(T.p[n],mul(I(choose(n,i)),mul(pow(as,i),pow(av,n-i)))),moment(P,i,n-i)));
 for(const[c,e,k]of T.hinges){let bs=mul(e,as),bv=mul(e,av),C=clip(P,bs,bv,neg(k));
 ans=add(ans,mul(c,sub(add(mul(bs,moment(C,1,0)),mul(bv,moment(C,0,1))),mul(k,moment(C,0,0)))));}
 return ans;
}
function directK(F,G,x,r,J){
 const P=[[z,sub(r,x)],[x,r],[z,add(r,x)]],pos=clip(P,z,o,z),negative=clip(P,z,neg(o),z);
 const plus=add(ridgeIntegral(F,o,o,pos),ridgeIntegral(G,neg(o),o,pos));
 const minus=add(ridgeIntegral(F,o,neg(o),negative),ridgeIntegral(G,neg(o),neg(o),negative));
 return mul(h,add(plus,mul(I(J),minus)));
}
function intProfile(P,alpha,beta,l,u){
 if(eq(alpha,z))return mul(sub(u,l),val(P,beta));
 return div(sub(val(P,add(mul(alpha,u),beta),1),val(P,add(mul(alpha,l),beta),1)),alpha);
}
function endpointDerivs(F,G,x,r,J){
 const a=add(r,x),q=sub(r,x),t=neg(q);
 const upper=add(intProfile(F,z,a,z,x),intProfile(G,I(-2),a,z,x));let lower;
 if(sgn(q)>=0)lower=add(intProfile(F,I(2),q,z,x),intProfile(G,z,q,z,x));
 else lower=add(mul(I(J),add(intProfile(F,z,t,z,t),intProfile(G,I(-2),t,z,t))),add(intProfile(F,I(2),q,t,x),intProfile(G,z,q,t,x)));
 return[mul(h,add(upper,lower)),mul(h,sub(upper,lower))];
}
const zero=profile(),one=profile([1]),sources=[];
for(let n=0;n<=4;n++){let p=Array(n+1).fill(0);p[n]=1;sources.push([profile(p),zero],[zero,profile(p)]);}
sources.push([profile([2,-3,1,2]),profile([-1,4,-2,3])]);
const Fp=profile([],[[1,1,1],[2,1,I(5,4)]]),Gp=profile([],[[3,1,1],[-2,-1,1]]);
sources.push([Fp,Gp],[profile([1,-2],[[2,1,1]]),profile([-3,1],[[-1,-1,1]])]);
let integral=0,derivative=0,characteristic=0,joins=0,boundary=0,initial=0,wedge=0;
const points=[];for(let i=0;i<=6;i++)for(let j=0;j<=7;j++)points.push([I(i,4),I(j,4)]);
points.push([I(3,2),I(1,4)],[I(5,7),I(2,9)],[I(7,8),I(7,8)]);
for(let si=0;si<sources.length;si++)for(const J of[-1,1]){
 const[F,G]=sources[si];
 for(const[x,r]of points){
  let k=kernel(F,G,x,r,J),v=directK(F,G,x,r,J),e=endpointDerivs(F,G,x,r,J),c=charDerivs(F,G,x,r,J),tag=si+":"+J+":"+str(x)+":"+str(r);
  need(eq(k[0],v),"integral "+tag);integral++;
  need(eq(k[1],e[0])&&eq(k[2],e[1]),"endpoint derivative "+tag);derivative+=2;
  need(eq(c[0],e[0])&&eq(c[1],e[1]),"characteristic derivative "+tag);characteristic+=2;
 }
 for(let i=0;i<=12;i++){const t=I(i,8),kp=kernel(F,G,t,t,J,true),kn=kernel(F,G,t,t,J,false);need(kp.every((v,j)=>eq(v,kn[j])),"C1 join");joins+=3;
 const k=kernel(F,G,t,z,J);if(J===-1){need(eq(k[0],z)&&eq(k[1],z),"Dirichlet shell");boundary+=2;}else{need(eq(k[2],z),"Neumann shell");boundary++;}
 const ki=kernel(F,G,z,t,J);need(ki.every(v=>eq(v,z)),"zero Cauchy data");initial+=3;}
}
for(const J of[-1,1])for(let i=0;i<=8;i++)for(let j=0;j+i<=8;j++){const x=I(i,8),r=I(j,8),k=kernel(Fp,Gp,x,r,J);need(k.every(v=>eq(v,z))&&eq(directK(Fp,Gp,x,r,J),z),"P1 exact zero wedge");wedge++;}
const oddCorrect=kernel(one,zero,o,z,-1)[0],oddWrong=kernel(one,zero,o,z,1)[0];
need(eq(oddCorrect,z)&&!eq(oddWrong,oddCorrect),"wrong reflection rejected");
const correct=directK(one,zero,o,I(2),1),missingHalf=mul(I(2),correct);need(!eq(correct,missingHalf),"missing factor half rejected");
return{status:"PASS_EXACT_REFLECTED_RIDGE_KERNEL_CONTROLS",arithmetic:"normalized BigInt fractions; exact polygon clipping and Green monomial integration",source_pairs:sources.length,points_per_source_parity:points.length,exact_integral_equalities:integral,autodifferentiated_first_derivative_equalities:derivative,characteristic_first_derivative_equalities:characteristic,C1_join_equalities:joins,shell_equalities:boundary,initial_zero_equalities:initial,continuous_P1_zero_wedge_cases:wedge,negative_controls:{wrong_reflection:{correct:str(oddCorrect),wrong:str(oddWrong),rejected:true},missing_retarded_half:{correct:str(correct),wrong:str(missingHalf),rejected:true}},limits:["Ordinary integral identity and finite algebra controls only.","No retained PDE/source premise or evolved total-response certificate is established.","No physical observable, nonlinear family or novelty claim."]};

}
process.stdout.write(JSON.stringify(run(), null, 2) + '\n');
