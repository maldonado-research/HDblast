"use strict";
function independentCoherenceControls(){
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]],add=(a,b)=>[a[0]+b[0],a[1]+b[1]],scale=(a,x)=>a.map(v=>v*x),norm=a=>a[0]*a[0]+a[1]*a[1],dot=(a,b)=>a[0]*b[0]+a[1]*b[1];
let max=0,maxn=0,maxW=0,rows=0;
for(const a of [.3,1,4])for(const H of [-.2,0,.8])for(const k of [.2,2,8])for(const n of [0,.03,2])for(const theta of [.2,1.7]){
const G=3,x=.4,vphi=-.7,p2=(k/a)**2,m2=(G*x)**2,om=Math.sqrt(p2+m2),od=(-H*p2+G*G*x*vphi)/om,b=1/Math.sqrt(2*a**3*om),al=[Math.sqrt(1+n),0],be=[Math.sqrt(n)*Math.cos(theta),Math.sqrt(n)*Math.sin(theta)],f=scale(add(al,be),b),vel=scale(mul([0,-om],add(al,scale(be,-1))),b),acc=add(scale(vel,-3*H),scale(f,-om*om)),C=al[0]*be[0],E=.5*(norm(vel)+om*om*norm(f)),P=.5*(norm(vel)-(p2/3+m2)*norm(f)),J=G*G*x*norm(f);
const E2=om*(n+.5)/a**3,P2=p2*(n+.5)/(3*a**3*om)-(2*p2/3+m2)*C/(a**3*om),J2=G*G*x*(n+.5+C)/(a**3*om);
max=Math.max(max,...[E-E2,P-P2,J-J2].map(Math.abs));const dE=dot(vel,acc)+om*od*norm(f)+om*om*dot(f,vel),ndot=a**3/om*(dE+(3*H-od/om)*E),ndot2=(3*H+od/om)*C;maxn=Math.max(maxn,Math.abs(ndot-ndot2));const W=2*a**3*(f[1]*vel[0]-f[0]*vel[1]);maxW=Math.max(maxW,Math.abs(W-1));rows++;
}
const boundary=[];for(const H of [.1,.8])for(const L of [2,5]){const rho=L**4/(16*Math.PI**2),pressure=rho/3,continuity=3*H*(rho+pressure),flux=H*L**4/(4*Math.PI**2);boundary.push({H,Lambda:L,relative_error:Math.abs(continuity-flux)/flux});}
if(max>1e-10||maxn>1e-10||maxW>1e-12||boundary.some(x=>x.relative_error>1e-14))throw Error("Independent identity control failed");
return {status:"PASS",scope:"independent algebraic/coherence and moving-band controls; no archived quantum-mode run",rows,max_stress_current_absolute_error:max,max_occupation_derivative_absolute_error:maxn,max_Wronskian_error:maxW,conformal_moving_physical_cutoff:boundary};
}
console.log(JSON.stringify(independentCoherenceControls(),null,2));
