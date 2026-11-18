// Validated initial-data shooting bracket. Build with IEEE binary64 and no fast-math.
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
using std::max;
using std::min;

struct Interval {
    double lo, hi;
    Interval(double v = 0) : lo(v), hi(v) {}
    Interval(double a, double b) : lo(a), hi(b) {
        if (!(a <= b) || !std::isfinite(a) || !std::isfinite(b))
            throw std::runtime_error("Invalid interval");
    }
};
using I = Interval;
double down(double x) { return std::nextafter(x,-INFINITY); }
double up(double x) { return std::nextafter(x,INFINITY); }
I decimal(double x) { return {down(x),up(x)}; }
I span(I a, I b) { return {min(a.lo,b.lo),max(a.hi,b.hi)}; }
I operator+(I a,I b) { return {down(a.lo+b.lo),up(a.hi+b.hi)}; }
I operator-(I a) { return {-a.hi,-a.lo}; }
I operator-(I a,I b) { return a+(-b); }
I operator*(I a,I b) {
    double w=a.lo*b.lo, x=a.lo*b.hi, y=a.hi*b.lo, z=a.hi*b.hi;
    return {down(min({w,x,y,z})),up(max({w,x,y,z}))};
}
I operator/(I a,I b) {
    if (b.lo<=0 && b.hi>=0) throw std::runtime_error("Division by zero");
    return a*I(down(1./b.hi),up(1./b.lo));
}
I sq(I a) {
    double m=max(std::abs(a.lo),std::abs(a.hi));
    double n=a.lo>0 ? a.lo : a.hi<0 ? -a.hi : 0;
    return {n==0?0:down(n*n),up(m*m)};
}
I root(I a) {
    if (a.lo<=0) throw std::runtime_error("Square root outside positive branch");
    return {down(std::sqrt(a.lo)),up(std::sqrt(a.hi))};
}
bool contains(I box,I x) { return box.lo<=x.lo && x.hi<=box.hi; }
I inflate(I x,double pad) { return x+I(-pad,pad); }

// Alternating exp(-x) series after reduction to x in [0,1/8]. Terms
// through odd order 15 bound from below; even order 14 bounds from above.
I exp_minus_point(double x) {
    if (x<0) throw std::runtime_error("Invalid exp argument");
    if (x>128) return {0,exp_minus_point(128).hi};
    int halves=0;
    while (x>0.125) { x*=0.5; ++halves; }
    I term=1, sum=1, even;
    for (int i=1;i<=15;i++) {
        term=term*I(-x)/i;
        sum=sum+term;
        if (i==14) even=sum;
    }
    I result(sum.lo,even.hi);
    for (int i=0;i<halves;i++) result=sq(result);
    return result;
}
I exp_minus(I x) {
    return {exp_minus_point(max(0.,x.hi)).lo,exp_minus_point(max(0.,x.lo)).hi};
}
// Endpoints s in [0,1]. f=exp(1-1/s), g=f/s^2. The g maximum is
// at s=1/2. Below s=1/128 both f and g increase, including at s=0.
I bump_value(double s,bool weighted) {
    if (s<=0) return I(0);
    if (s<1./128) return {0,bump_value(1./128,weighted).hi};
    I t=I(1)/s;
    I f=exp_minus(t-1);
    return weighted ? sq(t)*f : f;
}
struct Pulse { I f, dz; };
Pulse bump(I t) {
    I x=6-50*t;
    I s=1-sq(x);
    s=I(max(0.,s.lo),min(1.,max(0.,s.hi)));
    I f=span(bump_value(s.lo,false),bump_value(s.hi,false));
    I g=span(bump_value(s.lo,true),bump_value(s.hi,true));
    if (s.lo<=0.5 && s.hi>=0.5) g=span(g,bump_value(.5,true));
    return {f,-100*x*g};
}

I W(I f) { return 1-f+f*sq(f)/3; }
I W1(I f) { return sq(f)-1; }
I U(I f) { return sq(W1(f))/2-2*sq(W(f))/3; }
I U1(I f) { return 2*f*W1(f)-4*W(f)*W1(f)/3; }
I U2(I f) { return 2*W1(f)+4*sq(f)-4*(sq(W1(f))+2*f*W(f))/3; }
I potential_difference(I f,I h) {
    I dw=h*(W1(f)+f*h+sq(h)/3);
    I dw1=h*(2*f+h);
    return dw1*(2*W1(f)+dw1)/2-2*dw*(2*W(f)+dw)/3;
}

using State=std::array<I,4>; // E=exp(B0), Q=A0_Z, phi0, r=phi0_Z
State slope(State y) {
    auto [e,q,f,r]=y;
    return {-e*q,-1+sq(q)+sq(r)/6+sq(e)*U(f)/3,-r,-sq(e)*U1(f)+3*q*r};
}
State second(State y) {
    auto [e,q,f,r]=y;
    auto [de,dq,df,dr]=slope(y);
    return {-de*q-e*dq,2*q*dq+r*dr/3+(2*e*de*U(f)+sq(e)*U1(f)*df)/3,
            -dr,-2*e*de*U1(f)-sq(e)*U2(f)*df+3*dq*r+3*q*dr};
}
// Verify a Picard box covering the whole time step, followed by the
// second-order Taylor enclosure centered at the left endpoint.
State quiet_step(State y,I dt,State &box) {
    bool ok=false;
    State initial=slope(y);
    for (double factor=1.05;factor<10;factor*=2) {
        for (int j=0;j<4;j++) {
            double vmax=max(std::abs(initial[j].lo),std::abs(initial[j].hi));
            box[j]=inflate(y[j],factor*(vmax+.001)*dt.hi);
        }
        auto v=slope(box);
        ok=true;
        for (int j=0;j<4;j++) ok &= contains(box[j],y[j]+I(0,dt.hi)*v[j]);
        if (ok) break;
    }
    if (!ok) throw std::runtime_error("Quiet Picard box failed");
    State v=slope(y), a=second(box), next;
    for (int j=0;j<4;j++) next[j]=y[j]+dt*v[j]+sq(dt)*a[j]/2;
    return next;
}
struct Coeff { I a,b; };
Coeff mass_coeff(State box,Pulse bump,I eta) {
    auto [e,q,f,r]=box;
    I h=eta*bump.f, v=eta*bump.dz;
    I k=(q*(sq(v)+sq(r+v))+2*v*(r+v))/(3*(1-sq(q)));
    I b=2*sq(e)*(q+1)*(r*v+sq(v))/3+2*sq(sq(e))*q*potential_difference(f,h)/3;
    return {-k,-b}; // t=-Z increases towards the interior
}
I mass_step(I y,I dt,Coeff coeff,I &tube) {
    I initial=coeff.a*y+coeff.b;
    double vmax=max(std::abs(initial.lo),std::abs(initial.hi));
    for (double factor=1.05;factor<10;factor*=2) {
        tube=inflate(y,factor*(vmax+.001)*dt.hi);
        I v=coeff.a*tube+coeff.b;
        if (contains(tube,y+I(0,dt.hi)*v)) return y+dt*v;
    }
    throw std::runtime_error("Mass Picard box failed");
}
// Parameter derivatives on [0,.003] share the interval bounds of the
// linear mass equation. On a step, enclose D first, then D_eta, then
// D_eta_eta. Each previously enclosed tube bounds its forcing.
void variation_step(std::array<I,3> &d,I dt,State box,Pulse pulse,I eta) {
    auto [e,q,f,r]=box;
    I g=pulse.dz, h=eta*pulse.f;
    I a1=-2*(q+1)*(r*g+2*eta*sq(g))/(3*(1-sq(q)));
    I a2=-4*(q+1)*sq(g)/(3*(1-sq(q)));
    I factor=-2*sq(sq(e))*q/3;
    I b1=-2*sq(e)*(q+1)*(r*g+2*eta*sq(g))/3+factor*U1(f+h)*pulse.f;
    I b2=-4*sq(e)*(q+1)*sq(g)/3+factor*U2(f+h)*sq(pulse.f);
    Coeff c=mass_coeff(box,pulse,eta);
    I y,v1,v2;
    d[0]=mass_step(d[0],dt,c,y);
    d[1]=mass_step(d[1],dt,{c.a,a1*y+b1},v1);
    d[2]=mass_step(d[2],dt,{c.a,2*a1*v1+a2*y+b2},v2);
}
State wall() {
    I f=decimal(.0006170897672752272),c=decimal(.5975949350280132);
    I d=decimal(-3.106933495673783), delta=I(1)/10;
    I sigma=2*W(f)+delta*(1+c*f+d*sq(f)/2);
    I sig1=2*W1(f)+delta*(c+d*f);
    I rho=1/root(sq(sigma)/36-sq(sig1)/48+U(f)/6);
    return {rho,rho*sigma/6,f,-rho*sig1/2};
}
int main(int argc,char**argv) {
    static_assert(std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53,"Requires binary64");
    if (std::fegetround()!=FE_TONEAREST) throw std::runtime_error("Use round to nearest");
    int n=argc==2?std::stoi(argv[1]):16384; // cells per 0.01 radial distance
    if (n<32 || n>65536) throw std::runtime_error("Subdivision count out of range");
    int den=100*n;
    I dt=I(1)/den;
    State y=wall(), box;
    std::array<I,3> derivatives={0,0,0};
    std::array<I,5> defects={0,0,0,0,0};
    std::array<I,5> eta={I(1)/1000,I(3)/1000,span(I(1)/1000,I(3)/1000),
                         I(201)/100000,I(203)/100000};
    double q_min=INFINITY, mu_max=-INFINITY;
    const int early=(129*n+9)/10, witness=13*n;
    I phi_early, phi_late, phi_packet, x_packet;
    double x_quiet_upper=-INFINITY,r_quiet_min=INFINITY;
    for (int i=0;i<14*n;i++) {
        if (i==early) phi_early=y[2];
        y=quiet_step(y,dt,box);
        if (early<=i && i<witness) {
            x_quiet_upper=max(x_quiet_upper,(sq(box[3])/sq(box[0])).hi);
            r_quiet_min=min(r_quiet_min,box[3].lo);
        }
        q_min=min(q_min,box[1].lo);
        I base_mass=sq(box[0])*(1-sq(box[1]));
        if (i<10*n) { mu_max=max(mu_max,base_mass.hi); continue; }
        Pulse packet=bump(I(i,i+1)/den);
        for (int j=0;j<5;j++) {
            I tube;
            defects[j]=mass_step(defects[j],dt,mass_coeff(box,packet,eta[j]),tube);
            mu_max=max(mu_max,(base_mass+tube).hi);
        }
        variation_step(derivatives,dt,box,packet,span(I(0),eta[1]));
        if (i+1==witness) {
            Pulse sample=bump(I(witness)/den);
            phi_late=y[2];
            phi_packet=y[2]+eta[2]*sample.f;
            I mu0=sq(y[0])*(1-sq(y[1]));
            x_packet=(1+defects[2]/mu0)*(sq(y[3])+2*y[3]*eta[2]*sample.dz)/sq(y[0]);
        }
    }
    std::printf("cells per 0.01: %d\n",n);
    for (int j=0;j<5;j++)
        std::printf("defect[%d] %.17g %.17g\n",j,defects[j].lo,defects[j].hi);
    std::printf("Q lower bound %.17g; mu upper bound %.17g\n",q_min,mu_max);
    std::printf("quiet phi early [%.17g, %.17g]; late [%.17g, %.17g]\n",phi_early.lo,phi_early.hi,phi_late.lo,phi_late.hi);
    std::printf("packet phi [%.17g, %.17g]; quiet r lower %.17g\n",phi_packet.lo,phi_packet.hi,r_quiet_min);
    std::printf("quiet X upper %.17g; packet X [%.17g, %.17g]\n",x_quiet_upper,x_packet.lo,x_packet.hi);
    std::printf("endpoint second parameter derivative on [0,.003]: [%.17g, %.17g]\n", derivatives[2].lo, derivatives[2].hi);
    bool physical=phi_late.hi<phi_packet.lo && phi_packet.hi<phi_early.lo &&
                  r_quiet_min>0 && x_packet.lo>x_quiet_upper;
    return !(defects[3].lo>0 && defects[4].hi<0 && derivatives[2].hi<0 && physical && defects[0].lo>0 && defects[1].hi<0 && q_min>1 && mu_max<0);
}
