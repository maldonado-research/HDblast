// Axisymmetric pair scattering on an expanding wall. Use one reflected ensemble for
// both signs. Each molecular-chaos pair update represents two charge-mirrored events.
#include <algorithm>
#include "spectra.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

constexpr double pi=3.14159265358979323846, mV2=8300.*8300., width=.30, Q=1e7, Kwall=2.54;
constexpr double majorant=48.*pi/(mV2-8.); // bound on m² sigma times Moller speed
struct Vector {double x,y,z;};
Vector operator+(Vector a,Vector b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vector operator-(Vector a,Vector b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
Vector operator*(Vector a,double b){return {a.x*b,a.y*b,a.z*b};}
double dot(Vector a,Vector b){return a.x*b.x+a.y*b.y+a.z*b.z;}
Vector cross(Vector a,Vector b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm(Vector a){return std::sqrt(dot(a,a));}
struct Particle {double qperp,qz; int hits;};
struct Point {double t,a,Iref,D,vref,e,Gref,Bref,Xref;};
struct Birth {double t,a,pstar,c,d,target;};
struct Moments {double v=0.,e=0.,pz=0.,zz=0.,trace=0.,unhit=0.;};
struct RNG {
    std::mt19937_64 r;
    explicit RNG(unsigned seed):r(seed) {}
    double unit(){return ((r()>>11)+.5)*0x1.0p-53;}
    int index(int n){return std::min(n-1,static_cast<int>(unit()*n));}
};
struct Wall {
    double w,density_scale,weight,interaction;
    bool respond;
    std::vector<Particle> ps;
    RNG random;
    double I=0.,Y=0.,Z=0.,work=0.,redshift=0.,injection=0.,event_residual=0., exact_work=0., exact_red=0.;
    long long hits=0, proposed=0;
    Moments eval(double a,double impulse,bool detail) const {
        Moments m;
        for(const auto& p:ps) {
            const double z=(p.qz+impulse)/a, b=p.qperp/a, e=std::sqrt(1.+z*z+b*b);
            m.v+=z/e;
            m.trace+=(z*z+b*b)/e;
            if(detail) {m.e+=e;m.zz+=z*z/e;m.pz+=z;m.unhit+=(p.hits==0);}
        }
        return m;
    }
    double dy(const Point& p,double y,double z,double v_sum) const {
        return -4.*Kwall*z/p.a-2.*y/Kwall-Q*(w*v_sum-p.D*p.vref)/p.a;
    }
    void ledger(const Particle& p,double al,double ar,double il,double ir) {
        const double am=.5*(al+ar), z0=p.qz+il,z1=p.qz+ir;
        const double q0=p.qperp*p.qperp+z0*z0,q1=p.qperp*p.qperp+z1*z1;
        const double l=1./(al*al),m=1./(am*am),r=1./(ar*ar);
        const double em0=std::sqrt(1.+q0*m),em1=std::sqrt(1.+q1*m);
        exact_work+=w*(z0+z1)*(ir-il)*m/(em0+em1);
        exact_red+=w*(q0*(l-m)/(std::sqrt(1.+q0*l)+em0)+q1*(m-r)/(em1+std::sqrt(1.+q1*r)));
    }
    void advance(const Point& left,const Point& mid,const Point& right) {
        const double h=right.t-left.t;
        const double y_mid=respond ? Y+.5*h*dy(left,Y,Z,eval(left.a,I,false).v) : 0.;
        const double z_mid=Z+.5*h*Y/(left.a*Kwall);
        const double i_mid=I+mid.Iref-left.Iref+.5*h*Q*Y/(left.a*Kwall);
        const auto mm=eval(mid.a,i_mid,false);
        const double dI=(right.Iref-left.Iref)+h*Q*y_mid/(mid.a*Kwall);
        if(respond) Y+=h*dy(mid,y_mid,z_mid,mm.v);
        Z+=h*y_mid/(mid.a*Kwall);
        work+=w*mm.v*dI/mid.a;
        redshift+=w*mm.trace*(right.a-left.a)/mid.a;
        for(const auto& p:ps) ledger(p,left.a,right.a,I,I+dI);
        I+=dI;
    }
    void birth(const Birth& b,double a,double impulse) {
        double u;
        do {
            u=2.*random.unit()-1.;
            if(b.d==0. || b.target<=0.) break;
        } while(random.unit()*b.target>std::abs(b.c-b.d*u));
        ps.push_back({a*b.pstar*std::sqrt(1.-u*u),a*b.pstar*u-impulse,0});
        injection+=w*std::hypot(1.,b.pstar);
    }
    void collide(double a,double dt) {
        const int n=ps.size();
        if(n<2 || interaction==0.) return;
        const double density=density_scale*w*n/(a*a*a);
        const double rate=.5*n*density*majorant*interaction;
        double now=0.;
        for(;;) {
            now-=std::log(random.unit())/rate;
            if(now>dt) break;
            const int i=random.index(n);
            int j=random.index(n-1); if(j>=i) ++j;
            const auto pa=ps[i],pb=ps[j];
            const double angle=2.*pi*random.unit();
            const Vector p={pa.qperp/a,0.,(pa.qz+I)/a};
            const Vector q={pb.qperp*std::cos(angle)/a,pb.qperp*std::sin(angle)/a,-(pb.qz+I)/a};
            const double e=std::sqrt(1.+dot(p,p)), f=std::sqrt(1.+dot(q,q));
            // Avoid subtracting near-parallel large momenta when computing s-4.
            const double pq=dot(p,q),prod=e*f;
            const Vector cr=cross(p,q);
            // (E1 E2)^2 - (1+p.q)^2 = |p-q|^2 + |p x q|^2.
            const double q2=pq>0. ? 2.*(dot(p-q,p-q)+dot(cr,cr))/(prod+1.+pq) : 2.*(prod-pq-1.);
            const double s=4.+std::max(0.,q2), beta2=q2/s, g=width*(s+2.);
            const double k=6*pi*std::sqrt(std::max(0.,beta2))*g*g/((s-mV2)*(s-mV2)+beta2*g*g)/prod;
            ++proposed;
            if(k>majorant*(1.+1e-12)) throw std::runtime_error("Collision majorant failed");
            if(random.unit()*majorant>=k) continue;
            const Vector total=p+q;
            const double W=std::sqrt(s), kstar=.5*std::sqrt(q2);
            Vector axis=p-total*((e+.5*W)/(e+f+W));
            axis=axis*(1./norm(axis));
            Vector side=cross(axis,std::abs(axis.z)<.8 ? Vector{0,0,1}:Vector{1,0,0});
            side=side*(1./norm(side));
            double c0, A=3.-2.*beta2,B=beta2*beta2;
            do {c0=2.*random.unit()-1.;} while(random.unit()*(A+B)>A+B*c0*c0);
            const double phi=2*pi*random.unit(), sn=std::sqrt(1.-c0*c0);
            const Vector scat=(axis*c0+side*(sn*std::cos(phi))+cross(axis,side)*(sn*std::sin(phi)))*kstar;
            const Vector pnew=scat+total*(.5+dot(scat,total)/(W*(e+f+W))),qnew=total-pnew;
            const double er=std::abs((std::sqrt(1.+dot(pnew,pnew))+std::sqrt(1.+dot(qnew,qnew)))/(e+f)-1.);
            event_residual=std::max(event_residual,er);
            Vector dest={0,0,1};
            // Equal branching into three axes. A different final flavor accelerates
            // along its own perpendicular electric field after this collision.
            if(weight<1. && random.unit()>1./3.) {const double ph=2*pi*random.unit();dest={std::cos(ph),std::sin(ph),0};}
            const double pz=dot(pnew,dest),qz=dot(qnew,dest);
            ps[i]={a*norm(cross(pnew,dest)),a*pz-I,pa.hits+1};
            ps[j]={a*norm(cross(qnew,dest)),-a*qz-I,pb.hits+1};
            ++hits;
        }
    }
    double phase_area=0., phase_time=0., phase_pressure=0.;
    void snapshot(const Point& p,int ix,double dt=0.) {
        auto m=eval(p.a,I,true); int n=ps.size();
        if(n<1) throw std::runtime_error("No particles at marker");
        const double v=m.v/n, z=m.pz/n, E=m.e/n, pz=m.zz/n, tr=m.trace/n;
        const double gamma=1./std::sqrt((1.-v)*(1.+v));
        const double own=gamma*gamma*gamma*(E-2*v*z+v*v*pz);
        // Landau rest frame has zero energy flux. Its proper number density
        // is gamma_L*n_lab*(1-v_L*v_N) even if N and energy flow differ.
        const double l=2*z/(E+pz+std::sqrt(std::max(0.,(E+pz-2*z)*(E+pz+2*z))));
        const double gamL=1./std::sqrt((1.-l)*(1.+l));
        const double ownL=gamL*(E-2*l*z+l*l*pz)/(1.-l*v);
        const double d=w*n/(p.a*p.a*p.a), newgas=d*E;
        const double xdiff=Y/(p.a*p.a*Kwall), bdiff=2*Z/(p.a*p.a);
        const double adddot=p.e-p.Gref+newgas+Kwall*(p.Xref*xdiff+.5*xdiff*xdiff+p.Bref*bdiff+.5*bdiff*bdiff);
        // Replace gas and electric/magnetic wall energies on the prescribed background.
        if(dt>0.) {
            phase_area+=dt*d*(ownL-1.)/adddot; phase_time+=dt; phase_pressure+=dt*pz/tr;
            return;
        }
        std::vector<distribution::Observation> list;
        for(const auto& particle:ps) list.push_back({particle.qperp/p.a,(particle.qz+I)/p.a});
        const auto shape=distribution::evaluate(list,v,own);
        std::cout<<ix<<' '<<p.t<<' '<<n<<' '<<v<<' '<<l<<' '<<E<<' '<<own<<' '<<ownL
                 <<' '<<d*(own-1.)<<' '<<d*(ownL-1.)<<' '<<adddot<<' '<<xdiff<<' '
                 <<w*m.e<<' '<<work<<' '<<redshift<<' '<<injection<<' '<<event_residual<<' '
                 <<double(hits)*2/n<<' '<<m.unhit/n<<' '<<pz/tr<<' '<<double(proposed)<<' '<<exact_work<<' '<<exact_red<<' '<<shape.temperature<<' '<<shape.ks<<' '<<shape.longitudinal_share
                 <<' '<<shape.tube<<' '<<shape.thermal_tube<<' '<<(phase_time>0.?phase_area/phase_time:d*(ownL-1.)/adddot)
                 <<' '<<(phase_time>0.?phase_pressure/phase_time:pz/tr)<<'\n';
        phase_time=phase_area=phase_pressure=0.;
    }
};
int main(int argc,char** argv) {
    if(argc!=6)throw std::runtime_error("trace seed feedback interaction flavor_weight");
    std::ifstream st(argv[1],std::ios::binary);
    auto read=[&](auto& a){st.read(reinterpret_cast<char*>(&a),sizeof a);};
    int32_t nrow,n,markers; double w,density;
    read(nrow);read(n);read(markers);read(w);read(density);
    std::vector<Point> rows(nrow);
    for(auto& p:rows) read(p);
    std::vector<Birth> births(n);
    for(auto& b:births) read(b);
    std::vector<int32_t> indices(markers);for(auto& a:indices)read(a);
    if(!st||nrow%4!=1)throw std::runtime_error("Bad trace header");
    const double br=std::stod(argv[5]);
    Wall wall{w,density,br,std::stod(argv[4])*br,std::stoi(argv[3])!=0,{},RNG(std::stoul(argv[2]))};
    wall.I=rows[0].Iref;
    size_t born=0,mark=0;
    std::cout<<std::setprecision(16);
    for(int i=0;i+4<nrow;i+=4) {
        const auto& mid=rows[i+2];
        wall.advance(rows[i],rows[i+1],mid);
        // Births represented at their own kinetic momenta; correct for transport
        // from their true time to the operator-split midpoint.
        while(born<births.size() && births[born].t<=rows[i+4].t) {
            const auto& b=births[born++];
            const double frac=(b.t-rows[i].t)/(rows[i+4].t-rows[i].t);
            const double tide=(1.-frac)*rows[i].Iref+frac*rows[i+4].Iref;
            const double I_birth=wall.I + tide-mid.Iref + (b.t-mid.t)*Q*wall.Y/(mid.a*Kwall);
            wall.birth(b,b.a,I_birth);
            wall.ledger(wall.ps.back(),b.a,mid.a,I_birth,wall.I);
        }
        wall.collide(mid.a,rows[i+4].t-rows[i].t);
        wall.advance(mid,rows[i+3],rows[i+4]);
        if(mark<indices.size()) {
            if(mark>=2 && rows[i+4].t>rows[indices[mark]].t-.25)
                wall.snapshot(rows[i+4],mark,rows[i+4].t-std::max(rows[i].t,rows[indices[mark]].t-.25));
            if(i+4==indices[mark]) {wall.snapshot(rows[i+4],mark);++mark;}
        }
    }
    if(mark!=indices.size()) throw std::runtime_error("Missing markers");
}
