// Particle, gyro rotation and source-free Fourier evolution from HDblast #196, particles/patch.cpp.
// Wave energy here concerns the daughter-induced perturbation only.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

constexpr double Q=1e7, K=2.54, pi=3.14159265358979323846;
constexpr int max_modes=7;
using Field=std::array<double,2*max_modes>;
struct Particle {
    double x,px,py,pz, impulse=0;
    Field basis{};
    double by=0,bz=0,ey=0,ez=0;
};
struct Modes { Field ey{},ez{},by{},bz{}; };
struct Seed {double t,a,pz,pt,phi;};
struct Step {double a0,am,a1,deta,x,b,e_wave,b_wave;int births;};
double norm(double a,double x,double y,double z) {return std::sqrt(a*a+x*x+y*y+z*z);}
double dot(const Field& a,const Field& b,int n) {
    double sum=0.; for(int j=0;j<2*n;++j) sum+=a[j]*b[j]; return sum;
}
double field_energy(const Modes& f,int n) {
    return .25*K*(dot(f.ey,f.ey,n)+dot(f.ez,f.ez,n)+dot(f.by,f.by,n)+dot(f.bz,f.bz,n));
}
void waves(Modes& f,double phase,int n) {
    // For each polarization dE/deta=sign*dB/dx, dB/deta=sign*dE/dx.
    auto pair=[&](Field& e,Field& b,int sign) {
        for(int j=0;j<n;++j) {
            int c=2*j,s=c+1;
            double co=std::cos((j+1)*phase),si=sign*std::sin((j+1)*phase);
            double ec=e[c],es=e[s],bc=b[c],bs=b[s];
            e[c]=co*ec+si*bs;e[s]=co*es-si*bc;
            b[c]=co*bc+si*es;b[s]=co*bs-si*ec;
        }
    };
    pair(f.ez,f.by,1);pair(f.ey,f.bz,-1);
}
void gyro(Particle& p,double a,double by,double bz,double duration) {
    // Cayley rotation of the comoving momentum at fixed particle energy.
    const double scale=.5*Q*duration/norm(a,p.px,p.py,p.pz);
    const double ty=scale*by,tz=scale*bz,s=2./(1.+ty*ty+tz*tz);
    const double x=p.px+p.py*tz-p.pz*ty, y=p.py-p.px*tz,z=p.pz+p.px*ty;
    p.px+=s*(y*tz-z*ty);p.py-=s*x*tz;p.pz+=s*x*ty;
}
void basis(Particle& p,int n) {
    const double c=std::cos(p.x),s=std::sin(p.x); double a=1.,b=0.;
    for(int j=0;j<n;++j) {
        const double older=a;a=a*c-b*s;b=b*c+older*s;
        p.basis[2*j]=a;p.basis[2*j+1]=b;
    }
}
// Symmetric split for the homogeneous Fourier-Maxwell equation with outgoing electric leakage.
double propagate(Modes& f,double k,double h,double a,int n){
    const double damping=std::exp(-a*h/K);
    double loss=0.;
    for(int pass=0;pass<2;++pass){
        loss+=.25*K*(-std::expm1(-2.*a*h/K))*(dot(f.ey,f.ey,n)+dot(f.ez,f.ez,n));
        for(int j=0;j<2*n;++j){f.ey[j]*=damping; f.ez[j]*=damping;}
        if(pass==0)waves(f,k*h,n);
    }
    return loss;
}

int main(int argc,char** argv) {
    if(argc!=7)throw std::runtime_error("trace retained_harmonics spatial_samples physical_k source_multiplier feedback");
    std::ifstream in(argv[1]);
    int nsteps,nbirth,angular;
    double weight,a_start,theta;
    in>>nsteps>>nbirth>>angular>>weight>>a_start>>theta;
    std::vector<Seed> births(nbirth);
    for(auto& birth:births)in>>birth.t>>birth.a>>birth.pz>>birth.pt>>birth.phi;
    const int n=std::stoi(argv[2]), nx=std::stoi(argv[3]);
    const double k=std::stod(argv[4])*a_start, amplitude=std::stod(argv[5]);
    const bool feedback=std::stoi(argv[6]);
    if(!in || n<1 || n>max_modes || nx<=2*n || angular<=0)throw std::runtime_error("Invalid grid");
    const double w=weight/nx;
    std::vector<Particle> particles;
    particles.reserve(nbirth*nx);
    Modes induced;
    double injected=0.,axial_work=0.,mass_work=0.,drive_work=0.,induced_work=0.,loss=0.;
    double peak_response=0.,peak_density=0.;
    int next=0;
    for(int i=0;i<nsteps;++i){
        Step d;
        in>>d.a0>>d.am>>d.a1>>d.deta>>d.x>>d.b>>d.births>>d.e_wave>>d.b_wave;
        if(!in)throw std::runtime_error("Truncated input");
        while(next<d.births){
            const auto& b=births.at(next);
            double shift=std::fmod((next/angular+1)*.61803398874989485,1.);
            for(int j=0;j<nx;++j){
                const double phase=2*pi*(j+shift)/nx;
                particles.push_back({phase,b.a*b.pt*std::cos(b.phi),b.a*b.pt*std::sin(b.phi),b.a*b.pz});
                injected+=w*b.a*std::hypot(1.,std::hypot(b.pt,b.pz));
            }
            ++next;
        }
        // The supplied E_y and B_z are a particular solution over the entire history.
        // Our Fourier unknowns have zero initial fields and evolve only in response to the daughters.
        const double half=.5*d.deta;
        const double meanE=d.am*d.am*d.x,guide=d.am*d.am*d.b;
        const double driveE=d.am*d.am*d.e_wave*amplitude,driveB=d.am*d.am*d.b_wave*amplitude;
        loss+=propagate(induced,k,half,d.am,n);
        for(auto& p:particles){
            const double g=norm(d.am,p.px,p.py,p.pz);
            mass_work+=w*(d.am-d.a0)*(d.am+d.a0)/(g+norm(d.a0,p.px,p.py,p.pz));
            p.x+=half*k*p.px/g;
            basis(p,n);
            p.by=dot(p.basis,induced.by,n);
            p.bz=dot(p.basis,induced.bz,n)+driveB*p.basis[1];
            gyro(p,d.am,p.by,guide+p.bz,half);
        }
        Field ey=induced.ey,ez=induced.ez,jy{},jz{};
        // Implicit midpoint in the wave-current pair with the exact discrete energy gradient velocity.
        if(feedback)for(int iter=0;iter<24;++iter){
            jy.fill(0.);jz.fill(0.);
            for(auto& p:particles){
                const double py=p.py+d.deta*Q*(driveE*p.basis[0]+dot(p.basis,ey,n));
                const double pz=p.pz+d.deta*Q*(meanE+dot(p.basis,ez,n));
                const double denom=norm(d.am,p.px,p.py,p.pz)+norm(d.am,p.px,py,pz);
                const double vy=(p.py+py)/denom,vz=(p.pz+pz)/denom;
                for(int j=0;j<2*n;++j){jy[j]+=p.basis[j]*vy;jz[j]+=p.basis[j]*vz;}
            }
            double diff=0.;
            for(int j=0;j<2*n;++j){
                const double y=induced.ey[j]-d.deta*Q*w/K*jy[j],z=induced.ez[j]-d.deta*Q*w/K*jz[j];
                diff=std::max({diff,std::abs(y-ey[j]),std::abs(z-ez[j])});ey[j]=y;ez[j]=z;
            }
            if(diff<2.e-13)break;
            if(iter==23)throw std::runtime_error("Electric iteration failed");
        }
        for(int j=0;j<2*n;++j){induced.ey[j]=2*ey[j]-induced.ey[j];induced.ez[j]=2*ez[j]-induced.ez[j];}
        for(auto& p:particles){
            const double gy=dot(p.basis,ey,n),gz=dot(p.basis,ez,n),ex=driveE*p.basis[0];
            const double g0=norm(d.am,p.px,p.py,p.pz),py=p.py,pz=p.pz;
            p.py+=d.deta*Q*(ex+gy); p.pz+=d.deta*Q*(meanE+gz);
            const double denom=g0+norm(d.am,p.px,p.py,p.pz);
            const double vy=(py+p.py)/denom,vz=(pz+p.pz)/denom;
            axial_work+=w*d.deta*Q*meanE*vz;
            drive_work+=w*d.deta*Q*ex*vy;
            induced_work+=w*d.deta*Q*(gy*vy+gz*vz);
            gyro(p,d.am,p.by,guide+p.bz,half);
            const double g=norm(d.am,p.px,p.py,p.pz);
            p.x=std::remainder(p.x+half*k*p.px/g,2*pi);
            mass_work+=w*(d.a1-d.am)*(d.a1+d.am)/(g+norm(d.a1,p.px,p.py,p.pz));
        }
        loss+=propagate(induced,k,half,d.am,n);
        const double response=field_energy(induced,n);
        peak_response=std::max(peak_response,std::sqrt(2.*response/K)/(d.a1*d.a1));
        if(!particles.empty()){
            double c=0.,s=0.;
            for(const auto& p:particles){c+=p.basis[0];s+=p.basis[1];}
            peak_density=std::max(peak_density,2.*std::hypot(c,s)/particles.size());
        }
        if(i+1==nsteps){
            double gas=0.,sqm=0.,axis=0.,sum_stress=0.,sx=0.,sy=0.,syz=0.,flow=0.;
            std::vector<double> pt;
            pt.reserve(particles.size());
            for(const auto& p:particles){
                const double g=norm(d.a1,p.px,p.py,p.pz),r=std::hypot(p.px,p.py)/d.a1;
                gas+=w*g; sqm+=r*r;axis+=p.pz*p.pz/g;
                sx+=p.px*p.px/g;sy+=p.py*p.py/g;syz+=p.py*p.pz/g;flow+=w*p.px;
                sum_stress+=(p.px*p.px+p.py*p.py+p.pz*p.pz)/g;
                pt.push_back(r);
            }
            std::sort(pt.begin(),pt.end());
            const double reference=std::cyl_bessel_k(2.,1./theta);
            double ks=0.,witness=0.;
            for(size_t j=0;j<pt.size();++j){
                const double wp=std::hypot(1.,pt[j]);
                const double predicted=1.-wp*wp*std::cyl_bessel_k(2.,wp/theta)/reference;
                const double gap=std::max(std::abs(predicted-double(j)/pt.size()),std::abs(predicted-double(j+1)/pt.size()));
                if(gap>ks){ks=gap;witness=pt[j];}
            }
            const double eig_split=std::hypot(sy-axis,2.*syz);
            std::array<double,3> eig={sx/sum_stress,.5*(sy+axis-eig_split)/sum_stress,.5*(sy+axis+eig_split)/sum_stress};
            std::sort(eig.begin(),eig.end());
            const double thermal=3.*theta+std::cyl_bessel_k(1.,1./theta)/reference;
            const double mean_gamma=gas/(w*pt.size()*d.a1);
            const double err=gas-injected-mass_work-axial_work-drive_work-induced_work;
            std::cout<<std::setprecision(16)<<particles.size()<<' '<<ks<<' '<<witness<<' '<<std::sqrt(sqm/pt.size())
                <<' '<<pt[pt.size()/2]<<' '<<axis/sum_stress<<' '<<sx/sum_stress<<' '<<sy/sum_stress
                <<' '<<syz/sum_stress<<' '<<gas<<' '<<injected<<' '<<mass_work<<' '<<axial_work<<' '<<drive_work
                <<' '<<induced_work<<' '<<response<<' '<<loss<<' '<<peak_response<<' '<<peak_density
                <<' '<<err/gas<<' '<<(response+loss+induced_work)/gas<<' '<<eig[0]<<' '<<eig[1]<<' '<<eig[2]<<' '<<flow/gas<<' '<<mean_gamma<<' '<<mean_gamma/thermal<<'\n';
        }
    }
    if(next!=nbirth)throw std::runtime_error("Incomplete injection");
}
