// Nonzero Fourier modes on a periodic, locally flat expanding wall patch.
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
struct Step {double a0,am,a1,deta,x,b;int births;};
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
int main(int argc,char** argv) {
    if(argc!=8)throw std::runtime_error("input modes spatial_samples physical_k seed_amplitude feedback polarization");
    std::ifstream file(argv[1]);int nstep,nbirth,angular;
    double weight,x0,tref;file>>nstep>>nbirth>>angular>>weight>>x0>>tref;
    std::vector<Seed> birth(nbirth);
    for(auto& b:birth)file>>b.t>>b.a>>b.pz>>b.pt>>b.phi;
    const int n=std::stoi(argv[2]),nx=std::stoi(argv[3]),pol=std::stoi(argv[7]);
    const double k=std::stod(argv[4])*x0,amp=std::stod(argv[5]);
    const bool feedback=std::stoi(argv[6]);
    if(n<1||n>max_modes||nx<=2*n||angular<=0)throw std::runtime_error("Invalid spectral grid");
    const double w=weight/nx;
    std::vector<Particle> parts;parts.reserve(nbirth*nx);
    Modes f;
    if(pol==0)f.by[0]=x0*x0*amp;
    else if(pol==1)f.ey[0]=x0*x0*amp;
    else if(pol==2)f.ez[0]=x0*x0*amp;
    else if(pol==3) {f.ez[0]=x0*x0*amp;f.by[0]=-f.ez[0];} // right-moving vacuum wave for analytic control
    else throw std::runtime_error("Invalid polarization");
    const double initial=field_energy(f,n);
    double injected=0.,external=0.,mass_work=0.,plasma_work=0.,max_defect=0.,peak_field=0.,peak_diff=0.,peak_contrast=0.;
    Modes vac=f;
    int next=0;
    for(int i=0;i<nstep;++i) {
        Step d; file>>d.a0>>d.am>>d.a1>>d.deta>>d.x>>d.b>>d.births;
        if(!file)throw std::runtime_error("Truncated trajectory");
        // The input splits a time step at each quantile birth. All new tags exist at its left edge.
        while(next<d.births) {
            const auto& b=birth[next];
            const int cohort=next/angular;
            for(int l=0;l<nx;++l){
                const double phase=2*pi*(l+std::fmod((cohort+1)*.61803398874989485,1.))/nx;
                parts.push_back({phase,b.a*b.pt*std::cos(b.phi),b.a*b.pt*std::sin(b.phi),b.a*b.pz});
                injected+=w*b.a*std::hypot(1.,std::hypot(b.pt,b.pz));
            }++next;
        }
        const double half=.5*d.deta;
        waves(f,k*half,n);waves(vac,k*d.deta,n);
        const double meanE=d.am*d.am*d.x,guide=d.am*d.am*d.b;
        for(auto& p:parts){
            double g=norm(d.am,p.px,p.py,p.pz);
            mass_work+=w*(d.am-d.a0)*(d.am+d.a0)/(g+norm(d.a0,p.px,p.py,p.pz));
            p.x+=half*k*p.px/g;basis(p,n);
            p.by=dot(p.basis,f.by,n);p.bz=dot(p.basis,f.bz,n);
            gyro(p,d.am,p.by,guide+p.bz,half);
        }
        const double before=field_energy(f,n);
        Field ey=f.ey,ez=f.ez,jy{},jz{};
        // Solve the simultaneous electric pushes and mode-current kick. The
        // ratio of sums gives the exact discrete derivative of particle energy.
        for(int it=0;it<24;++it){
            jy.fill(0);jz.fill(0);
            for(auto& p:parts){
                p.ey=dot(p.basis,ey,n);p.ez=dot(p.basis,ez,n);
                double py=p.py+d.deta*Q*p.ey,pz=p.pz+d.deta*Q*(meanE+p.ez);
                double denom=norm(d.am,p.px,p.py,p.pz)+norm(d.am,p.px,py,pz);
                double vy=(p.py+py)/denom,vz=(p.pz+pz)/denom;
                if(feedback)for(int j=0;j<2*n;++j){jy[j]+=p.basis[j]*vy;jz[j]+=p.basis[j]*vz;}
            }
            double diff=0.;
            for(int j=0;j<2*n;++j) {
                double y=f.ey[j]-d.deta*Q*w/K*jy[j],z=f.ez[j]-d.deta*Q*w/K*jz[j];
                diff=std::max({diff,std::abs(y-ey[j]),std::abs(z-ez[j])}); ey[j]=y;ez[j]=z;
            }
            if(diff<2e-13*(1+std::abs(amp*x0*x0))) break;
            if(it==23)throw std::runtime_error("Electric current iteration failed");
        }
        for(int j=0;j<2*n;++j){f.ey[j]=2*ey[j]-f.ey[j];f.ez[j]=2*ez[j]-f.ez[j];}
        double self_work=0.;
        for(auto& p:parts){
            p.ey=dot(p.basis,ey,n);p.ez=dot(p.basis,ez,n);
            double g0=norm(d.am,p.px,p.py,p.pz),py=p.py,pz=p.pz;
            p.py+=d.deta*Q*p.ey;p.pz+=d.deta*Q*(meanE+p.ez);
            double denom=g0+norm(d.am,p.px,p.py,p.pz),vy=(py+p.py)/denom,vz=(pz+p.pz)/denom;
            external+=w*d.deta*Q*meanE*vz;
            self_work+=w*d.deta*Q*(p.ey*vy+p.ez*vz);
            // Guide-aligned perturbations only turn transverse momentum in its plane.
            // |E_y|+|v_z*B_y| bounds the rate of increase of its magnitude.
            p.impulse+=d.deta*Q*(std::abs(p.ey)+std::abs(vz*p.by));
            gyro(p,d.am,p.by,guide+p.bz,half);
            double g=norm(d.am,p.px,p.py,p.pz);
            p.x=std::remainder(p.x+half*k*p.px/g,2*pi);
            mass_work+=w*(d.a1-d.am)*(d.a1+d.am)/(g+norm(d.a1,p.px,p.py,p.pz));
        }
        plasma_work+=self_work;
        if(feedback)max_defect=std::max(max_defect,std::abs(field_energy(f,n)-before+self_work));
        waves(f,k*half,n);
        peak_field=std::max(peak_field,std::sqrt(2*field_energy(f,n)/K)/(d.a1*d.a1));
        Modes diff=f;
        for(int j=0;j<2*n;++j) {
            diff.ey[j]-=vac.ey[j];diff.ez[j]-=vac.ez[j];diff.by[j]-=vac.by[j];diff.bz[j]-=vac.bz[j];
        }
        peak_diff=std::max(peak_diff,std::sqrt(2*field_energy(diff,n)/K)/(d.a1*d.a1));
        if(!parts.empty()){
            double c=0.,s=0.;
            for(const auto& p:parts){c+=p.basis[0];s+=p.basis[1];}
            peak_contrast=std::max(peak_contrast,2.*std::hypot(c,s)/parts.size());
        }
        if(i+1==nstep) {
            double gas=0.,p2=0.,tube=0.,j2=0.,zt=0.,tr=0.,p50=0.,p90=0.;
            std::vector<double> perps; perps.reserve(parts.size());
            for(const auto& p:parts) {
                const double g=norm(d.a1,p.px,p.py,p.pz),pt=std::hypot(p.px,p.py)/d.a1;
                gas+=w*g;p2+=pt*pt;tube+=(pt<=std::sqrt(.1025));j2+=p.impulse*p.impulse/(d.a1*d.a1);
                zt+=p.pz*p.pz/g;tr+=(p.px*p.px+p.py*p.py+p.pz*p.pz)/g;
                perps.push_back(pt);
            }
            std::sort(perps.begin(),perps.end()); p50=perps[perps.size()/2];p90=perps[(perps.size()*9)/10];
            double ks=0., reference=std::cyl_bessel_k(2.,1./tref);
            for(size_t j=0;j<perps.size();++j){
                double wperp=std::hypot(1.,perps[j]);
                double F=1.-wperp*wperp*std::cyl_bessel_k(2.,wperp/tref)/reference;
                ks=std::max({ks,std::abs(F-double(j)/perps.size()),std::abs(F-double(j+1)/perps.size())});
            }
            double resid=gas-injected-external-mass_work-plasma_work;
            std::cout<<std::setprecision(16)<<parts.size()<<' '<<gas<<' '<<field_energy(f,n)<<' '<<initial<<' '<<external<<' '<<plasma_work
                 <<' '<<resid/std::max(1e-16,gas)<<' '<<max_defect<<' '<<peak_field<<' '<<peak_diff<<' '<<peak_contrast<<' '<<std::sqrt(p2/parts.size())
                 <<' '<<p50<<' '<<p90<<' '<<std::sqrt(j2/parts.size())<<' '<<tube/parts.size()<<' '<<zt/tr<<' '<<mass_work<<' '<<injected<<' '<<ks<<'\n';
        }
    }
    if(next!=nbirth)throw std::runtime_error("Uninserted births");
}
