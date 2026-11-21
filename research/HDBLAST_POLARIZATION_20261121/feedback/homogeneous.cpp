// Two unit-RMS horizontal Hopf fields on the round wall. Momentum components are comoving.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

constexpr double K=2.54, Q=1.e7, kappa=2.;
struct V {double x=0.,y=0.,z=0.;};
V operator+(V a,V b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
V operator-(V a,V b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
V operator*(double s,V a){return {s*a.x,s*a.y,s*a.z};}
double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double length(V a){return std::sqrt(dot(a,a));}
V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double energy(double a,V u){return std::sqrt(a*a+dot(u,u));}
V lateral(V v){return {v.x,v.y,0.};}
template<int n> std::vector<std::array<double,n>> read(const char* name){
    std::ifstream in(name,std::ios::binary);
    if(!in)throw std::runtime_error("open input");
    std::vector<std::array<double,n>> rows;std::array<double,n> r;
    while(in.read(reinterpret_cast<char*>(r.data()),sizeof(r)))rows.push_back(r);
    if(!in.eof()||in.gcount()||rows.empty())throw std::runtime_error("partial input");
    return rows;
}

struct Particle{V u;double sign,birth,birth_a,factor=0.,g0=0.;};
struct Wave{V e,b;};
double propagate(Wave& w,double h,double a){
    // Symmetric electric leakage and the exact curl rotation for curl=-2.
    double s=std::sin(kappa*h),c=std::cos(kappa*h);
    double attenuation=std::exp(-a*h/K);
    double factor=.5*K*(-std::expm1(-2*a*h/K));
    double loss=factor*dot(w.e,w.e);
    V e=attenuation*w.e,b=w.b;
    w.e=c*e-s*b;w.b=c*b+s*e;
    loss+=factor*dot(w.e,w.e);
    w.e=attenuation*w.e;
    return loss;
}
V rotate(V u,V magnetic,double a,double q_duration){
    V t=(.5*q_duration/energy(a,u))*magnetic;
    V intermediate=u+cross(u,t);
    return u+cross(intermediate,(2./(1.+dot(t,t)))*t);
}
int main(int argc,char** argv){
    if(argc!=7)throw std::runtime_error("fields seed net_mass feedback raw history");
    const auto f=read<8>(argv[1]),births=read<8>(argv[2]);
    const double weight=std::stod(argv[3])/(2*births.size());const bool active=std::stoi(argv[4])!=0;
    std::ofstream raw(argv[5],std::ios::binary),hist(argv[6],std::ios::binary);
    if(!raw||!hist||f.size()<2)throw std::runtime_error("invalid output or field table");
    std::vector<Particle> p;
    for(const auto& b:births){
        for(int sign:{1,-1})p.push_back({sign*V{b[5],b[6],b[7]},double(sign),b[0],b[1]});
    }
    std::stable_sort(p.begin(),p.end(),[](const Particle& a,const Particle& b){return a.birth<b.birth;});
    Wave wave;size_t n=0;
    double loss=0.,birth_energy=0.,expansion=0.,axial=0.,drive=0.,induced=0.,bound=0.;
    int max_iterations=0;
    for(size_t j=0;j+1<f.size();++j){
        const auto& l=f[j];const auto& r=f[j+1];double dt=r[0]-l[0];
        auto mid=[&](int k){return .5*(l[k]+r[k]);};
        double a=mid(1),a2=a*a,deta=dt/a;
        V external=a2*V{mid(4),mid(6),mid(2)},external_b=a2*V{mid(5),mid(7),mid(3)};
        loss+=propagate(wave,.5*deta,a);
        const size_t old=n;
        while(n<p.size()&&p[n].birth<=r[0]){
            birth_energy+=weight*energy(p[n].birth_a,p[n].u);++n;
        }
        double extent=0.;
        for(size_t i=0;i<n;++i){
            auto& v=p[i];const double start_a=i<old?l[1]:v.birth_a;
            double h=(r[0]-std::max(v.birth,l[0]))/a;
            extent+=h;
            v.factor=Q*v.sign*h;
            double g=energy(a,v.u);
            expansion+=weight*(a-start_a)*(a+start_a)/(g+energy(start_a,v.u));
            v.u=rotate(v.u,external_b+wave.b,a,.5*v.factor);
            v.g0=energy(a,v.u);
        }
        bound+=Q*weight*extent/K;
        V avg=wave.e,current;
        if(active){
            for(int iter=1;iter<=28;++iter){
                current={};
                for(size_t i=0;i<n;++i){
                    const auto& v=p[i];V next=v.u+v.factor*(external+avg);
                    double denom=v.g0+energy(a,next);
                    current=current+(weight*v.factor/denom)*lateral(v.u+next);
                }
                V updated=wave.e-(.5/K)*current;
                const double residual=length(updated-avg);
                avg=updated;max_iterations=std::max(max_iterations,iter);
                if(residual<3.e-15*(1.+length(avg)))break;
                if(iter==28)throw std::runtime_error("implicit current failed");
            }
        }
        wave.e=2.*avg-wave.e;
        V used_current;
        for(size_t i=0;i<n;++i){
            auto& v=p[i];V next=v.u+v.factor*(external+avg);
            V flux=(weight*v.factor/(v.g0+energy(a,next)))*(v.u+next);
            used_current=used_current+lateral(flux);
            axial+=flux.z*external.z;drive+=dot(flux,lateral(external));induced+=dot(flux,avg);
            v.u=rotate(next,external_b+wave.b,a,.5*v.factor);
            expansion+=weight*(r[1]-a)*(r[1]+a)/(energy(a,v.u)+energy(r[1],v.u));
        }
        loss+=propagate(wave,.5*deta,a);
        const std::array<double,10> trace={r[0],wave.e.x,wave.e.y,wave.b.x,wave.b.y,
                    used_current.x/deta,used_current.y/deta,bound,double(n)*weight,loss};
        hist.write(reinterpret_cast<const char*>(trace.data()),sizeof(trace));
    }
    if(n!=p.size())throw std::runtime_error("incomplete births");
    const double a=f.back()[1];double gas=0.;
    for(const auto& v:p){
        gas+=weight*energy(a,v.u);
        std::array<double,8> r={v.sign,1.,0.,0.,0.,v.u.x/a,v.u.y/a,v.u.z/a};
        raw.write(reinterpret_cast<const char*>(r.data()),sizeof(r));
    }
    double field=.5*K*(dot(wave.e,wave.e)+dot(wave.b,wave.b));
    std::cout<<std::setprecision(17)<<gas<<" "<<birth_energy<<" "<<expansion<<" "<<axial<<" "<<drive<<" "
         <<induced<<" "<<field<<" "<<loss<<" "<<max_iterations<<"\n";
    if(!raw||!hist)throw std::runtime_error("output write");
}
