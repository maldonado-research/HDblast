// Homogeneous relativistic Boltzmann collision box, one distinguishable sign pair, m=1.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

constexpr double pole=8300., width=.30, pi=3.14159265358979323846;
struct Vec { double x=0., y=0., z=0.; };
Vec operator+(Vec a, Vec b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vec operator-(Vec a, Vec b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vec operator*(Vec a,double b) { return {a.x*b,a.y*b,a.z*b}; }
double dot(Vec a, Vec b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec cross(Vec a,Vec b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
double length(Vec v) { return std::sqrt(dot(v,v)); }
struct Particle { Vec p; double e; int hits=0; };

struct Random {
    std::mt19937_64 engine;
    explicit Random(uint64_t seed): engine(seed) {}
    double unit() { return ((engine()>>11)+.5)*0x1.0p-53; }
    int pick(int n) { return std::min(n-1, int(n*unit())); }
};

// Energy-weighted tree used only when the cold off-shell majorant is more efficient.
struct Tree {
    std::vector<double> nodes;
    double total=0.;
    explicit Tree(int n): nodes(n+1,0.) {}
    void add(int at,double amount) {
        total+=amount;
        for(int i=at+1;i<int(nodes.size());i+=i&-i) nodes[i]+=amount;
    }
    int pick(double u) const {
        double target=u*total;
        int j=0,step=1;
        while(2*step<int(nodes.size())) step*=2;
        for(;step;step/=2) {
            int k=j+step;
            if(k<int(nodes.size()) && target>=nodes[k]) { target-=nodes[k]; j=k; }
        }
        return std::min(j,int(nodes.size())-2);
    }
};

struct Kernel { double s, flux_xs; };
Kernel kernel(const Particle& a,const Particle& b) {
    // p1.p2 in mostly-minus metric. Double precision with a floor at the exact threshold.
    const double s=std::max(4.,2.+2.*(a.e*b.e-dot(a.p,b.p)));
    const double beta=std::sqrt(1.-4./s), imag=width*(s+2.)*beta;
    const double sig=12.*pi*std::pow(width*(s+2.),2)/(s*(std::pow(s-pole*pole,2)+imag*imag));
    return {s, .5*s*beta*sig/(a.e*b.e)};
}

struct Bath {
    int n; double energy, p0, sigma0, initial_p, weighted_bound, constant_bound;
    bool weighted;
    std::vector<Particle> pos,neg;
    Tree plus_energy,minus_energy;
    Random rng,probe;
    uint64_t hits=0, attempts=0;
    double angular_residual=0., worst_event_energy=0.;
    Bath(int count,double e,std::string shape,uint64_t seed)
        : n(count),energy(e),p0(std::sqrt((e-1.)*(e+1.))),plus_energy(n),minus_energy(n),
          rng(seed),probe(seed^0x51300dafeedULL) {
        Particle f{{0,0,p0},e},b{{0,0,-p0},e};
        sigma0=kernel(f,b).flux_xs;
        initial_p=0.;
        if (n%2 || !(e>1.) || (shape!="mono" && shape!="split" && shape!="thermal")) throw std::runtime_error("Bad initial ensemble");
        for(int i=0;i<n;++i) {
            const double ei=shape!="split"?e:1.+(i<n/2?.15:1.85)*(e-1.);
            const double z=std::sqrt((ei-1.)*(ei+1.));
            pos.push_back({{0,0,z},ei}); neg.push_back({{0,0,-z},ei});
            plus_energy.add(i,ei); minus_energy.add(i,ei); initial_p+=z/n;
        }
        const double m2=pole*pole, q=width*width*(1.-400./m2);
        weighted_bound=96.*pi/(m2*m2)*std::max(2.25*width*width/(.99*.99),
                                     width*width*(1.+q)/q*std::pow(1.+200./m2,2));
        constant_bound=std::max(108.*pi*width*width/m2,48.*pi/(m2-8.));
        weighted=weighted_bound*e*e<constant_bound;
    }

    void scatter(int i,int j,double s) {
        auto &a=pos[i],&b=neg[j];
        const double energy_sum=a.e+b.e, w=std::sqrt(s), k=std::sqrt(.25*s-1.);
        const Vec sum=a.p+b.p;
        Vec incoming=a.p-sum*((a.e+.5*w)/(energy_sum+w));
        Vec u=incoming*(1./length(incoming));
        Vec side=cross(u, std::abs(u.z)<.8?Vec{0,0,1}:Vec{1,0,0});
        side=side*(1./length(side));
        const double beta2=(s-4.)/s, A=3.-2.*beta2, B=beta2*beta2;
        double c;
        do { c=2.*rng.unit()-1.; } while (rng.unit()*(A+B)>A+B*c*c);
        angular_residual+=c*c-(A/3.+B/5.)/(A+B/3.);
        const double phi=2.*pi*rng.unit(), sn=std::sqrt(1.-c*c);
        const Vec q=(u*c + side*(sn*std::cos(phi)) + cross(u,side)*(sn*std::sin(phi)))*k;
        const Vec p=q+sum*(.5+dot(q,sum)/(w*(energy_sum+w)));
        const double ea=std::sqrt(1.+dot(p,p)), eb=std::sqrt(1.+dot(sum-p,sum-p));
        worst_event_energy=std::max(worst_event_energy,std::abs((ea+eb)/energy_sum-1.));
        if(weighted) { plus_energy.add(i,ea-a.e); minus_energy.add(j,eb-b.e); }
        a={p,ea,a.hits+1}; b={sum-p,eb,b.hits+1};
        ++hits;
    }

    double proposal_rate() const {
        return weighted?weighted_bound*plus_energy.total*minus_energy.total/(n*sigma0):
                        n*constant_bound/sigma0;
    }
    void propose() {
        const int i=weighted?plus_energy.pick(rng.unit()):rng.pick(n);
        const int j=weighted?minus_energy.pick(rng.unit()):rng.pick(n);
        const auto K=kernel(pos[i],neg[j]);
        const double bound=weighted?weighted_bound*pos[i].e*neg[j].e:constant_bound;
        ++attempts;
        if(K.flux_xs>bound*(1.+1e-12)) throw std::runtime_error("Collision majorant exceeded");
        if(rng.unit()*bound<K.flux_xs) scatter(i,j,K.s);
    }
    void batch(double dt) {
        std::vector<int> idx(n);
        std::iota(idx.begin(),idx.end(),0);
        for(int j=n-1;j>0;--j) std::swap(idx[j],idx[rng.pick(j+1)]);
        for(int i=0;i<n;++i) {
            const auto K=kernel(pos[i],neg[idx[i]]);
            if(rng.unit() < -std::expm1(-dt*K.flux_xs/sigma0)) scatter(i,idx[i],K.s);
        }
    }
};

struct Maxwell {
    double theta,xmax; std::vector<double> cdf;
    double prob(double e) const {
        const double u=std::clamp((e-1.)/theta,0.,xmax)*(cdf.size()-1)/xmax;
        const int j=std::min(int(cdf.size())-2,int(u));
        return cdf[j]+(u-j)*(cdf[j+1]-cdf[j]);
    }
    double ks(std::vector<double> e) const {
        std::sort(e.begin(),e.end());
        double d=0.; const double n=e.size();
        for(size_t i=0;i<e.size();++i) {
            const double p=prob(e[i]);
            d=std::max({d,std::abs(i/n-p),std::abs((i+1)/n-p)});
        }
        return d;
    }
};

// Draw the equilibrium reference independently, then remove its tiny sample COM
// velocity and rescale all spatial momenta to the exact target energy.
void thermal_start(Bath& b,const Maxwell& mj) {
    double total_e=0.; Vec total_p{};
    for(auto* side:{&b.pos,&b.neg}) for(auto& a:*side) {
        const double u=b.rng.unit();
        int j=std::clamp(int(std::lower_bound(mj.cdf.begin(),mj.cdf.end(),u)-mj.cdf.begin()),1,int(mj.cdf.size())-1);
        const double f=(u-mj.cdf[j-1])/(mj.cdf[j]-mj.cdf[j-1]);
        const double x=mj.xmax*(j-1+f)/(mj.cdf.size()-1);
        const double e=1.+mj.theta*x, p=std::sqrt((e-1.)*(e+1.)), z=2.*b.rng.unit()-1.;
        const double phi=2.*pi*b.rng.unit(), pt=p*std::sqrt(1-z*z);
        a={{pt*std::cos(phi),pt*std::sin(phi),p*z},e}; total_p=total_p+a.p; total_e+=e;
    }
    const double w=std::sqrt(total_e*total_e-dot(total_p,total_p));
    for(auto* side:{&b.pos,&b.neg}) for(auto& a:*side) {
        a.p=a.p+total_p*((dot(a.p,total_p)/(total_e+w)-a.e)/w);
    }
    double lo=.5,hi=2.;
    for(int k=0;k<53;++k) {
        double mid=.5*(lo+hi), avg=0.;
        for(auto* side:{&b.pos,&b.neg}) for(const auto& a:*side) avg+=std::sqrt(1+mid*mid*dot(a.p,a.p))/(2.*b.n);
        if(avg>b.energy) hi=mid; else lo=mid;
    }
    b.plus_energy=Tree(b.n); b.minus_energy=Tree(b.n);
    for(int j=0;j<2;++j) for(int i=0;i<b.n;++i) {
        auto& a=(j?b.neg:b.pos)[i]; a.p=a.p*(.5*(lo+hi)); a.e=std::sqrt(1.+dot(a.p,a.p));
        (j?b.minus_energy:b.plus_energy).add(i,a.e);
    }
    b.initial_p=b.p0;
}

void snapshot(Bath& bath,const Maxwell& mj,double tau) {
    std::array<std::vector<double>,2> energies;
    double z=0., stress=0., tr=0., l4=0., e_sum=0., e2=0., untouched=0.;
    Vec psum{};
    for(int side=0;side<2;++side) {
        const auto& particles=side==0?bath.pos:bath.neg;
        energies[side].reserve(bath.n);
        for(const auto &a:particles) {
            const double p2=dot(a.p,a.p);
            z+=(side==0?1.:-1.)*a.p.z;
            stress+=(1.5*a.p.z*a.p.z-.5*p2)/a.e;
            tr+=p2/a.e; psum=psum+a.p;
            if(p2>0.) {const double u2=a.p.z*a.p.z/p2;l4+=(35*u2*u2-30*u2+3.)/8.;}
            e_sum+=a.e; e2+=a.e*a.e;
            untouched+=(a.hits==0);
            energies[side].push_back(a.e);
        }
    }
    const double ks_plus=mj.ks(energies[0]),ks_minus=mj.ks(energies[1]);
    energies[0].insert(energies[0].end(),energies[1].begin(),energies[1].end());
    const int Q=32768;
    double k=0.,drag=0.;
    for(int l=0;l<Q;++l) {
        const auto &a=bath.pos[bath.probe.pick(bath.n)], &b=bath.neg[bath.probe.pick(bath.n)];
        const double r=kernel(a,b).flux_xs/bath.sigma0;
        k+=r;drag+=r*(a.p.z-b.p.z)/(2.*bath.initial_p);
    }
    const double N=2.*bath.n;
    std::cout<<tau<<' '<<z/(N*bath.initial_p)<<' '<<stress/tr<<' '<<l4/N<<' '
             <<mj.ks(energies[0])<<' '<<std::max(ks_plus,ks_minus)<<' '
             <<e2/N/(bath.energy*bath.energy)<<' '<<untouched/N<<' '
             <<double(bath.hits)/bath.n<<' '<<k/Q<<' '<<drag/Q<<' '
             <<std::abs(e_sum/N/bath.energy-1.)<<' '<<length(psum)/(N*bath.p0)<<' '
             <<bath.worst_event_energy<<' '<<bath.angular_residual/std::max<uint64_t>(1,bath.hits)<<'\n';
}

int main(int argc,char** argv) {
    if(argc!=7) throw std::runtime_error("COUNT ENERGY SHAPE SEED TIMESTEP CONFIG");
    Bath bath(std::stoi(argv[1]),std::stod(argv[2]),argv[3],std::stoull(argv[4]));
    const double dt=std::stod(argv[5]);
    std::ifstream file(argv[6]);
    int n,m;
    Maxwell mj;file>>mj.theta>>mj.xmax>>n>>m;
    mj.cdf.resize(n);std::vector<double> times(m);
    for(double& x:mj.cdf) file>>x;
    for(double& x:times) file>>x;
    if(!file||times[0]!=0.||!std::is_sorted(times.begin(),times.end())) throw std::runtime_error("Bad snapshot config");
    if(std::string(argv[3])=="thermal") thermal_start(bath,mj);
    std::cout<<std::setprecision(12);
    double t=0.,next=-std::log(bath.rng.unit())/bath.proposal_rate();
    for(double target:times) {
        if(dt==0.) {
            while(next<=target) {
                t=next; bath.propose();
                next=t-std::log(bath.rng.unit())/bath.proposal_rate();
            }
        } else {
            while(t+1e-11<target) {
                const double step=std::min(target-t,dt);
                bath.batch(step); t+=step;
            }
        }
        snapshot(bath,mj,target);
    }
}
