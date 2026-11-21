// Opposite-charge elastic vector-pole collisions, in units m=c=1.
// This header follows the #198 invariant-frame vector and particle definitions.
#include <cstdint>
#include <limits>
#include <random>

class PoleCollisions {
    std::mt19937_64 rng;
    std::vector<unsigned> contacts;
    const double density_factor, kernel_max=1.118569e-6;
    const bool isotropic;
    uint64_t proposals=0, events=0;
    double max_energy_error=0., max_momentum_error=0., max_acceptance=0.;
    double mu_residual=0., mu_variance=0.;

    double uniform() {return (static_cast<double>(rng() >> 12) + .5) * 0x1p-52;}

    void encounter(Particle& first, Particle& second, double radius, unsigned i, unsigned j) {
        const V p=(1./radius)*first.u, q=(1./radius)*second.u;
        const double ep=energy(1.,p), eq=energy(1.,q), lp=length(p), lq=length(q);
        const V P=p+q;
        // Rationalized scalar product remains accurate for nearly collinear streams.
        double minkowski=(1.+lp*lp+lq*lq)/(ep*eq+lp*lq);
        if (lp>0. && lq>0.) {
            V d=(1./lp)*p-(1./lq)*q;
            minkowski+=.5*lp*lq*dot(d,d);
        }
        const double s=std::max(4.,2.*(1.+minkowski));
        const double beta2=1.-4./s, beta=std::sqrt(beta2);
        const double coupling=.30*(s+2.), delta=s-8300.*8300.;
        const double kernel=6.*std::acos(-1.)*beta*coupling*coupling/
            (ep*eq*(delta*delta+beta2*coupling*coupling));
        const double accept=kernel/kernel_max;
        if (!(accept>=0. && accept<=1.)) throw std::runtime_error("Vector-pole envelope exceeded");
        max_acceptance=std::max(max_acceptance,accept);
        if (uniform()>=accept) return;
        ++events;
        ++contacts[i]; ++contacts[j];

        const double M=std::sqrt(s), total=ep+eq;
        const V incoming=p-((ep+.5*M)/(total+M))*P;
        const V axis=(1./length(incoming))*incoming;
        const V ref=std::abs(axis.z)<.8 ? V{0.,0.,1.} : V{1.,0.,0.};
        const V perp=cross(axis,ref);
        const V e1=(1./length(perp))*perp, e2=cross(axis,e1);
        // Contraction of the two spin-averaged vector-current tensors in the CM frame.
        const double A=isotropic?1.:3.-2.*beta2, B=isotropic?0.:beta2*beta2;
        double mu;
        do {mu=2.*uniform()-1.;} while(uniform()*(A+B)>A+B*mu*mu);
        const double phi=2.*std::acos(-1.)*uniform();
        const V side=std::cos(phi)*e1+std::sin(phi)*e2;
        const V star=(.5*M*beta)*(mu*axis+std::sqrt(1.-mu*mu)*side);
        const V next=star+(.5+dot(P,star)/(M*(total+M)))*P;
        const V other=P-next;
        max_energy_error=std::max(max_energy_error,std::abs(energy(1.,next)+energy(1.,other)-total)/total);
        max_momentum_error=std::max(max_momentum_error,length(next+other-P)/total);
        first.u=radius*next;
        second.u=radius*other;
        const double mean=(A/3.+B/5.)/(A+B/3.);
        mu_residual+=mu*mu-mean;
        mu_variance+=(A/5.+B/7.)/(A+B/3.)-mean*mean;
    }

public:
    PoleCollisions(size_t count, double factor, uint64_t seed, bool isotropic):
        rng(seed),contacts(count,0),density_factor(factor),isotropic(isotropic) {}

    void step(std::vector<Particle>& p, size_t n, double weight,
              double left, double right, double radius) {
        if (n==0 || density_factor==0.) return;
        if (n%2!=0) throw std::runtime_error("Unpaired births");
        const size_t pairs=n/2;
        // n*weight is total accepted comoving rest, including both signs and all
        // orientations in the symmetric density closure. This box represents one flavor.
        const double rate=pairs*density_factor*(n*weight)*kernel_max/(radius*radius*radius);
        double offset=0.;
        while ((offset-=std::log(uniform())/rate)<right-left) {
            const size_t i=2*static_cast<size_t>(uniform()*pairs);
            const size_t j=2*static_cast<size_t>(uniform()*pairs)+1;
            if (p[i].sign!=1. || p[j].sign!=-1.) throw std::runtime_error("Charge-pair order");
            if (p[i].birth>left+offset || p[j].birth>left+offset) continue;
            ++proposals;
            encounter(p[i],p[j],radius,i,j);
        }
    }

    uint64_t event_count() const { return events; }

    void save(const char* file) const {
        const size_t hit=std::count_if(contacts.begin(),contacts.end(),[](unsigned v){return v>0;});
        const double fraction=static_cast<double>(hit)/contacts.size();
        std::ofstream out(file);
        out << std::setprecision(17) << proposals << " " << events << " " << fraction << " "
            << max_energy_error << " " << max_momentum_error << " " << max_acceptance << " "
            << mu_residual << " " << mu_variance << "\n";
        if (!out) throw std::runtime_error("Missing collision diagnostics");
    }
};
