// Independent read-only Dirac pair annihilation estimator. Momenta and energies in units of m.
#include <cstdint>
#include <random>

double spin_half_kernel(V p,V q,double ep,double eq) {
    const double r=length(p), s=length(q);
    double d=(r-s)*(r-s)/(ep*eq+r*s+1.);
    if (r>0. && s>0.) {
        V delta=(1./r)*p-(1./s)*q;
        d+=.5*r*s*dot(delta,delta);
    }
    const double root=std::sqrt(d*(d+2.));
    const double t=d==0. ? 1. : std::log1p(d+root)/root;
    const double value=((d*(d+6.)+6.)*t-(d+4.))/((d+2.)*ep*eq);
    if (!(value>=0. && value<=2.125+1.e-12)) throw std::runtime_error("Annihilation invariant bound");
    return value;  // m^2 sigma v_M divided by pi alpha^2
}

void check_annihilation_kernel(){
    auto check=[](double value,double target){
        if(std::abs(value-target)>3e-10*target)throw std::runtime_error("Annihilation cross-section control");
    };
    check(spin_half_kernel({},{},1.,1.),1.);
    for(double u:{.01,1.,30.,1e4}){
        V p={0.,0.,u};
        const double e=energy(1.,p),beta=u/e;
        check(spin_half_kernel(p,p,e,e),1./(e*e));
        // Independent CM expression via detailed balance from Breit-Wheeler photon pair creation.
        const double L=2.*std::asinh(u);
        const double factor=((3.-std::pow(beta,4))*L-2.*beta*(2.-beta*beta))/(2.*beta*e*e);
        check(spin_half_kernel(p,(-1.)*p,e,e),factor);
    }
}

class Annihilation {
    static constexpr int batches=4;
    static constexpr int strata=256;
    std::array<std::mt19937_64,batches> seeds;
    std::array<std::array<double,2>,batches> sums{};
    double count_cap=0.,energy_cap=0.;
    std::vector<V> before;
    std::vector<double> positive,negative;
    static double unit(std::mt19937_64& rng) {return double(rng()>>11)*0x1p-53;}
public:
    Annihilation(size_t n):before(n),positive(n/2),negative(n/2) {
        check_annihilation_kernel();
        for(int i=0;i<batches;++i)seeds[i].seed(4799921+124817*i);
    }
    void begin(const std::vector<Particle>& p,size_t live) {
        for(size_t i=0;i<live;++i)before[i]=p[i].u;
    }
    void add(const std::vector<Particle>& p,size_t live,double w,double left,double right,double al,double ar) {
        if(live==0)return;
        const size_t pairs=live/2;
        if(live!=2*pairs)throw std::runtime_error("Incomplete opposite-charge pair");
        const double dt=right-left,D=w*live,amin=std::min(al,ar),amax=std::max(al,ar);
        const double pref=D*D/strata;
        double maximum_conformal=0.;
        double sums_positive=0.,sums_negative=0.;
        for(size_t i=0;i<live;++i){
            maximum_conformal+=std::max(energy(amax,before[i]),energy(amax,p[i].u));
            // Half uniform and half inverse midpoint gamma resolves rare cold tags.
            const double inverse=1./energy(1.,(1./(al+ar))*(before[i]+p[i].u));
            if(i%2==0)positive[i/2]=(sums_positive+=inverse);
            else negative[i/2]=(sums_negative+=inverse);
        }
        for(size_t k=0;k<pairs;++k){
            const double uniform_mass=double(k+1)/pairs;
            positive[k]=.5*(uniform_mass+positive[k]/sums_positive);
            negative[k]=.5*(uniform_mass+negative[k]/sums_negative);
        }
        positive[pairs-1]=negative[pairs-1]=1.;
        count_cap+=dt*2.125*D*D/(amin*amin*amin);
        energy_cap+=dt*2.125*D*w*maximum_conformal/(amin*amin*amin);
        for(int b=0;b<batches;++b){
            auto& rng=seeds[b];
            for(int k=0;k<strata;++k){
                const double u=std::min((k+unit(rng))/strata,std::nextafter(1.,0.));
                const size_t pindex=std::upper_bound(positive.begin(),positive.begin()+pairs,u)-positive.begin();
                const size_t nindex=std::upper_bound(negative.begin(),negative.begin()+pairs,unit(rng))-negative.begin();
                const size_t i=2*pindex,j=2*nindex+1;
                const double pprob=positive[pindex]-(pindex ? positive[pindex-1] : 0.);
                const double nprob=negative[nindex]-(nindex ? negative[nindex-1] : 0.);
                const double correction=1./(double(pairs)*pprob*double(pairs)*nprob);
                if(p[i].sign!=1. || p[j].sign!=-1.)throw std::runtime_error("Unexpected charge order");
                const double start=std::max(left,std::max(p[i].birth,p[j].birth));
                if(start>=right)continue;
                // Independent uniform quadrature time on the overlap of the two lifetimes.
                const double t=start+(right-start)*unit(rng), f=(t-left)/dt;
                const double a=(1-f)*al+f*ar;
                auto interpolate=[&](size_t idx){
                    const double t0=std::max(left,p[idx].birth);
                    const double frac=(t-t0)/(right-t0);
                    return (1./a)*((1-frac)*before[idx]+frac*p[idx].u);
                };
                V pi=interpolate(i),pj=interpolate(j);
                const double ei=energy(1.,pi),ej=energy(1.,pj);
                const double net=correction*pref*(right-start)*spin_half_kernel(pi,pj,ei,ej)/(a*a*a);
                sums[b][0]+=net;
                sums[b][1]+=.5*a*(ei+ej)*net;
            }
        }
    }
    void save(const char* filename)const {
        std::ofstream out(filename);out<<std::setprecision(17);
        for(const auto& row:sums)out<<row[0]<<" "<<row[1]<<"\n";
        out<<count_cap<<" "<<energy_cap<<"\n";
        if(!out)throw std::runtime_error("Annihilation output");
    }
};
