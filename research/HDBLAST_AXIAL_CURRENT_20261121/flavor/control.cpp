// A two-particle resonant gas sees a constant inclusive rate even as both labels change.
#define main shared_shell_main
#include "shared.cpp"
#undef main

int main() {
    constexpr unsigned trials=8192;
    const double mass=8300., beta=std::sqrt(1.-4./(mass*mass));
    const double peak=24.*std::acos(-1.)/(mass*mass*beta);
    const V p={0.,0.,std::sqrt(.25*mass*mass-1.)};
    unsigned events=0,hits=0;
    std::array<unsigned,3> destinations{};
    for(unsigned i=0;i<trials;++i) {
        std::vector<Particle> pair={{p,1.,0.,1.},{(-1.)*p,-1.,0.,1.}};
        PoleCollisions pole(2,2.4/(3.*peak),430071+i,false,true);
        pole.step(pair,2,.5,0.,1.,1.);
        events+=pole.event_count();
        hits+=pole.event_count()>0;
        auto counts=pole.flavor_counts();
        for(unsigned f=0;f<3;++f)destinations[f]+=counts[f];
    }
    // Compare with the explicit physical bases in three Cartesian boxes. The
    // outgoing pair has nonzero total momentum so a wrong inverse is detectable.
    const std::array<std::array<unsigned,3>,3> coordinates={{{0,1,2},{2,0,1},{1,2,0}}};
    auto rotate_explicit=[&](V q,unsigned f){
        const std::array<double,3> c={q.x,q.y,q.z};
        return V{c[coordinates[f][0]],c[coordinates[f][1]],c[coordinates[f][2]]};
    };
    double error=0.,charge_error=0.;
    for(unsigned seed=0;seed<300;++seed) {
        const V in_plus={.3,-2.4,9.},in_minus={-1.1,3.,-.25};
        V plus=in_plus,minus=in_minus;
        FlavorJump jumps(seed,true);
        jumps.project(plus,minus,in_plus+in_minus,1.,energy(1.,plus)+energy(1.,minus));
        const auto count=jumps.counts();
        const auto shift=std::find(count.begin(),count.end(),1)-count.begin();
        if(shift>=3)throw std::runtime_error("Missing destination");
        for(unsigned f=0;f<3;++f){
            unsigned g=(f+shift)%3;
            error=std::max({error,length(rotate_explicit(plus,g)-rotate_explicit(in_plus,f)),
                            length(rotate_explicit(minus,g)-rotate_explicit(in_minus,f))});
            std::array<int,3> plus_change{},minus_change{};
            plus_change[f]--;plus_change[g]++;minus_change[f]--;minus_change[g]++;
            for(unsigned i=0;i<3;++i)charge_error=std::max(charge_error,double(std::abs(plus_change[i]-minus_change[i])));
        }
    }
    double chi=0.; for(unsigned count:destinations)chi+=std::pow(count-events/3.,2)/(events/3.);
    std::cout<<std::setprecision(17)<<trials<<" "<<events<<" "<<hits<<" "<<destinations[0]<<" "
             <<destinations[1]<<" "<<destinations[2]<<" "<<error<<" "<<charge_error<<" "<<chi<<"\n";
}
