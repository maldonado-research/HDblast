// Test absolute collision normalization against an exactly solvable two-particle gas.
#define main simulated_shell_main
#include "evolve.cpp"
#undef main

int main() {
    const size_t trials=8192;
    const double mass=8300., beta=std::sqrt(1.-4./(mass*mass));
    const double k_res=24.*std::acos(-1.)/(mass*mass*beta);
    const V incoming={0.,0.,std::sqrt(.25*mass*mass-1.)};
    for (int mode=0;mode<3;++mode) {
        const double birth=mode==1 ? .5 : 0.;
        const V other=mode==2 ? incoming : (-1.)*incoming;
        unsigned hit=0;
        uint64_t events=0;
        for (size_t i=0;i<trials;++i) {
            std::vector<Particle> p={{incoming,1.,birth,1.},{other,-1.,birth,1.}};
            PoleCollisions pole(2,.8/k_res,2026112100ULL+i,false);
            pole.step(p,2,.5,0.,1.,1.);
            events+=pole.event_count();
            hit+=pole.event_count()>0;
        }
        std::cout<<mode<<" "<<trials<<" "<<events<<" "<<hit<<"\n";
    }
}
