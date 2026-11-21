// All momenta use one flavor's local axes. Its two physical replicas are cyclic rotations.
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <random>

class FlavorJump {
    std::mt19937_64 stream;
    bool convert;
    std::array<uint64_t,3> destinations{};
    double p_residual=0., energy_residual=0.;

public:
    FlavorJump(uint64_t seed, bool enabled) : stream(seed ^ 0x4ec8d231ab5910f3ULL), convert(enabled) {}
    static V to_global(V p, unsigned f) {
        if(f==1)return {p.z,p.x,p.y};
        if(f==2)return {p.y,p.z,p.x};
        return p;
    }

    void project(V& plus, V& minus, V incident_p, double radius, double e_pair) {
        // Reject the single unused integer so the three choices have exactly equal probability.
        uint64_t bits;
        do {bits=stream();} while(bits==UINT64_MAX);
        const unsigned shift=convert?unsigned(bits%3):0;
        ++destinations[shift];
        plus=to_global(plus,(3-shift)%3);
        minus=to_global(minus,(3-shift)%3);
        const V outgoing=(1./radius)*to_global(plus+minus,shift);
        p_residual=std::max(p_residual,length(outgoing-incident_p)/e_pair);
        energy_residual=std::max(energy_residual,std::abs((energy(radius,plus)+energy(radius,minus))/radius-e_pair)/e_pair);
    }

    std::array<uint64_t,3> counts() const {return destinations;}

    void save(const char* fn) const {
        std::ofstream out(fn);
        out<<std::setprecision(17)<<destinations[0]<<" "<<destinations[1]<<" "<<destinations[2]<<" "
           <<p_residual<<" "<<energy_residual<<"\n";
        if(!out)throw std::runtime_error("Flavor transfer report");
    }
};
