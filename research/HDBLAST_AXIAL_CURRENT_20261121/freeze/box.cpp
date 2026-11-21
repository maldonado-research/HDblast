#include "vector.hpp"
#include "pole.hpp"

double PoleCollisions::attempt(std::vector<Particle>& particles) {
    const size_t pairs=particles.size()/2;
    const double lapse=-std::log(uniform())/pairs;
    const size_t i=2*static_cast<size_t>(uniform()*pairs);
    const size_t j=2*static_cast<size_t>(uniform()*pairs)+1;
    ++proposals;
    encounter(particles[i],particles[j],1.,i,j);
    return lapse;
}

int main(int argc,char** argv) {
    if(argc!=7) throw std::runtime_error("data, output, seed, flat, max-count and checks required");
    const auto raw=read<4>(argv[1]);
    std::vector<Particle> particles;
    for (const auto& v:raw) particles.push_back({V{v[1],v[2],v[3]},v[0],0.,1.});
    const size_t total=particles.size(), pairs=total/2;
    if (total<2 || total%2!=0) throw std::runtime_error("Unpaired population");
    for(size_t i=0;i<total;i++) if (particles[i].sign!=(i%2?-1.:1.)) throw std::runtime_error("Wrong sign");
    PoleCollisions pole(total,1.,std::stoull(argv[3]),false);
    pole.flat=std::stoi(argv[4])!=0;
    std::ofstream output(argv[2],std::ios::binary);
    double clock=0;
    const size_t end=std::stoul(argv[5]);
    for(size_t count: {0,1,2,4,8,16,32}) {
        if (count>end) break;
        while(pole.event_count()<pairs*count) clock+=pole.attempt(particles);
        output.write(reinterpret_cast<const char*>(&clock),sizeof(double));
        for(size_t i=0;i<total;++i) {
            const V u=particles[i].u;
            const std::array<double,4> out={u.x,u.y,u.z,double(pole.hits(i))};
            output.write(reinterpret_cast<const char*>(out.data()),sizeof(out));
        }
    }
    pole.save(argv[6]);
    if(!output) throw std::runtime_error("Writing snapshots failed");
}
