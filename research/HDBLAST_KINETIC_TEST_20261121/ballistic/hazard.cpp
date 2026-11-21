#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

struct Birth {
    double t, axial, transverse;
};
struct Target {
    double z, q, E;
};
struct State {
    double t, dt, a, impulse, density; // density = n_lab L^3 / (mL)^2
};
template<class T> void read_array(std::ifstream& input, std::vector<T>& v) {
    input.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size()*sizeof(T)));
}
int main(int argc, char** argv) {
    if (argc!=5) return 2;
    const int steps=std::stoi(argv[2]), n=std::stoi(argv[3]), targets=std::stoi(argv[4]);
    std::ifstream file(argv[1],std::ios::binary);
    file.exceptions(std::ios::failbit|std::ios::badbit);
    std::vector<Birth> tags(n), incoming(targets);
    std::vector<double> cosphi(targets), depth(n);
    std::vector<Target> pairs(targets);
    read_array(file,tags);
    read_array(file,cosphi);
    constexpr double M2=8300.*8300.;
    const double factor=6.*std::acos(-1.);
    for(int j=0;j<steps;++j) {
        State h;
        file.read(reinterpret_cast<char*>(&h),sizeof(h));
        read_array(file,incoming);
        const double inv_a=1/h.a;
        for(int k=0;k<targets;++k) {
            const double z=(h.impulse+incoming[k].axial)*inv_a;
            const double q=incoming[k].transverse*inv_a;
            pairs[k]={z,q,std::sqrt(1.+z*z+q*q)};
        }
        for(int i=0;i<n;++i) {
            const double alive_dt = h.dt*0.5 + h.t - tags[i].t;
            if(alive_dt<=0.) continue;
            const double z=(h.impulse+tags[i].axial)*inv_a, q=tags[i].transverse*inv_a;
            const double E=std::sqrt(1.+z*z+q*q);
            double kernel=0.;
            for(int k=0;k<targets;++k) {
                // The negative stream is the momentum mirror of an independent positive-stream draw.
                const auto& p=pairs[k];
                const double product=E*p.E;
                const double dp=std::max(0.,product+z*p.z+q*p.q*cosphi[k]-1.);
                const double s=4.+2.*dp, beta=std::sqrt(dp/(2.+dp));
                const double g=.30*(s+2.), w=beta*g, delta=s-M2;
                kernel+=factor*beta*g*g/((delta*delta+w*w)*product);
            }
            depth[i]+=std::min(h.dt,alive_dt)*h.density*(kernel/targets);
        }
    }
    std::cout<<std::setprecision(17);
    for(double tau:depth) std::cout<<tau<<'\n';
}
