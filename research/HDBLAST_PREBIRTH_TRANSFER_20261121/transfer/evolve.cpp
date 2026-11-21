// A single local Fourier mode of the expanding-wall Maxwell boundary equation.
// Canonical variables p=K a^2 E_cos and b=a^2 B_sin, using wall proper time.
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Row = std::array<double,6>; // t, a, K, i, i_delayed_short, i_delayed_long
struct Mode {
    double k, p=0., b=0., work=0., leaked=0., scalar=0., peak_energy=0.;
    int kind;
    Mode(double wavenumber, int drive): k(wavenumber), kind(drive) {}
    double energy(double K) const { return .25*(p*p/K+K*b*b); }
    void step(double h, double a, double K0, double K1, double incident) {
        double K=.5*(K0+K1);
        // Changing K at fixed canonical fields, then a constant-coefficient implicit midpoint step.
        scalar += .25*(K-K0)*(b*b-p*p/(K*K0));
        double half_angle=.5*h*k/a;
        double source=4.*a*a*incident;
        double pm=(p-half_angle*K*b+.5*h*source)/(1.+h/K+half_angle*half_angle);
        double bm=b+half_angle*pm/K;
        p=2.*pm-p; b=2.*bm-b;
        work+=.5*h*(pm/K)*source;
        leaked+=h*(pm/K)*(pm/K);
        scalar+=.25*(K1-K)*(b*b-p*p/(K1*K));
        peak_energy=std::max(peak_energy,energy(K1));
    }
};

int main(int argc,char** argv) {
    if(argc!=6) throw std::runtime_error("expected trajectory modes birth_time output steps_file");
    std::ifstream input(argv[1],std::ios::binary), specs(argv[2]), targets(argv[5]);
    std::ofstream output(argv[4],std::ios::binary);
    if(!input||!specs||!targets||!output) throw std::runtime_error("file open failed");
    std::vector<Mode> modes;
    double k; int kind; while(specs>>k>>kind) modes.emplace_back(k,kind);
    if(modes.empty()) throw std::runtime_error("no modes");
    double birth=std::stod(argv[3]);
    std::vector<int> keep; int idx; while(targets>>idx) keep.push_back(idx);
    if(keep.empty()||keep.front()<0) throw std::runtime_error("invalid output nodes");
    std::size_t mark=0;
    Row l{},r{};
    if(!input.read(reinterpret_cast<char*>(&l),sizeof(Row))) throw std::runtime_error("empty trajectory");
    auto save=[&](int i) {
        if(mark<keep.size() && keep[mark]==i) {
            for(const auto& m:modes) {
                std::array<double,7> o{m.p/(l[2]*l[1]*l[1]),m.b/(l[1]*l[1]),m.energy(l[2]),
                                       m.work,m.leaked,m.scalar,m.peak_energy};
                output.write(reinterpret_cast<char*>(&o),sizeof(o));
            }
            ++mark;
        }
    };
    int i=0; save(i);
    while(mark<keep.size() && input.read(reinterpret_cast<char*>(&r),sizeof(Row))) {
        double h=r[0]-l[0], a=.5*(r[1]+l[1]);
        if(h<=0. || a<=0. || l[2]<=0. || r[2]<=0.) throw std::runtime_error("nonpositive time/coefficients");
        if(l[0]<birth && birth<r[0]) throw std::runtime_error("birth knot missing from input grid");
        for(auto& mode:modes) {
            bool active=!(mode.kind==1&&l[0]>=birth) && !(mode.kind==2&&r[0]<=birth);
            int col=mode.kind==3?4:(mode.kind==4?5:3);
            mode.step(h,a,l[2],r[2],active*.5*(l[col]+r[col]));
        }
        l=r; save(++i);
    }
    if(!output || mark!=keep.size()) throw std::runtime_error("unwritten output knots");
}
