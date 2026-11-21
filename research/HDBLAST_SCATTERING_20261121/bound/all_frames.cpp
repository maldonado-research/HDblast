// Enclose the energy-threshold discrepancy for any constant three-velocity.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using Four = std::array<double,4>;
using Three = std::array<double,3>;
constexpr int cuts = 96;
constexpr double spacing = 0.125;
constexpr double rounding = 1.e-7;

struct Cover {
    std::vector<Four> sample;
    double theta;
    double beta;
    double requested;
    std::vector<double> sizes;
    long long visited=0, accepted=0, exterior=0;
    int max_depth=0;
    double certified=1.;
    Three weakest{};
    double weakest_halfwidth=0.;
    int weakest_cut=0;
    int weakest_direction=0;

    // Primitive of exp(-x/k)*(1+x/k+(x/k)^2/2)*(B - h*offset + h*k).
    static double primitive(double x, double k, double speed, double slope, double offset) {
        if (k<=0.) return 0.;
        return std::exp(-x/k)*((speed-slope*offset)*(k+.5*x)+.5*slope*k*(k+x));
    }

    static double massless_cdf(double x, double speed, double slope, double offset) {
        if (speed==0.) {
            if (offset<=0.) return 1.;
            double z=x/offset;
            return 1.-std::exp(-z)*(1.+z+.5*z*z);
        }
        double upper=primitive(x,offset+speed,speed,slope,offset);
        double lower=primitive(x,offset-speed,speed,slope,offset);
        return 1.-(upper-lower)/(2.*speed*speed);
    }

    static int event_bin(double value) {
        return std::clamp(static_cast<int>(std::ceil(value/spacing)-1.),0,cuts);
    }

    void visit(Three center, double radius, int depth) {
        ++visited;
        double smallest_squared=0.;
        for (double coordinate:center) {
            double gap=std::max(0.,std::abs(coordinate)-radius);
            smallest_squared+=gap*gap;
        }
        if (smallest_squared>1.+rounding) { ++exterior; return; }
        max_depth=std::max(max_depth,depth);
        std::array<int,cuts+1> empirical_upper{},empirical_lower{};
        for (std::size_t i=0;i<sample.size();++i) {
            const auto &p=sample[i];
            double mid=p[0]+center[0]*p[1]+center[1]*p[2]+center[2]*p[3];
            double margin=radius*sizes[i]+rounding;
            ++empirical_lower[event_bin(mid+margin)];
            ++empirical_upper[event_bin(mid-margin)];
        }
        double speed=std::hypot(center[0],center[1],center[2]);
        double delta=std::sqrt(3.)*radius;
        double slope=speed>0.?beta*center[2]/speed:0.;
        double epsilon=1./theta;
        double density_margin=.5*std::expm1(epsilon)+epsilon/8.+rounding;
        int count_upper=0,count_lower=0,best_cut=0,best_dir=0;
        double bound=0.;
        for (int k=0;k<cuts;++k) {
            count_upper+=empirical_upper[k];
            count_lower+=empirical_lower[k];
            double threshold=spacing*(k+1);
            double high=massless_cdf(threshold,speed,slope,1.-delta)+density_margin;
            double low=massless_cdf(threshold-epsilon,speed,slope,1.+delta)-density_margin;
            double obs_excess=double(count_lower)/sample.size()-high;
            double target_excess=low-double(count_upper)/sample.size();
            if (obs_excess>bound) { bound=obs_excess; best_cut=k+1; best_dir=1; }
            if (target_excess>bound) { bound=target_excess; best_cut=k+1; best_dir=-1; }
        }
        if (bound>=requested) {
            ++accepted;
            if (bound<certified) {
                certified=bound; weakest=center; weakest_halfwidth=radius;
                weakest_cut=best_cut; weakest_direction=best_dir;
            }
            return;
        }
        if (depth>=10) throw std::runtime_error("Unresolved boost region");
        double next=radius*.5;
        for (int child=0;child<8;++child) {
            Three midpoint=center;
            for (int j=0;j<3;++j) midpoint[j]+=(child&(1<<j))?next:-next;
            visit(midpoint,next,depth+1);
        }
    }
};

int main(int argc, char **argv) {
    if (argc==2 && std::string(argv[1])=="--check-cdf") {
        double x,speed,slope,offset;
        std::cout<<std::setprecision(16);
        while (std::cin>>x>>speed>>slope>>offset)
            std::cout<<Cover::massless_cdf(x,speed,slope,offset)<<"\n";
        return 0;
    }
    if (argc!=5) throw std::runtime_error("usage: all_frames binary_four_rows theta beta required_gap");
    std::ifstream in(argv[1],std::ios::binary);
    Cover c;
    c.theta=std::stod(argv[2]); c.beta=std::stod(argv[3]); c.requested=std::stod(argv[4]);
    Four row;
    while(in.read(reinterpret_cast<char*>(row.data()),sizeof(row))) {
        c.sample.push_back(row);
        c.sizes.push_back(std::abs(row[1])+std::abs(row[2])+std::abs(row[3]));
    }
    if (c.sample.empty() || !in.eof() || in.gcount()!=0) throw std::runtime_error("Bad snapshot bytes");
    c.visit({},1.,0);
    if (8*(c.visited-c.accepted-c.exterior)+1!=c.visited) throw std::runtime_error("Octree count");
    std::cout<<std::setprecision(14)<<c.certified<<" "<<c.visited<<" "<<c.accepted<<" "<<c.exterior<<" "
        <<c.max_depth<<" "<<c.weakest[0]<<" "<<c.weakest[1]<<" "<<c.weakest[2]<<" "
        <<c.weakest_halfwidth<<" "<<c.weakest_cut<<" "<<c.weakest_direction<<"\n";
}
