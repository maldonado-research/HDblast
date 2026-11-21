// Sum all unordered test-particle pairs. Positive and negative signs alternate in the verified tape.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

template <size_t count> auto rows(const char* path) {
    std::vector<std::array<double,count>> data;
    std::ifstream in(path,std::ios::binary);
    if(!in)throw std::runtime_error("Missing table");
    std::array<double,count> r;
    while(in.read(reinterpret_cast<char*>(r.data()),sizeof(r)))data.push_back(r);
    if(!in.eof()||in.gcount()||data.size()<2)throw std::runtime_error("Bad table");
    return data;
}
struct State{double x,y,z,p,e,boost,invboost;int hot;};
struct Sum{
    double transported=0.,hot=0.,cross_section=0.,fast_boost=0.,hot_fast_boost=0.,small_flux=0.;
    void add(double area,double flux,int top,bool big,bool slow) {
        const double w=area*flux;
        transported+=w;hot+=top*w;cross_section+=area;
        if(big){fast_boost+=w;hot_fast_boost+=top*w;}
        if(slow)small_flux+=w;
    }
};
int main(int argc,char** argv){
    if(argc!=6)throw std::runtime_error("states table bins_per_log cut output");
    auto particles=rows<4>(argv[1]); auto table=rows<7>(argv[2]);
    const double density=std::stod(argv[3]);const double cut=std::stod(argv[4]);
    std::vector<State> v;
    for(auto& r:particles){
        const double p=std::hypot(std::hypot(r[0],r[1]),r[2]),e=std::hypot(1.,p),b=e+p;
        double inv=p==0.?0.:1./p;
        v.push_back({r[0]*inv,r[1]*inv,r[2]*inv,p,e,b,1./b,int(r[3])});
    }
    std::array<std::array<Sum,7>,2> totals{};
    unsigned long long counts[2][2]={};
    for(size_t i=0;i<v.size();i++){
        const auto& a=v[i];
        // One loop per kind to minimize work in the inner pair census.
        for(int like=0;like<2;like++){
            const int first=like?5:0,last=like?7:5;
            for(size_t j=i+1+like;j<v.size();j+=2){
                const auto& b=v[j];
                const double dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;
                const double u=a.boost-b.boost;
                const double rel=.5*u*u*a.invboost*b.invboost+.5*a.p*b.p*(dx*dx+dy*dy+dz*dz);
                const double s=4.+2.*rel;
                const double flux=std::sqrt(rel*(rel+2.))/(a.e*b.e);
                const double position=density*std::log1p(.5*rel);
                const size_t left=size_t(position);const double f=position-left;
                if(left+1>=table.size())throw std::runtime_error("Kernel outside table");
                auto& out=totals[s>cut];++counts[s>cut][like];
                const bool big=(a.e+b.e)*(a.e+b.e)>100.*s,slow=flux<.1;
                for(int col=first;col<last;col++){
                    const double y=table[left][col]+f*(table[left+1][col]-table[left][col]);
                    out[col].add(y,flux,a.hot+b.hot,big,slow);
                }
            }
        }
    }
    std::ofstream out(argv[5]);out<<std::setprecision(17);
    for(int cut=0;cut<2;cut++){
        out<<counts[cut][0]<<" "<<counts[cut][1]<<"\n";
        for(auto& r:totals[cut])out<<r.transported<<" "<<r.hot<<" "<<r.cross_section<<" "<<r.fast_boost<<" "<<r.hot_fast_boost<<" "<<r.small_flux<<"\n";
    }
    if(!out)throw std::runtime_error("Missing output");
}
