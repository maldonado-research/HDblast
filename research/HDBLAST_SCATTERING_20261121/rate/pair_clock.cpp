// Integrate two independent unscattered daughters on the prescribed homogeneous fields.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Vec { double x,y,z; };
Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vec operator*(double s,Vec a){return {s*a.x,s*a.y,s*a.z};}
double dot(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
Vec cross(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}

template<size_t N> std::vector<std::array<double,N>> read(const char* path){
    std::ifstream f(path,std::ios::binary);
    if (!f) throw std::runtime_error("Missing input file");
    std::vector<std::array<double,N>> rows;
    std::array<double,N> row;
    while(f.read(reinterpret_cast<char*>(row.data()),sizeof(row))) rows.push_back(row);
    if (rows.empty() || !f.eof() || f.gcount()) throw std::runtime_error("Incomplete input table");
    return rows;
}

Vec rotate(Vec p,Vec magnetic,double a,double duration){
    const Vec t=(.5*duration/std::sqrt(a*a+dot(p,p)))*magnetic;
    const Vec v=p+cross(p,t);
    return p+cross(v,(2./(1.+dot(t,t)))*t);
}

Vec step(Vec p,Vec e,Vec b,double a,double duration){
    p=rotate(p,b,a,.5*duration);
    p=p+duration*e;
    return rotate(p,b,a,.5*duration);
}

double pair_kernel(Vec p,Vec q){
    const double energies=std::sqrt((1.+dot(p,p))*(1.+dot(q,q)));
    const double s=std::max(4.,2.*(1.+energies-dot(p,q)));
    const double beta=std::sqrt((s-4.)/s),g=.30*(s+2.);
    const double delta=s-8300.*8300.,width=beta*g;
    return (6.*std::acos(-1.))*beta*g*g/(energies*(delta*delta+width*width));
}

int main(int argc,char** argv){
    if(argc!=3) throw std::runtime_error("Expected field and pair-tag tables");
    const auto field=read<8>(argv[1]), births=read<8>(argv[2]);
    std::vector<Vec> plus,minus;
    for(const auto& r:births){
        plus.push_back({r[1],r[2],r[3]});
        minus.push_back({r[5],r[6],r[7]});
    }
    std::cout<<std::setprecision(17);
    for(size_t j=1;j<field.size();j++){
        const auto& l=field[j-1];
        const auto& r=field[j];
        const double t=.5*(l[0]+r[0]);
        auto mean=[&](int k){return .5*(l[k]+r[k]);};
        const double a=mean(1),a2=a*a, qstep=1.e7/a;
        Vec e=a2*Vec{mean(4),mean(6),mean(2)}, b=a2*Vec{mean(5),mean(7),mean(3)};
        double kernel=0.,maximum=0.;
        size_t live_pairs=0;
        for(size_t i=0;i<births.size();i++){
            const double b1=births[i][0],b2=births[i][4];
            if(b1<t) plus[i]=step(plus[i],e,b,a,qstep*(t-std::max(b1,l[0])));
            if(b2<t) minus[i]=step(minus[i],e,b,a,-qstep*(t-std::max(b2,l[0])));
            if(b1<=t && b2<=t){
                const double k=pair_kernel((1./a)*plus[i],(1./a)*minus[i]);
                kernel+=k;
                maximum=std::max(maximum,k);
                ++live_pairs;
            }
            if(b1<r[0]) plus[i]=step(plus[i],e,b,a,qstep*(r[0]-std::max(b1,t)));
            if(b2<r[0]) minus[i]=step(minus[i],e,b,a,-qstep*(r[0]-std::max(b2,t)));
        }
        std::cout<<kernel/births.size()<<" "<<maximum<<" "<<double(live_pairs)/births.size()<<"\n";
    }
}
