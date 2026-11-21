#include <iomanip>
#include <iostream>
#include "shell.hpp"
int main(int argc,char** argv) {
    if(argc==2 && std::string(argv[1])=="probe") {
        double t,u;
        while(std::cin>>t>>u) std::cout<<std::setprecision(17)<<shellbath::rate(t,u,24)<<"\n";
        return 0;
    }
    const int scale=argc>1 ? std::stoi(argv[1]) : 1;
    if(scale!=1 && scale!=2) throw std::invalid_argument("Invalid grid resolution");
    // All sampled stream velocities and all positive gas temperatures above 1e-5
    // are contained here. Cooler states use the underlying integral directly.
    const int nt=1+88*scale,ny=1+400*scale;
    const double dt=.125/scale,dy=.05/scale;
    std::cout<<nt<<' '<<ny<<" -5 "<<dt<<' '<<dy<<'\n'<<std::setprecision(13);
    for(int i=0;i<nt;i++) {
        double t=std::pow(10.,-5+i*dt), s=std::sqrt(t/(1.+t));
        for(int j=0;j<ny;j++) {
            double c=shellbath::rate(t,s*std::sinh(j*dy),scale==2?24:12);
            if(!(c>0.)||!std::isfinite(c)) throw std::runtime_error("Invalid drag");
            if(j) std::cout<<' ';
            std::cout<<std::log(c);
        }
        std::cout<<'\n';
    }
}
