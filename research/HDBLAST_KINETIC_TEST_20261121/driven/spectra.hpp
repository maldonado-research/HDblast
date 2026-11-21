#pragma once
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace distribution {
struct Observation { double transverse, longitudinal; };
struct Shape {double temperature, ks, longitudinal_share, tube, thermal_tube;};
struct Reference {
    static constexpr int bins=1536;
    std::vector<double> cdf;
    double mean;
    explicit Reference(double t):cdf(bins+1) {
        double total=0.,energy=0.;
        // x=(E-1)/Theta, sample at midpoints in sqrt(x) to regularize the endpoint.
        for(int i=0;i<bins;++i) {
            const double s=(i+.5)*8./bins,x=s*s;
            const double weight=x*(1.+t*x)*std::sqrt(2.+t*x)*std::exp(-x);
            total+=weight; energy+=t*x*weight;
            cdf[i+1]=total;
        }
        mean=energy/total;
        for(auto& f:cdf) f/=total;
    }
    double probability(double x) const {
        const double b=bins/8.*std::sqrt(std::max(0.,x));
        if(b>=bins)return 1.;
        const int i=int(b);
        return cdf[i]+(b-i)*(cdf[i+1]-cdf[i]);
    }
};
inline Shape evaluate(const std::vector<Observation>& particles,double v,double proper_energy) {
    const double thermal=proper_energy-1.,gamma=1./std::sqrt((1.-v)*(1.+v));
    double lo=thermal/3.,hi=2.*thermal/3.;
    for(int j=0;j<30;++j) {
        const double t=.5*(lo+hi);
        if(Reference(t).mean>thermal) hi=t; else lo=t;
    }
    const double t=.5*(lo+hi);
    Reference gas(t);
    double sum=0.,trace=0.,longitudinal=0.,in_tube=0.;
    std::vector<std::pair<double,double>> sample;
    for(const auto& p:particles) {
        const double transverse=p.transverse, z=p.longitudinal;
        const double e=std::hypot(std::hypot(1.,transverse),z);
        const double ez=gamma*(z-v*e), eprime=std::hypot(std::hypot(1.,transverse),ez);
        const double w=eprime/e;
        sum+=w;trace+=w*(ez*ez+transverse*transverse)/eprime;
        longitudinal+=w*ez*ez/eprime;
        in_tube+=transverse*transverse<=.1025?1.:0.;
        sample.push_back({(eprime-1.)/t,w});
    }
    std::sort(sample.begin(),sample.end());
    double cum=0.,ks=0.;
    for(const auto& entry:sample) {
        const double F=gas.probability(entry.first);
        ks=std::max(ks,std::abs(cum/sum-F));
        cum+=entry.second;
        ks=std::max(ks,std::abs(cum/sum-F));
    }
    // MJ transverse marginal is boost invariant. Integrate the probability of being
    // in an axial cylinder at fixed |p|, using 1-sqrt(1-y)=y/(1+sqrt(1-y)).
    double thermal_tube=0.;
    for(int i=0;i<Reference::bins;++i) {
        const double s=(i+.5)*8./Reference::bins,ex=t*s*s;
        const double y=std::min(1.,.1025/(ex*(2.+ex)));
        thermal_tube+=(gas.cdf[i+1]-gas.cdf[i])*y/(1.+std::sqrt(1.-y));
    }
    return {t,ks,longitudinal/trace,in_tube/particles.size(),thermal_tube};
}
}
