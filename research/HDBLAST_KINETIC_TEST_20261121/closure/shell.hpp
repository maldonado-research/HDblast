#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// Units m=1. An angular average of F sigma, where F is the invariant two-body flux.
namespace shellbath {
constexpr double width = .30, pole2 = 8300.*8300.;
inline double kernel(double q) {  // q = s - 4; entries are nonnegative
    const double beta2=q/(4.+q), g=width*(6.+q), delta=4.+q-pole2;
    return 6.*std::acos(-1.)*std::sqrt(beta2)*g*g/(delta*delta+beta2*g*g);
}
struct Gauss {
    std::vector<double> x,w;
    explicit Gauss(int n): x(n),w(n) {
        for(int i=0;i<(n+1)/2;++i) {
            double z=std::cos(std::acos(-1.)*(i+.75)/(n+.5)), d;
            for(int j=0;j<16;++j) {
                double a=1.,b=z;
                for(int k=2;k<=n;++k){double c=((2*k-1)*z*b-(k-1)*a)/k; a=b;b=c;}
                d=n*(z*b-a)/(z*z-1.);
                const double step=b/d;
                z-=step;
                if(std::abs(step)<2e-16) break;
            }
            x[i]=-z; x[n-1-i]=z;
            w[i]=w[n-1-i]=2./((1.-z*z)*d*d);
        }
    }
};
struct Primitive {
    // log integral of K(q)dq, tabulated using its known derivative for Hermite interpolation.
    static constexpr double low=-70., step=.0125;
    static constexpr int cells=16000;
    std::vector<double> log_i, slope;
    Primitive():log_i(cells+1),slope(cells+1) {
        const Gauss gl(8);
        double q=std::exp(low), area=(2./3.)*q*kernel(q);
        for(int i=0;i<=cells;++i) {
            const double l=low+i*step;
            q=std::exp(l);
            log_i[i]=std::log(area); slope[i]=q*kernel(q)/area;
            if(i==cells) break;
            for(size_t j=0;j<gl.x.size();++j) {
                const double p=std::exp(l+.5*step*(1.+gl.x[j]));
                area+=.5*step*gl.w[j]*p*kernel(p);
            }
        }
    }
    double log_area(double q) const {
        const double l=std::log(q), z=(l-low)/step;
        if(z<0.) return std::log(2./3.)+l+std::log(kernel(q));
        if(z>=cells) throw std::out_of_range("Shell primitive argument outside table");
        const int i=static_cast<int>(z);
        const double f=z-i, g=1.-f;
        return g*g*(1.+2.*f)*log_i[i]+f*f*(3.-2.*f)*log_i[i+1]
               +step*f*g*(g*slope[i]-f*slope[i+1]);
    }
    double mean(double q0,double q1) const {
        if(q1==0.) return 0.;
        if(q1-q0 < 1e-5*q1) return kernel(.5*(q0+q1));
        const double high=log_area(q1);
        const double fraction=q0>0. ? -std::expm1(log_area(q0)-high) : 1.;
        return std::exp(high)*fraction/(q1-q0);
    }
};
inline const Primitive& primitive() {static const Primitive p; return p;}
inline double squared(double v){return v*v;}
inline double rate(double theta,double drift_u,int order=12) {
    const double u=std::abs(drift_u), gamma=std::hypot(1.,u), v=u/gamma, r=std::asinh(u);
    if(theta==0.) return kernel(4.*u*u)/gamma;
    if(theta<0.) throw std::invalid_argument("Negative shell temperature");
    const double eth=1.5*theta+2.25*theta*theta/(1.+std::hypot(1.,1.5*theta));
    const double E=1.+eth, h=E+theta, k=std::sqrt(eth*(2.+eth)), t=std::asinh(k);
    if(theta<1e-10 && v>1e4*k) return kernel(4.*u*u)/gamma;
    const auto& prim=primitive();
    if(v<1e-5*k/E) {
        const double q=4.*k*k;
        return (prim.mean(0.,q)+2.*kernel(q))/(3.*E*h);
    }
    const double hyper=2.*gamma*u, lo=std::abs(t-2*r),hi=t+2*r;
    const double jac=1./(k*hyper);
    const double center=2.*std::acosh(4150.);
    std::vector<double> edges{lo,hi};
    for(double z:{std::abs(t-center),t+center}) {
        for(double d:{-.45,0.,.45})
            if(z+d>lo && z+d<hi) edges.push_back(z+d);
    }
    std::sort(edges.begin(),edges.end());
    static const Gauss standard(12), accurate(24);
    const auto& gl=order==12?standard:accurate;
    double sum=0.;
    for(size_t b=1;b<edges.size();++b) {
        const double span=edges[b]-edges[b-1];
        const int parts=std::max(1,static_cast<int>(std::ceil(span/1.25)));
        const double half=span/(2.*parts);
        for(int i=0;i<parts;++i) {
            const double mid=edges[b-1]+(2*i+1)*half;
            for(size_t j=0;j<gl.x.size();++j) {
                const double a=mid+half*gl.x[j];
                const double z=1.-2.*jac*std::sinh((a+lo)/2)*std::sinh((a-lo)/2);
                const double q0=4.*squared(std::sinh((t-a)/2));
                const double q1=4.*squared(std::sinh((t+a)/2));
                sum+=half*gl.w[j]*(E-k*z/v)*jac*std::sinh(a)*prim.mean(q0,q1);
            }
        }
    }
    return .5*sum/(gamma*h*E*E);
}
}
