#pragma once
#include <algorithm>
#include <cmath>
#include "synge_table.inc"
namespace synge {
struct Thermal { double specific, derivative; };
inline Thermal thermal(double t) {   // random rest-frame energy and its temperature derivative
    if(t<=1e-4) {
        const double b=1.875;
        return {t*(1.5+t*(b-b*t)), 1.5+t*(2*b-3*b*t)};
    }
    if(t>=1e4) return {3*t-1.+.5/t,3.-.5/(t*t)};
    const double inverse_step=50./std::log(10.); // 50 intervals per decade
    const double z=std::log(t)*inverse_step+200.;
    const int i=std::clamp(static_cast<int>(z),0,399);
    const double f=z-i, g=1.-f;
    const double a=synge_nodes[i],b=synge_nodes[i+1];
    const double sa=synge_slopes[i],sb=synge_slopes[i+1];
    const double cv=g*g*(1.+2*f)*a+f*f*(3.-2.*f)*b+f*g*(g*sa-f*sb);
    const double diff=6*f*g*(b-a)+g*(1.-3*f)*sa+f*(3*f-2)*sb;
    return {cv*t,cv+inverse_step*diff};
}
// Solve in tracked excess energy to retain precision when streams are cold and boosted.
inline double invert(double ratio,double extra,double guess) {
    if(extra<=0.) return 0.;
    const double cold=std::hypot(1.,ratio);
    double t=guess;
    for(int j=0;j<15;++j) {
        const auto e=thermal(t);
        const double hm1=e.specific+t, h=1.+hm1, dh=1.+e.derivative;
        const double w=std::hypot(h,ratio), invgamma=h/w;
        const double model=hm1*(2.+hm1)/(w+cold)-t*invgamma;
        const double derivative=((h-t)*dh-h+t*invgamma*invgamma*dh)/w;
        const double delta=(model-extra)/derivative;
        const double next=t-delta;
        if(std::abs(delta)<2e-13*std::max(t,1e-290)) return t;
        t=next>=0. ? next : .5*t;
    }
    return t;
}
}
