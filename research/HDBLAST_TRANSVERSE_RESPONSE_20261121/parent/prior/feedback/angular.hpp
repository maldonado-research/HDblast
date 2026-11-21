// Angular source moments for the charge sign accelerated by +X.
// u is the daughter's physical direction relative to the positive triad axis.
struct BirthMoments {
    double accept=1., first=0.;
};

BirthMoments hard_moments(double c, double d, double target) {
    const double field=std::abs(d);
    if(target<=0. || std::abs(c)>=field+target) return {};
    if(field>=std::abs(c)+target) {
        const double lost=target/(2.*field);
        return {1.-lost, -(c/d)*lost};
    }
    if(d==0.) return {std::min(1.,std::abs(c)/target),0.};
    std::array<double,5> edges{-1.,1.,std::clamp((c-target)/d,-1.,1.),
                          std::clamp(c/d,-1.,1.),std::clamp((c+target)/d,-1.,1.)};
    std::sort(edges.begin(),edges.end());
    BirthMoments sum{0.,0.};
    for(size_t i=1;i<edges.size();++i) {
        const double lo=edges[i-1],hi=edges[i],width=hi-lo;
        // A full saturated piece has value one at both ends, including a
        // rounded root. This keeps its large area independent of root cancellation.
        const bool saturated=std::abs(c-d*(.5*(lo+hi)))>=target;
        const double left=saturated ? 1. : std::min(1.,std::abs(c-d*lo)/target);
        const double right=saturated ? 1. : std::min(1.,std::abs(c-d*hi)/target);
        const double area=width*(left+right)/4.;
        sum.accept+=area;
        sum.first+=.5*(lo+hi)*area+width*width*(right-left)/24.;
    }
    return sum;
}

BirthMoments vacancy_moments(double c,double d,double target) {
    // f=h/(1+h) for equal independent vacancies df/dh=(1-f)^2 on a fresh passage.
    if(target<=0.) return {};
    if(d==0.) return {std::abs(c)/(std::abs(c)+target),0.};
    const double turn=std::clamp(c/d,-1.,1.);
    std::array<double,3> edges{-1.,turn,1.};
    BirthMoments sum{0.,0.};
    for(size_t i=1;i<edges.size();++i) {
        double left=edges[i-1],right=edges[i];
        if(left==right) continue;
        double lo=std::abs(c-d*left),hi=std::abs(c-d*right);
        if(lo>hi) {std::swap(lo,hi); std::swap(left,right);}
        const double width=std::abs(right-left),scale=target+lo;
        const double delta=(hi-lo)/scale, weight=target/scale;
        // Average A and t*A with t increasing linearly from low to high speed.
        const double base=lo/scale;
        double mean, moment;
        if(delta<1e-3) {
            mean=base+weight*delta*(.5+delta*(-1./3.+delta*(.25-delta*.2)));
            moment=.5*base+weight*delta*(1./3.+delta*(-.25+delta*(.2-delta/6.)));
        } else {
            const double log_ratio=std::log1p(delta)/delta;
            mean=1.-weight*log_ratio;
            moment=.5-weight*(1.-log_ratio)/delta;
        }
        sum.accept+=.5*width*mean;
        sum.first+=.5*width*(left*mean+(right-left)*moment);
    }
    return sum;
}

BirthMoments birth_moments(int law, double c, double d, double target) {
    if(law==0 || target<=0.) return {};
    if(law==1) {
        const double abs_c=std::abs(c),abs_d=std::abs(d);
        const double mean=abs_d<=abs_c ? abs_c : .5*(abs_d+abs_c*abs_c/abs_d);
        return {std::min(1.,mean/target),0.};
    }
    return law==2 ? hard_moments(c,d,target) : vacancy_moments(c,d,target);
}
