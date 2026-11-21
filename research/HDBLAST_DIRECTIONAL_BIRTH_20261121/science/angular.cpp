// Integrate sharp, direction-local pair supply on read-only native time steps.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

constexpr std::array<double,9> edges{-1.,-.5,-.1,-.01,0.,.01,.1,.5,1.};
constexpr double pi=3.14159265358979323846;

// Angles use the positive daughter's direction relative to its electric force.
// h and d are signed radial shell speed and absolute translation speed in units
// of requested normal speed. Pairing the negative daughter at -n gives the
// identical normal speed. Accepted pair and daughter number fractions coincide.
struct Angular {
    double loss=0., pooled_loss=0., crowded=0., p1=0., p2=0., p4=0., signed_p1=0.;
    std::array<double,8> bins{};
    void add(const Angular& x, double w) {
        loss+=w*x.loss; pooled_loss+=w*x.pooled_loss; crowded+=w*x.crowded;
        p1+=w*x.p1; p2+=w*x.p2; p4+=w*x.p4; signed_p1+=w*x.signed_p1;
        for(size_t i=0;i<bins.size();++i) bins[i]+=w*x.bins[i];
    }
};

Angular angular(double h, double signed_d) {
    const double d=std::abs(signed_d),c=std::abs(h);
    Angular out;
    const double mean=d<=c ? c : .5*(d+c*c/d);
    out.pooled_loss=std::max(0.,1.-mean);
    if(c>=d+1.)return out;
    if(d==0.) {
        out.loss=out.pooled_loss; out.crowded=1.;
        for(size_t j=0;j<out.bins.size();++j) out.bins[j]=out.loss*(edges[j+1]-edges[j])*.5;
        return out;
    }
    const double lower=std::clamp((h-1.)/d,-1.,1.),upper=std::clamp((h+1.)/d,-1.,1.);
    const double apex=std::clamp(h/d,lower,upper);
    out.crowded=.5*(upper-lower);
    // Three-point Gauss quadrature exactly integrates a linear tent times P4.
    constexpr double node=.77459666924148337704;
    for(size_t j=0;j<out.bins.size();++j) {
        const double low=std::max(lower,edges[j]),high=std::min(upper,edges[j+1]);
        if(low>=high)continue;
        const double pivot=std::clamp(apex,low,high);
        const std::array<double,3> pts{low,pivot,high};
        for(int i=0;i<2;++i) {
            const double half=.5*(pts[i+1]-pts[i]),mid=.5*(pts[i]+pts[i+1]);
            if(half<=0.)continue;
            for(int g=-1;g<=1;++g) {
                const double u=mid+g*node*half,u2=u*u;
                const double w=half*(g==0?4./9.:5./18.);
                const double lost=w*std::max(0.,1.-std::abs(h-d*u));
                out.bins[j]+=lost; out.loss+=lost; out.p1+=lost*u;
                out.p2+=lost*(3.*u2-1.)/2.; out.p4+=lost*((35.*u2-30.)*u2+3.)/8.;
            }
        }
    }
    out.signed_p1=std::copysign(1.,signed_d)*out.p1;
    return out;
}

void accumulate_segment(Angular& result, double h, double lo, double hi, double weight) {
    constexpr std::array<double,4> node{-.861136311594052575, -.339981043584856265,
                                         .339981043584856265, .861136311594052575};
    constexpr std::array<double,4> w{.347854845137453857, .652145154862546143,
                                     .652145154862546143,.347854845137453857};
    const double mid=.5*(lo+hi),half=.5*(hi-lo);
    for(size_t i=0;i<node.size();++i)result.add(angular(h,mid+half*node[i]),.5*weight*w[i]);
}

Angular step_average(double h, double d0, double d1) {
    if(d0==d1)return angular(h,d0);
    if(d1<d0)std::swap(d0,d1);
    const double range=d1-d0,scale=std::max(1.,std::abs(h));
    Angular avg;
    // Small relative steps already resolve angular features in the integrated
    // loss. Wide steps near a reversal need subdivisions at the physical breaks.
    if(range<.005*std::max({scale,std::abs(d0),std::abs(d1)})) {
        accumulate_segment(avg,h,d0,d1,1.);
        return avg;
    }
    std::vector<double> cuts{d0,d1};
    auto insert=[&](double x){if(x>d0 && x<d1)cuts.push_back(x);};
    insert(0.);
    for(double edge:edges) {
        if(edge==0.)continue;
        for(double shift: {-1.,0.,1.}) {
            const double x=(h-shift)/edge;
            insert(x);insert(-x);
        }
    }
    // The pooled bound changes sign where mean |h - d u| crosses 1.
    if(std::abs(h)<1.) {
        const double root=1.+std::sqrt((1.-h)*(1.+h));
        insert(root);insert(-root);
    }
    // Rational tails decay as powers of 1/|d|. Bound adjacent panel ratios
    // through a large reversal so narrow central losses remain resolved.
    const double maxd=std::max(std::abs(d0),std::abs(d1));
    for(double x=1.;x<maxd;x*=1.5) {insert(x);insert(-x);}
    std::sort(cuts.begin(),cuts.end());
    for(size_t i=1;i<cuts.size();++i) {
        const double width=cuts[i]-cuts[i-1];
        if(width>0.)accumulate_segment(avg,h,cuts[i-1],cuts[i],width/range);
    }
    return avg;
}

struct Window {
    double lower,upper,fraction=0.,peak_loss=0.,peak_tau=0.,peak_fraction=0.,peak_pooled=0.;
    Angular accumulated;
    Window(double low,double high):lower(low),upper(high){}
    void collect(const std::array<double,10>& row,double final_rest,double mass,double action) {
        const auto [t,dt,rest0,rest1,k0,k1,e0,e1,b0,b1]=row;
        const double step_rest=rest1-rest0;
        if(!(step_rest>0 && dt>0))throw std::runtime_error("Invalid step");
        const double u0=std::clamp((lower*final_rest-rest0)/step_rest,0.,1.);
        const double u1=std::clamp((upper*final_rest-rest0)/step_rest,0.,1.);
        if(u0>=u1)return;
        const double k=k0+.5*(u0+u1)*(k1-k0),ct=(k1-k0)/dt;
        const double threshold=2.*pi*pi*action*(step_rest/dt)/(12.*mass*k*k);
        const double h=ct/threshold,dl=(e0+u0*(e1-e0))/threshold,dh=(e0+u1*(e1-e0))/threshold;
        const double population=step_rest*(u1-u0)/final_rest;
        accumulated.add(step_average(h,dl,dh),population);
        fraction+=population;
        const double lo=std::min(dl,dh),hi=std::max(dl,dh);
        const double H=std::abs(h),optimal=H>.5 ? std::sqrt((1.+H)*(1.+H)-2.) : 0.;
        for(double d: {std::clamp(-optimal,lo,hi),std::clamp(optimal,lo,hi)}) {
            // For fixed h the angular integral peaks at |d|=optimal.
            const auto observation=angular(h,d);
            if(observation.loss>peak_loss) {
                peak_loss=observation.loss;
                const double where=dl==dh?0.:std::clamp((d-dl)/(dh-dl),0.,1.);
                const double f=u0+where*(u1-u0);
                peak_tau=t+dt*f; peak_fraction=(rest0+step_rest*f)/final_rest;
            }
        }
        const double near=std::abs(std::clamp(0.,lo,hi));
        const double mean=near<=H?H:.5*(near+H*H/near);
        if(mean>0.)peak_pooled=std::max(peak_pooled,1./mean);
    }
    std::map<std::string,double> values() const {
        std::map<std::string,double> result{
            {"fraction_tested",fraction},{"denied_over_all_requested",accumulated.loss},
            {"denied_over_tested",accumulated.loss/fraction},{"pooled_denied_over_all_requested",accumulated.pooled_loss},
            {"requests_in_crowded_directions_over_all",accumulated.crowded},
            {"peak_denied_over_instantaneous_request",peak_loss},{"peak_denied_at_tau",peak_tau},
            {"peak_denied_at_delivery_fraction",peak_fraction},{"pooled_peak_demand_over_capacity",peak_pooled}};
        const double accepted=fraction-accumulated.loss;
        result["accepted_positive_charge_mean_cosine_along_force"]=-accumulated.p1/accepted;
        result["accepted_positive_charge_mean_signed_cosine"]=-accumulated.signed_p1/accepted;
        result["accepted_partner_mean_cosine_on_same_axis"]=-result["accepted_positive_charge_mean_cosine_along_force"];
        result["accepted_mean_P2"]=-accumulated.p2/accepted;
        result["accepted_mean_P4"]=-accumulated.p4/accepted;
        result["accepted_orthogonal_triad_cubic_moment"]=.4*result["accepted_mean_P4"];
        for(size_t i=0;i<accumulated.bins.size();++i)
            result["bin_"+std::to_string(i)+"_denied_over_all"]=accumulated.bins[i];
        return result;
    }
};

int main(int argc, char** argv) {
    if(argc!=5)throw std::runtime_error("Usage: angular steps.bin final_simulated_rest daughter_mass S0");
    const double rest=std::stod(argv[2]),mass=std::stod(argv[3]),action=std::stod(argv[4]);
    std::ifstream stream(argv[1],std::ios::binary);
    if(!stream)throw std::runtime_error("Unable to read steps");
    Window central(.01,.99),extended(.00001,.99999);
    std::array<double,10> row{};
    while(stream.read(reinterpret_cast<char*>(row.data()),sizeof(row))) {
        extended.collect(row,rest,mass,action);central.collect(row,rest,mass,action);
    }
    if(!stream.eof() || stream.gcount()!=0)throw std::runtime_error("Truncated native readout");
    std::cout<<std::setprecision(17);
    for(auto const& [name,window]:std::map<std::string,Window>{{"central",central},{"extended",extended}}) {
        if(std::abs(window.fraction-(window.upper-window.lower))>1.e-8)
            throw std::runtime_error("Incomplete emission history");
        for(auto const& [key,value]:window.values())std::cout<<name<<" "<<key<<" "<<value<<"\n";
    }
}
