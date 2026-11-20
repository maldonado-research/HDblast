// Read-only diagnostics of local Weyl phase-space flux. All calculations use the
// proper wall clock. The signed electric acceleration is linear on each accepted step.
struct SweepIntegral {
    double mean=0., unmet=0., deficient_time=0.;
};

double directional_mean(double d, double c) {
    const double mag=std::abs(d);
    return mag<=c ? c : .5*(mag+c*c/mag);
}

SweepIntegral sweep_integral(double d0,double d1,double c,double target) {
    if(d0==d1) {
        const double mean=directional_mean(d0,c);
        return {mean,std::max(0.,target-mean),target>mean ? 1. : 0.};
    }
    if(d0>d1) std::swap(d0,d1);
    const double edge=target>c ? target+std::sqrt((target-c)*(target+c)) : 0.;
    double cuts[7]={d0,d1,0.,-c,c,-edge,edge};
    std::sort(cuts,cuts+7);
    SweepIntegral result;
    for(int j=0;j<6;++j) {
        const double l=std::max(d0,cuts[j]), h=std::min(d1,cuts[j+1]);
        if(h<=l) continue;
        const double range=(h-l)/(d1-d0);
        const double vmin=std::min(std::abs(l),std::abs(h)),vmax=std::max(std::abs(l),std::abs(h));
        double mean=c;
        if(.5*(vmin+vmax)>c) {
            mean=.25*(vmin+vmax);
            if(c>0.) mean+=.5*c*c*std::log1p((vmax-vmin)/vmin)/(vmax-vmin);
        }
        result.mean+=range*mean;
        if(target>mean) {
            result.unmet+=range*(target-mean);result.deficient_time+=range;
        }
    }
    return result;
}

void check_sweep_integral() {
    const auto k=sweep_integral(-3.,3.,1.,1.25);
    const auto z=sweep_integral(-2.,2.,0.,.5);
    if(std::abs(k.mean-(3.+.5*std::log(3.))/3.)>1.e-12 ||
       std::abs(k.unmet-(.75-.5*std::log(2.))/3.)>1.e-12 ||
       std::abs(k.deficient_time-2./3.)>1.e-12 ||
       std::abs(z.mean-.5)>1.e-12 || std::abs(z.unmet-.125)>1.e-12 ||
       std::abs(z.deficient_time-.5)>1.e-12)
        throw std::runtime_error("Electric-crossing integral disagrees with analytic controls");
}

struct SweepBin {
    explicit SweepBin(double w):width(w){}
    double width,request=0.,capacity=0.,coverage=0.,highest=0.,worst_start=0.;
    int index=-1,complete=0,violated=0;
    void finish() {
        if(coverage<width*(1.-1e-7)) return;
        ++complete;
        const double ratio=request/capacity;
        if(ratio>1.) ++violated;
        if(ratio>highest) {highest=ratio;worst_start=index*width;}
    }
    void add(double t,double dt,double d0,double d1,double c,double k,double request_rate) {
        double start=t;
        while(start<t+dt) {
            const int next=static_cast<int>(std::floor(start/width+1e-9));
            if(next!=index) {finish();index=next;request=capacity=coverage=0.;}
            const double end=std::min(t+dt,(index+1)*width);
            const double part=end-start;
            if(part<=0.) break;
            const double lower=d0+(d1-d0)*(start-t)/dt, upper=d0+(d1-d0)*(end-t)/dt;
            const auto integral=sweep_integral(lower,upper,c,0.);
            request+=part*request_rate;capacity+=part*12.*k*k*integral.mean;coverage+=part;
            start=end;
        }
    }
};

struct SweepCollector {
    double peak=0.,peak_time=0.,peak_fraction=0.,native_peak=0.,first_time=-1.,last_time=0.;
    double absorbed=0.,during=0.,delivered=0.,coarse_absorbed=0.,total_crossings=0.,hot_zeroes=0.;
    std::array<SweepBin,3> bins{SweepBin(.0001),SweepBin(.001),SweepBin(.01)};
    void collect(double t,double dt,double d0,double d1,double c,double k,double source,double drest,double f0,double f1) {
        if(first_time<0.) first_time=t;
        last_time=t+dt;
        const double target=source/(12.*k*k);
        const double u=d0==d1?0.:std::clamp(-d0/(d1-d0),0.,1.);
        const double min_d=(d0*d1<=0.)?0.:std::min(std::abs(d0),std::abs(d1));
        const double ratio=target/directional_mean(min_d,c);
        if(ratio>peak) {peak=ratio;peak_time=t+dt*u;peak_fraction=f0+(f1-f0)*u;}
        if(d0*d1<0.) {++total_crossings;if(target>c) ++hot_zeroes;}
        const double native=target/directional_mean(d1,c);
        native_peak=std::max(native_peak,native);
        coarse_absorbed+=drest*std::max(0.,1.-1./native);
        const auto integral=sweep_integral(d0,d1,c,target);
        absorbed+=drest*integral.unmet/target;
        during+=drest*integral.deficient_time;
        delivered+=drest;
        for(auto &bin:bins) bin.add(t,dt,d0,d1,c,k,source);
    }
    void save(std::map<std::string,double>& out,const std::string& prefix,double final_rest) {
        out[prefix+"start_tau"]=first_time;out[prefix+"end_tau"]=last_time;
        out[prefix+"peak"]=peak;out[prefix+"peak_tau"]=peak_time;out[prefix+"peak_fraction"]=peak_fraction;
        out[prefix+"native_peak"]=native_peak;out[prefix+"deficit_fraction"]=absorbed/final_rest;
        out[prefix+"emitted_fraction_during_deficit"]=during/final_rest;out[prefix+"fraction_tested"]=delivered/final_rest;
        out[prefix+"native_deficit_fraction"]=coarse_absorbed/final_rest;
        out[prefix+"zero_crossings"]=total_crossings;out[prefix+"zero_crossings_with_deficit"]=hot_zeroes;
        for(auto &b:bins) {
            b.finish();
            const std::string p=prefix+"width_"+std::to_string(b.width)+"_";
            out[p+"largest_ratio"]=b.highest;out[p+"start"]=b.worst_start;
            out[p+"complete_bins"]=b.complete;out[p+"over_capacity_bins"]=b.violated;
        }
    }
};
