// A conservative wall-energy variable shares the kinetic particle current and pressure
// with the field and gravity. The directly sampled lab energy is retained independently.
struct ParticleStress {
    double current=0., pressure=0., energy=0., intrinsic=0., eckart=0., average_v=0.;
    double longitudinal_share=0., landau_velocity=0., proper_particle_energy=1.;
    double lab_dispersion=0., susceptibility=0.;
};
kin::Wall* cloud=nullptr;
int kinetic_particles=0;
unsigned generation=0;
double last_a=-1., last_I=0.;
unsigned last_generation=0;
bool last_detailed=false;
ParticleStress cached;
ParticleStress stress(double a,double I,bool detailed=false) {
    if (a==last_a && I==last_I && generation==last_generation && (!detailed || last_detailed)) return cached;
    last_a=a; last_I=I; last_generation=generation; last_detailed=detailed;
    cached=ParticleStress{};
    if(cloud->ps.empty()) return cached;
    double e=0.,trace=0.,v=0.,ep=0.,em=0.,np=0.,nm=0.,long_difference=0.,transverse_trace=0.,susceptibility=0.,longitudinal_momentum=0.;
    const double scale=1./a;
    if(!detailed) {
        for(const auto& p:cloud->ps) {
            const double z=(p.qz+I)*scale,b=p.qperp*scale,q2=z*z+b*b;
            const double inv=1./std::sqrt(1.+q2);
            trace+=q2*inv;v+=z*inv;
        }
        const double den=cloud->w/(a*a*a);
        cached.current=den*charge_mass*v; cached.pressure=den*trace/3.;
        cached.average_v=v/cloud->ps.size();
        return cached;
    }
    for(const auto& p:cloud->ps) {
        const double z=(p.qz+I)*scale, b=p.qperp*scale, transverse=1.+b*b;
        const double energy=std::sqrt(transverse+z*z), inv=1./energy;
        const double big=energy+std::abs(z);
        const double small=transverse/big;
        const double plus=z>=0 ? big : small, minus=z>=0 ? small : big;
        ep+=plus*plus*inv; em+=minus*minus*inv;
        np+=plus*inv; nm+=minus*inv;
        long_difference+=transverse*inv;transverse_trace+=b*b*inv;
        susceptibility+=transverse*inv*inv*inv;longitudinal_momentum+=z;
        e+=energy; trace+=(b*b+z*z)*inv;v+=z*inv;
    }
    const double den=cloud->w/(a*a*a), n=cloud->ps.size();
    const double ex=std::sqrt(std::sqrt(ep/em));
    const double comoving_count=.5*(ex*nm+np/ex);
    const double proper_energy=.5*(long_difference+std::sqrt(ep*em));
    const double boost_n=std::sqrt(np/nm);
    const double e_eckart=.25*(em*boost_n*boost_n+ep/(boost_n*boost_n)+2*long_difference);
    const double n_eckart=std::sqrt(np*nm);
    cached.current=den*charge_mass*v; cached.pressure=den*trace/3.;
    cached.energy=den*e; cached.average_v=v/n;
    const double pzz=std::max(0.,.5*(std::sqrt(ep*em)-long_difference));
    cached.longitudinal_share=(pzz+transverse_trace)>0. ? pzz/(pzz+transverse_trace) : 0.;
    cached.landau_velocity=(ex*ex-1.)/(ex*ex+1.);
    cached.lab_dispersion=den*std::max(0.,e-std::hypot(n,longitudinal_momentum));
    cached.susceptibility=den*susceptibility;
    cached.proper_particle_energy=proper_energy/comoving_count;
    cached.intrinsic=den*n*std::max(0.,cached.proper_particle_energy-1.);
    cached.eckart=den*n*std::max(0.,e_eckart/n_eckart-1.);
    return cached;
}

struct Interval {
    double start=-1.,growth=1.,start_radius=0.,end_radius=0.,qualified_log=0.,qualified_time=0.;
    int intervals=0;
    void update(double dt,double a0,double a1,double h0,double h1,double phi0,double phi1,
                double hot0,double hot1,double q0=1./3.,double q1=1./3.) {
        double lo=0.,hi=1.;
        auto intersect=[&](double x,double y) {
            if(x<0. && y<0.) {lo=1.;hi=0.;}
            else if(x<0.) lo=std::max(lo,-x/(y-x));
            else if(y<0.) hi=std::min(hi,x/(x-y));
        };
        intersect(h0,h1);intersect(phi0,phi1);intersect(hot0-.9,hot1-.9);
        intersect(q0-(1./3.-.05),q1-(1./3.-.05));intersect(1./3.+.05-q0,1./3.+.05-q1);
        if(hi<=lo) {start=-1.; return;}
        if(kinetic_particles>=2048) {
            qualified_log+=std::log1p((hi-lo)*(a1-a0)/(a0+lo*(a1-a0)));
            qualified_time+=(hi-lo)*dt;
        }
        if(lo>0. || start<0.) {start=a0+lo*(a1-a0);++intervals;}
        const double end=a0+hi*(a1-a0);
        if(end/start>growth) {growth=end/start;start_radius=start;end_radius=end;}
        if(hi<1.) start=-1.;
    }
    void write(std::map<std::string,double>& out,const std::string& prefix) const {
        out[prefix+"_growth"]=growth;out[prefix+"_intervals"]=intervals;
        if(kinetic_particles>=2048) {
            out[prefix+"_qualified_log_growth"]=qualified_log;
            out[prefix+"_qualified_time"]=qualified_time;
        }
        if(intervals>0) {out[prefix+"_start_radius"]=start_radius;out[prefix+"_end_radius"]=end_radius;}
    }
};
Interval landau_angular,lab_hot,lab_angular,qualifying_landau,post_release;
void snapshot(std::map<std::string,double>& out,const std::string& name,double a,double impulse,double wall_energy) {
    const auto s=stress(a,impulse,true);
    if(cloud->ps.empty() || s.proper_particle_energy<=1.) return;
    std::vector<distribution::Observation> obs;
    obs.reserve(cloud->ps.size());
    for(const auto& p:cloud->ps) obs.push_back({p.qperp/a,(p.qz+impulse)/a});
    const auto shape=distribution::evaluate(obs,s.landau_velocity,s.proper_particle_energy);
    out[name+"_lab_dispersion_fraction"]=s.lab_dispersion/wall_energy;
    out[name+"_shape_ks"]=shape.ks;
    out[name+"_angular_trace_share"]=shape.longitudinal_share;
    out[name+"_fitted_temperature_over_m"]=shape.temperature;
    out[name+"_empirical_tube_share"]=shape.tube;
    out[name+"_shape_pressure_check"]=std::abs(shape.longitudinal_share-s.longitudinal_share);
}
