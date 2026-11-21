// The build inserts the pinned #198 vector types and Lorentz rotation before this file.
#include "collide.hpp"

struct Budget {
    V incoming{},outgoing{},current_work{},source_work{};
};
V squares(V a){return {a.x*a.x,a.y*a.y,a.z*a.z};}
V products(V a,V b){return {a.x*b.x,a.y*b.y,a.z*b.z};}
V sum_fields(Wave w){return (.5*K)*(squares(w.e)+squares(w.b));}

// Advance the wall leakage at fixed radial incident characteristic, analytically.
void impedance(V& electric,V inward,double duration,double radius,Budget& budget) {
    const double lambda=2.*radius/K;
    const double fraction=-std::expm1(-lambda*duration);
    const double fraction2=-std::expm1(-2.*lambda*duration);
    V difference=electric-2.*inward;
    V in=(duration*2.*radius)*squares(inward);
    V out=in+(4.*radius*fraction/lambda)*products(inward,difference)
            +(radius*fraction2/lambda)*squares(difference);
    budget.incoming=budget.incoming+in;
    budget.outgoing=budget.outgoing+out;
    electric=electric-fraction*difference;
}

void free_wall(Wave& wave, V inward,double duration,double radius,Budget& budget) {
    impedance(wave.e,inward,.5*duration,radius,budget);
    const double c=std::cos(kappa*duration),s=std::sin(kappa*duration);
    V old=wave.e;
    wave.e=c*old-s*wave.b;
    wave.b=c*wave.b+s*old;
    impedance(wave.e,inward,.5*duration,radius,budget);
}

V source_current(V predicted,double reference,int mode) {
    if(mode==0)return {0.,0.,reference};  // field replay
    if(mode==1)return {predicted.x,predicted.y,reference};
    return predicted;
}

int main(int argc,char** argv) {
    if(argc!=12)throw std::runtime_error("history births D mode stride density seed counts particles trace ledger");
    const auto table=read<14>(argv[1]);
    const auto births=read<8>(argv[2]);
    const double weight=std::stod(argv[3])/(2*births.size());
    const int mode=std::stoi(argv[4]), stride=std::stoi(argv[5]);
    const double density=std::stod(argv[6]); const uint64_t seed=std::stoull(argv[7]);
    const auto& initial=table.front();
    const double a0=initial[1];
    Wave wave{(a0*a0)*V{initial[2],initial[3],initial[4]},(a0*a0)*V{initial[5],initial[6],initial[7]}};
    V first_energy=sum_fields(wave); Budget budget;
    double injection=0.,expansion=0.,residual=0.;
    int largest=0;
    std::vector<Particle> particles;
    for(const auto& birth:births)
        for(int sign:{1,-1}) particles.push_back({sign*V{birth[5],birth[6],birth[7]},double(sign),birth[0],birth[1]});
    std::stable_sort(particles.begin(),particles.end(),[](const Particle& a,const Particle& b){return a.birth<b.birth;});
    PoleCollisions collisions(particles.size(),density,seed,false);
    size_t n=0;
    std::ofstream trace(argv[10],std::ios::binary);
    auto record=[&](double t,double a,V current,V used) {
        double a2=a*a;
        // The field is conformal on integration, physical in the output. Current is conformal.
        std::array<double,13> row={t,wave.e.x/a2,wave.e.y/a2,wave.e.z/a2,wave.b.x/a2,wave.b.y/a2,wave.b.z/a2,
                                   current.x,current.y,current.z,used.x,used.y,used.z};
        trace.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    };
    record(initial[0],a0,{},{});
    for(size_t j=1;j<table.size();++j){
        const auto& l=table[j-1];const auto& r=table[j];
        const double dt=r[0]-l[0],a=.5*(r[1]+l[1]),deta=dt/a;
        if(!(dt>0.))throw std::runtime_error("time grid");
        auto mid=[&](int k){return .5*(l[k]+r[k]);};
        const V inward=(a*a)*V{mid(8),mid(9),mid(10)};
        free_wall(wave,inward,.5*deta,a,budget);
        const size_t old=n;
        while(n<particles.size() && particles[n].birth<=r[0]){
            injection+=weight*energy(particles[n].birth_a,particles[n].u);++n;
        }
        for(size_t i=0;i<n;++i){
            auto& v=particles[i];
            double start=i<old?l[1]:v.birth_a;
            v.factor=Q*v.sign*(r[0]-std::max(v.birth,l[0]))/a;
            expansion+=weight*(a-start)*(a+start)/(energy(a,v.u)+energy(start,v.u));
            v.u=rotate(v.u,wave.b,a,.5*v.factor);
            v.g0=energy(a,v.u);
        }
        const double ref=deta*mid(11);
        V avg=wave.e,current{};
        for(int it=1;it<=100;++it){
            current={};
            for(size_t i=0;i<n;++i){
                const auto& v=particles[i];
                V after=v.u+v.factor*avg;
                current=current+(weight*v.factor/(v.g0+energy(a,after)))*(v.u+after);
            }
            const V next=wave.e-(.5/K)*source_current(current,ref,mode);
            const double mismatch=length(next-avg);
            largest=std::max(largest,it);
            avg=next;
            if(mismatch<=2.e-14*(1.+length(avg))){residual=std::max(residual,mismatch);break;}
            if(it==100)throw std::runtime_error("implicit Maxwell solve convergence");
        }
        wave.e=2.*avg-wave.e;
        current={};
        for(size_t i=0;i<n;++i){
            auto& v=particles[i];
            const V next=v.u+v.factor*avg;
            const V flux=(weight*v.factor/(v.g0+energy(a,next)))*(v.u+next);
            current=current+flux;
            v.u=rotate(next,wave.b,a,.5*v.factor);
            expansion+=weight*(r[1]-a)*(r[1]+a)/(energy(a,v.u)+energy(r[1],v.u));
        }
        const V used=source_current(current,ref,mode);
        budget.current_work=budget.current_work+products(current,avg);
        budget.source_work=budget.source_work+products(used,avg);
        free_wall(wave,inward,.5*deta,a,budget);
        collisions.step(particles,n,weight,l[0],r[0],r[1]);
        if(j%stride==0 || j+1==table.size())record(r[0],r[1],(1./deta)*current,(1./deta)*used);
    }
    if(n!=particles.size())throw std::runtime_error("Incomplete births");
    const double final_a=table.back()[1];
    double gas=0.;
    std::ofstream output(argv[9],std::ios::binary);
    for(const auto& v:particles){
        gas+=weight*energy(final_a,v.u);
        std::array<double,4> row={v.sign,v.u.x/final_a,v.u.y/final_a,v.u.z/final_a};
        output.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    }
    collisions.save(argv[8]);
    std::ofstream ledger(argv[11]); ledger<<std::setprecision(17);
    for(V v:{first_energy,sum_fields(wave),budget.incoming,budget.outgoing,budget.source_work,budget.current_work})
        ledger<<v.x<<" "<<v.y<<" "<<v.z<<"\n";
    ledger<<gas<<" "<<injection<<" "<<expansion<<" "<<largest<<" "<<residual<<"\n";
    if(!trace||!output||!ledger)throw std::runtime_error("Incomplete output");
}
