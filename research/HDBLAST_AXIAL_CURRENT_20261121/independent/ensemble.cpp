// Three flavors evolve independently in the same global invariant frame.
// Vector primitives, wall flux stepping, and the verified two-body kernel are supplied by the builder.
using Family = std::array<std::vector<size_t>,2>;
using Membership = std::array<Family,3>;

V rotate_flavor(V v, int index) {
    if(index==1)return {v.z,v.x,v.y};
    if(index==2)return {v.y,v.z,v.x};
    return v;
}
double components(V v){return v.x+v.y+v.z;}
double open_uniform(std::mt19937_64& rng){return (static_cast<double>(rng()>>12)+.5)*0x1p-52;}

struct Tagged { Particle ion; int flavor; };

class GasCollisions {
    std::mt19937_64 rng;
    PoleCollisions pair;
    double density_factor;
    int mode;
public:
    uint64_t events=0, candidates=0, peak_count_range=0;
    std::array<std::array<uint64_t,3>,3> transitions{};
    std::vector<unsigned> hits;
    GasCollisions(size_t count, uint64_t seed,double density,int collision_mode):
        rng(seed),pair(count,0.,seed^0x72617465ULL,false),density_factor(density),mode(collision_mode),hits(count,0) {}
    void advance(std::vector<Tagged>& all,Membership& lists,double duration,double radius,double left,double weight) {
        if(mode==0)return;
        double now=0.;
        while(true){
            std::array<double,3> pairs;
            for(int f=0;f<3;f++)pairs[f]=double(lists[f][0].size())*lists[f][1].size();
            double sum=pairs[0]+pairs[1]+pairs[2];
            if(sum==0.)break;
            double rate=sum*2.*weight*density_factor*1.118569e-6/(radius*radius*radius);
            now-=std::log(open_uniform(rng))/rate;
            if(now>=duration)break;
            const double draw=open_uniform(rng)*sum;
            const int f=draw<pairs[0]?0:(draw<pairs[0]+pairs[1]?1:2);
            auto &plus=lists[f][0], &minus=lists[f][1];
            const size_t ip=static_cast<size_t>(open_uniform(rng)*plus.size());
            const size_t im=static_cast<size_t>(open_uniform(rng)*minus.size());
            const size_t i=plus[ip],j=minus[im];
            if(all[i].ion.birth>left+now || all[j].ion.birth>left+now)continue;
            ++candidates;
            const uint64_t old=pair.event_count();
            pair.encounter(all[i].ion,all[j].ion,radius,i,j);
            if(old==pair.event_count())continue;
            ++events; ++hits[i]; ++hits[j];
            int dest=mode==1?f:static_cast<int>(open_uniform(rng)*3.);
            ++transitions[f][dest];
            if(dest==f)continue;
            // The pair retains the momenta just produced in the global frame.
            // Only its gauge flavor, and hence the field for subsequent evolution, changes.
            plus[ip]=plus.back();plus.pop_back();
            minus[im]=minus.back();minus.pop_back();
            all[i].flavor=all[j].flavor=dest;
            lists[dest][0].push_back(i);lists[dest][1].push_back(j);
            const auto n0=lists[0][0].size(),n1=lists[1][0].size(),n2=lists[2][0].size();
            peak_count_range=std::max<uint64_t>(peak_count_range,std::max({n0,n1,n2})-std::min({n0,n1,n2}));
        }
    }
    void report(const char* name)const {pair.save(name);}
};

int main(int argc,char** argv){
    if(argc!=10)throw std::runtime_error("table tags D density seed mode raw ledger counts");
    auto table=read<14>(argv[1]);auto births=read<9>(argv[2]);
    const double D=std::stod(argv[3]);
    const double w=3.*D/(2.*births.size()); // Effective current weight; total energy uses w/3 per particle.
    std::vector<Tagged> beam;
    for(auto b:births){
        int f=static_cast<int>(b[8]);
        if(f<0 || f>=3)throw std::runtime_error("Invalid gauge index");
        for(int sign:{1,-1})beam.push_back({{sign*rotate_flavor({b[5],b[6],b[7]},f),double(sign),b[0],b[1]},f});
    }
    std::stable_sort(beam.begin(),beam.end(),[](const Tagged& a,const Tagged& b){return a.ion.birth<b.ion.birth;});
    GasCollisions collisions(beam.size(),std::stoull(argv[5]),std::stod(argv[4]),std::stoi(argv[6]));
    const auto init=table.front();
    const double a02=init[1]*init[1];
    std::array<Wave,3> waves; std::array<Budget,3> budgets;
    double initial=0.,injected=0.,expansion=0.,charged_work=0.;
    for(int f=0;f<3;f++){
        waves[f]={a02*rotate_flavor({init[2],init[3],init[4]},f),a02*rotate_flavor({init[5],init[6],init[7]},f)};
        initial+=components(sum_fields(waves[f]))/3.;
    }
    size_t n=0;Membership lists;
    for(size_t j=1;j<table.size();j++){
        const auto left=table[j-1],right=table[j];
        auto mid=[&](int k){return .5*(left[k]+right[k]);};
        const double dt=right[0]-left[0],a=mid(1),de=dt/a;
        V in=(a*a)*V{mid(8),mid(9),mid(10)};
        for(int f=0;f<3;f++)free_wall(waves[f],rotate_flavor(in,f),.5*de,a,budgets[f]);
        const size_t old=n;
        while(n<beam.size() && beam[n].ion.birth<=right[0]){
            const auto& v=beam[n];
            injected+=(w/3.)*energy(v.ion.birth_a,v.ion.u);
            lists[v.flavor][v.ion.sign<0.].push_back(n++);
        }
        for(size_t i=0;i<n;i++){
            auto& v=beam[i].ion;
            double a_start=i<old?left[1]:v.birth_a;
            expansion+=(w/3.)*(a-a_start)*(a+a_start)/(energy(a,v.u)+energy(a_start,v.u));
            v.factor=Q*v.sign*(right[0]-std::max(v.birth,left[0]))/a;
            v.u=rotate(v.u,waves[beam[i].flavor].b,a,.5*v.factor);
            v.g0=energy(a,v.u);
        }
        std::array<V,3> avg={waves[0].e,waves[1].e,waves[2].e};
        for(int it=0;it<100;it++){
            std::array<V,3> current{};
            for(size_t i=0;i<n;i++){
                auto& tag=beam[i];const auto& v=tag.ion;
                const V next=v.u+v.factor*avg[tag.flavor];
                current[tag.flavor]=current[tag.flavor]+(w*v.factor/(v.g0+energy(a,next)))*(v.u+next);
            }
            double residual=0.;
            for(int f=0;f<3;f++){
                const V next=waves[f].e-(.5/K)*current[f];
                residual=std::max(residual,length(next-avg[f])/(1.+length(avg[f])));
                avg[f]=next;
            }
            if(residual<=2.e-14)break;
            if(it==99)throw std::runtime_error("No Maxwell midpoint convergence");
        }
        for(int f=0;f<3;f++)waves[f].e=2.*avg[f]-waves[f].e;
        for(size_t i=0;i<n;i++){
            auto& tag=beam[i];auto& v=tag.ion;const int f=tag.flavor;
            const V next=v.u+v.factor*avg[f];
            charged_work+=(w/3.)*dot((v.factor/(v.g0+energy(a,next)))*(v.u+next),avg[f]);
            v.u=rotate(next,waves[f].b,a,.5*v.factor);
            expansion+=(w/3.)*(right[1]-a)*(right[1]+a)/(energy(a,v.u)+energy(right[1],v.u));
        }
        for(int f=0;f<3;f++)free_wall(waves[f],rotate_flavor(in,f),.5*de,a,budgets[f]);
        collisions.advance(beam,lists,dt,right[1],left[0],w);
    }
    if(n!=beam.size())throw std::runtime_error("Missing births at marker");
    const double a=table.back()[1];double gas=0.,field=0.,flux_in=0.,flux_out=0.;
    std::ofstream raw(argv[7],std::ios::binary),out(argv[8]);
    for(size_t i=0;i<n;i++){
        const auto& tag=beam[i]; const auto& v=tag.ion;
        gas+=(w/3.)*energy(a,v.u);
        std::array<double,6> r={double(tag.flavor),v.sign,v.u.x/a,v.u.y/a,v.u.z/a,double(collisions.hits[i])};
        raw.write(reinterpret_cast<const char*>(r.data()),sizeof(r));
    }
    out<<std::setprecision(17);
    for(int f=0;f<3;f++){
        const auto& wave=waves[f];const auto &budget=budgets[f];
        field+=components(sum_fields(wave))/3.;
        flux_in+=components(budget.incoming)/3.;flux_out+=components(budget.outgoing)/3.;
        const V e=(1./(a*a))*rotate_flavor(wave.e,(3-f)%3);
        out<<e.x<<" "<<e.y<<" "<<e.z<<"\n";
    }
    out<<initial<<" "<<field<<" "<<flux_in<<" "<<flux_out<<" "<<charged_work<<" "<<injected<<" "<<expansion<<" "<<gas<<"\n";
    for(auto& row:collisions.transitions)out<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    out<<collisions.peak_count_range<<"\n";
    collisions.report(argv[9]);
    if(!raw || !out)throw std::runtime_error("Incomplete outputs");
}
