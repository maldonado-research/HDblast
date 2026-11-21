// Collisionless history of paired canonical labels in the sharp, isotropic Weyl approximation.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

struct Point { double t,k,b,D; };
Point mix(Point a, Point b, double s) {
    return {a.t+s*(b.t-a.t),a.k+s*(b.k-a.k),a.b+s*(b.b-a.b),a.D+s*(b.D-a.D)};
}
struct HitSum {
    double H=0., capped=0., capped2=0.;
    unsigned n=0;
    std::vector<std::pair<double,double>> visits;
    void add(double h,double time) {visits.emplace_back(time,h); double q=std::min(1.,h); H+=h;capped+=q;capped2+=q*q;++n; }
};
double extra_after_clearance(HitSum &sum,double rate) {
    std::sort(sum.visits.begin(),sum.visits.end());
    double occupied=0.,last=sum.visits.front().first,extra=0.;
    for(const auto &visit:sum.visits) {
        occupied*=std::exp(-rate*(visit.first-last));
        const double clipped=std::min(1.,visit.second);
        extra+=std::max(0.,occupied+clipped-1.);
        occupied=std::min(1.,occupied+clipped);last=visit.first;
    }
    return extra;
}
struct Node { int start,end,left=-1,right=-1; double klo,khi,blo,bhi; };
struct Event { int i; double s,k,u,axial,perp2,h; };

class History {
    std::vector<Point> x;
    std::vector<Node> tree;
    double load,local_average;
    std::vector<double> time_cdf;
    int make_node(int first, int end) {
        int id=tree.size();tree.push_back({first,end});
        Node n{first,end};
        if(end-first<=6) {
            n.klo=n.khi=x[first].k;n.blo=n.bhi=x[first].b;
            for(int i=first+1;i<=end;++i) {
                n.klo=std::min(n.klo,x[i].k);n.khi=std::max(n.khi,x[i].k);
                n.blo=std::min(n.blo,x[i].b);n.bhi=std::max(n.bhi,x[i].b);
            }
        } else {
            int mid=(first+end)/2;
            n.left=make_node(first,mid);n.right=make_node(mid,end);
            const Node &l=tree[n.left], &r=tree[n.right];
            n.klo=std::min(l.klo,r.klo);n.khi=std::max(l.khi,r.khi);
            n.blo=std::min(l.blo,r.blo);n.bhi=std::max(l.bhi,r.bhi);
        }
        tree[id]=n;return id;
    }
    void hits_in_segment(int i,const Event &e,HitSum &sum) const {
        const Point &p=x[i], &q=x[i+1];
        const long double db=q.b-p.b,dk=q.k-p.k,amount=load*(q.D-p.D);
        if(amount<=0.) return;
        const Point &eb=x[e.i], &ec=x[e.i+1];
        const long double center_delta=(p.b-eb.b)-static_cast<long double>(e.s)*(ec.b-eb.b);
        const long double radius_delta=(eb.k-p.k)+static_cast<long double>(e.s)*(ec.k-eb.k);
        const long double birth_k=p.k+radius_delta, z=center_delta+birth_k*e.u;
        const long double A=(db-dk)*(db+dk);
        const long double B=2*(z*db-p.k*dk);
        auto add=[&](long double s,long double deriv) {
            if(s>=0 && s<1.) {
                long double k=p.k+s*dk;
                sum.add(static_cast<double>(2*amount/(k*std::abs(deriv))),p.t+static_cast<double>(s)*(q.t-p.t));
            }
        };
        if(i==e.i) {
            // The sampled visit is known exactly. Factor it out before finding the other root.
            sum.add(e.h,p.t+e.s*(q.t-p.t));
            const long double ds=2*e.k*(db*e.u-dk);
            if(A!=0.) add(e.s-ds/A,-ds);
            return;
        }
        // Expand relative to the sampled birth radius to retain grazing angles far below one.
        const long double C=center_delta*(center_delta+2*birth_k*e.u)+radius_delta*(birth_k+p.k);
        if(A==0.) { if(B!=0.)add(-C/B,B);return; }
        const long double disc=B*B-4*A*C;
        if(disc<=0.) return; // A double root has zero measure in label space.
        const long double root=std::sqrt(disc), w=-.5L*(B+std::copysign(root,B));
        add(w/A,root);if(w!=0.)add(C/w,root);
    }
    void visit(int id,const Event &e,HitSum &sum) const {
        const auto &n=tree[id];
        const double lo=e.axial+n.blo,hi=e.axial+n.bhi;
        const double far=std::max(std::abs(lo),std::abs(hi));
        const double near=std::max({0.,lo,-hi});
        // Guard root membership at roundoff scale after absolute canonical center arithmetic.
        const double pad=1.e-8*(1.+n.khi);
        if(near*near+e.perp2>(n.khi+pad)*(n.khi+pad) ||
           far*far+e.perp2<std::max(0.,n.klo-pad)*std::max(0.,n.klo-pad))return;
        if(n.left>=0) {visit(n.left,e,sum);visit(n.right,e,sum);return;}
        for(int i=n.start;i<n.end;++i) hits_in_segment(i,e,sum);
    }
    static double local_loss(double T,double c,double d) {
        if(d==0.)return std::max(0.,1-c/T);
        double area=0.;
        for(double sign: {1.,-1.}) {
            const double l=std::max(0.,sign*c-d),r=std::min(T,sign*c+d);
            if(l<r)area+=(r-l)*(1.-(r+l)/(2*T));
        }
        return area/(2*d);
    }
public:
    History(std::vector<Point> p,double weight):x(std::move(p)),load(weight) {
        if(x.size()<2)throw std::runtime_error("Empty history");
        make_node(0,x.size()-1);
        local_average=0.;
        time_cdf.push_back(0.);
        for(size_t i=0;i+1<x.size();++i){
            const Point &p=x[i], &q=x[i+1];double dd=q.D-p.D;
            if(dd<=0)throw std::runtime_error("Birth history must be increasing");
            // Two-point Gauss rule in time within this segment.
            double sum=0.;
            for(double xi:{.2113248654051871,.7886751345948129}) {
                double k=p.k+xi*(q.k-p.k);
                sum+=local_loss(load*dd/(k*k),std::abs(q.k-p.k),std::abs(q.b-p.b));
            }
            local_average+=.5*dd*sum;
            time_cdf.push_back(local_average);
        }
        for(size_t i=0;i<time_cdf.size();++i) {
            const double fraction=(x[i].D-x.front().D)/(x.back().D-x.front().D);
            time_cdf[i]=local_average>0. ? .25*fraction+.75*time_cdf[i]/local_average : fraction;
        }
        local_average/=(x.back().D-x.front().D);
    }
    double local()const{return local_average;}
    HitSum at(int i, double s, double u) const {
        const Point &p=x[i], &q=x[i+1];
        double dk=q.k-p.k,db=q.b-p.b,k=p.k+s*dk,b=p.b+s*db;
        double h=load*(q.D-p.D)/(k*k*std::abs(dk-db*u));
        Event e{i,s,k,u,-b+k*u,k*k*(1.-u)*(1.+u),h};
        HitSum hits;visit(0,e,hits);return hits;
    }
    std::array<double,8> estimate(int N,double phase_time,int order)const {
        std::array<double,8> avg{};
        if(order!=4 && order!=8)throw std::invalid_argument("Choose 4 or 8 angular nodes per panel");
        size_t i=0;
        for(int j=0;j<N;++j) {
            const double quantile=(j+phase_time)/N;
            while(i+2<x.size() && time_cdf[i+1]<=quantile)++i;
            const auto &p=x[i], &q=x[i+1]; double dd=q.D-p.D,dk=q.k-p.k,db=q.b-p.b;
            const double probability=time_cdf[i+1]-time_cdf[i];
            const double time_weight=dd/((x.back().D-x.front().D)*probability);
            double s=(quantile-time_cdf[i])/probability,k=p.k+s*dk,b=p.b+s*db;
            const double T=load*dd/(k*k);
            auto add_direction=[&](double u,double w) {
                double h=T/std::max(1e-300,std::abs(dk-db*u));
                Event e{static_cast<int>(i),s,k,u,-b+k*u,k*k*(1.-u)*(1.+u),h};
                HitSum hits;visit(0,e,hits);
                if(!hits.n || hits.H<.9999999999*h)throw std::runtime_error("Known birth root missing");
                double cap=std::min(1.,h);
                // For each label, sum(min(1,h_i)) - 1 after local clipping measures additional loss;
                // assign this excess to its visits in proportion to squared clipped loads.
                double extra=std::max(0.,hits.capped-1.)*cap*cap/(hits.capped2*h);
                avg[0]+=w;
                avg[1]+=w*std::max(0.,1.-1./h);
                avg[2]+=w*extra;
                avg[3]+=w*(hits.n>1);
                avg[4]+=w*std::max(0.,1.-1./hits.H);
                // Illustrative one-way rate prescriptions on globally visited labels.
                const double H=hits.H;
                const double correlated_loss=H<1.e-5 ? H*(.5-H/6.+H*H/24.) : 1.+std::expm1(-H)/H;
                avg[5]+=w*correlated_loss;
                avg[6]+=w*H/(1.+H);
                avg[7]+=w*extra_after_clearance(hits,18.)*cap*cap/(hits.capped2*h);
            };
            std::vector<double> edges{-1.,1.};
            if(db!=0.) {
                const double center=dk/db, width=T/std::abs(db);
                if(-1.<center && center<1.)edges.push_back(center);
                for(double gap=.25*width; gap<=2*(1+std::abs(center)); gap*=2.) {
                    if(gap==0)break;
                    for(double side: {-1.,1.}) {
                        const double boundary=center+side*gap;
                        if(-1.<boundary && boundary<1.)edges.push_back(boundary);
                    }
                }
            }
            std::sort(edges.begin(),edges.end());
            const bool high_order=order==8;
            const double g4[4]={-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526};
            const double w4[4]={.3478548451374539,.6521451548625461,.6521451548625461,.3478548451374539};
            const double g8[8]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,
                                 .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
            const double w8[8]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,
                                 .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
            const double *gauss=high_order?g8:g4, *wg=high_order?w8:w4;
            for(size_t edge=1;edge<edges.size();++edge) {
                double l=edges[edge-1],r=edges[edge];
                if(r<=l)continue;
                for(int n=0;n<(high_order?8:4);++n)
                    add_direction(.5*(l+r+(r-l)*gauss[n]),.25*time_weight*(r-l)*wg[n]);
            }
        }
        for(auto &v:avg)v/=N;
        return avg;
    }
};

std::vector<Point> cut(const std::vector<Point>& x,double lo,double hi) {
    if(x.front().D>lo || x.back().D<hi)throw std::runtime_error("Trace misses requested window");
    std::vector<Point> out;
    for(int i=0;i+1<static_cast<int>(x.size());++i) {
        const auto &p=x[i], &q=x[i+1];
        if(q.D<lo || p.D>=hi)continue;
        if(out.empty())out.push_back(mix(p,q,(lo-p.D)/(q.D-p.D)));
        if(q.D>=hi) {out.push_back(mix(p,q,(hi-p.D)/(q.D-p.D)));break;}
        if(q.t>out.back().t)out.push_back(q);
    }
    return out;
}


void validate() {
    History outward({{0.,1.,0.,0.},{1.,2.,0.,1.}},.8);
    if(outward.at(0,.3,.7).n!=1 || std::abs(outward.local())>1.e-14)
        throw std::runtime_error("Monotone shell control failed");
    History same({{0.,1.,0.,0.},{1.,2.,0.,1.},{2.,1.,0.,2.}},.8);
    auto two=same.at(0,.2,.3);
    if(two.n!=2 || std::abs(two.H-1.6/1.44)>1.e-12)
        throw std::runtime_error("Double radial passage control failed");
    const double no_reset=extra_after_clearance(two,0.);
    if(std::abs(no_reset-std::max(0.,two.capped-1.))>1.e-12 || extra_after_clearance(two,100.)>1.e-10)
        throw std::runtime_error("Clearance controls failed");
    double upper=std::sqrt(1.6), expected=upper-1.-(upper*upper*upper-1.)/4.8;
    auto a=same.estimate(1024,.37,4);
    if(std::abs(a[2]-expected)>1.e-4 || std::abs(a[1])>1.e-13)
        throw std::runtime_error("Integrated retraced shell control failed");
    History translated({{0.,1.,0.,0.},{1.,1.,4.,1.}},.6);
    auto pair=translated.at(0,.5,.2);
    if(pair.n!=2 || std::abs(pair.H-1.5)>1.e-12)
        throw std::runtime_error("Entering and leaving a translating shell failed");
}

struct Native {double t,dt,D0,D1,k0,k1,E0,E1,b0,b1;};
void add_piece(std::vector<Point>& p,const Native& row,double mu,double l,double r,int depth,
               int max_depth,double load) {
    const double dk=(row.k1-row.k0)/mu,db=(row.b1-row.b0)/mu;
    const double v0=row.E0*row.dt/mu,v1=row.E1*row.dt/mu;
    const double a=db-v0,b=v1-db;
    auto shift=[&](double f){return f*db+f*(f-1.)*(a+(b-a)*f);};
    auto deriv=[&](double f){return db+(2.*f-1.)*(a+(b-a)*f)+f*(f-1.)*(b-a);};
    auto point=[&](double f){return Point{row.t+row.dt*f,row.k0/mu+f*dk,row.b0/mu+shift(f),
                row.D0+f*(row.D1-row.D0)};};
    const double mid=.5*(l+r), k=row.k0/mu+mid*dk;
    const double slope=std::abs(deriv(mid));
    // Give the interpolated electric reversal fine angular-speed resolution
    // relative to either the expansion speed or the required per-angle source.
    const double required=load*(row.D1-row.D0)/(k*k);
    if(depth<max_depth && std::abs(deriv(l)-deriv(r))>
       .07*std::max({1.e-18,std::abs(dk),required,slope})) {
        add_piece(p,row,mu,l,mid,depth+1,max_depth,load);
        add_piece(p,row,mu,mid,r,depth+1,max_depth,load);
    } else {
        if(p.empty())p.push_back(point(l));
        p.push_back(point(r));
    }
}
int main(int argc,char **argv){
    validate();
    if(argc!=10)throw std::runtime_error("binary trace, mu, load, Dlo, Dhi, subdivision depth, N, time phase, angular quadrature order");
    std::ifstream f(argv[1],std::ios::binary);
    const double mu=std::stod(argv[2]),load=std::stod(argv[3]);
    const int depth=std::stoi(argv[6]),N=std::stoi(argv[7]);
    std::vector<Point> p;Native q{};
    while(f.read(reinterpret_cast<char*>(&q),sizeof(q)))add_piece(p,q,mu,0.,1.,0,depth,load);
    History h(cut(p,std::stod(argv[4]),std::stod(argv[5])),load);
    auto a=h.estimate(N,std::stod(argv[8]),std::stoi(argv[9]));
    std::cout<<std::setprecision(14)<<h.local();for(auto v:a)std::cout<<" "<<v;std::cout<<"\n";
}
