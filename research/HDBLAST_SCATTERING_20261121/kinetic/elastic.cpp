// Spatially homogeneous equal-mass Boltzmann collisions on a frozen wall time slice.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct Vec { double x, y, z; };
Vec operator+(Vec a, Vec b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
Vec operator-(Vec a, Vec b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
Vec operator*(double k, Vec a) { return {k*a.x, k*a.y, k*a.z}; }
double dot(Vec a, Vec b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec cross(Vec a, Vec b) { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
Vec unit(Vec a) {
    double length=std::sqrt(dot(a,a));
    if (!(length>0)) throw std::domain_error("Zero direction");
    return (1.0/length)*a;
}
double energy(Vec p) { return std::sqrt(1.0 + dot(p,p)); }

class Random {
    std::mt19937_64 engine;
public:
    explicit Random(uint64_t seed):engine(seed) {}
    double draw() { return (static_cast<double>(engine()>>12)+0.5)*0x1p-52; }
    size_t index(size_t n) { return static_cast<size_t>(draw()*n); }
};

class Group {
    // Binary sum tree permits a weighted draw and two weight changes per proposal.
    size_t leaves=1;
    std::vector<double> tree;
    Vec axis;
    double weight(Vec p) const {
        const double e=energy(p), along=dot(p,axis);
        if (along<=0 || dot(axis,axis)==0) return 1.0-along/e;
        const Vec perpendicular=cross(p,axis);
        return (1.0+dot(perpendicular,perpendicular))/(e*(e+along));
    }
public:
    std::vector<size_t> particles;
    Group(const std::vector<Vec>& p, int sign):axis{0,0,0} {
        Vec sum{0,0,0};
        for (size_t i=0;i<p.size();++i) {
            if (sign==0 || ((i%2==0)?1:-1)==sign) {
                particles.push_back(i);
                sum=sum+p[i];
            }
        }
        if (sign!=0 && dot(sum,sum)>0) axis=unit(sum);
        while (leaves<particles.size()) leaves*=2;
        tree.assign(2*leaves,0.0);
        for (size_t j=0;j<particles.size();++j) tree[leaves+j]=weight(p[particles[j]]);
        for (size_t j=leaves-1;j>0;--j) tree[j]=tree[2*j]+tree[2*j+1];
    }
    double eta(size_t index) const { return tree[leaves+index]; }
    double proposal_rate(size_t total) const {
        return 2.0*tree[1]*static_cast<double>(particles.size()-1)/(total-1);
    }
    size_t weighted_index(Random& rng) const {
        double mass=tree[1]*rng.draw();
        size_t i=1;
        while (i<leaves) {
            i*=2;
            if (mass>=tree[i]) { mass-=tree[i]; ++i; }
        }
        if (i-leaves>=particles.size()) throw std::runtime_error("Invalid weighted index");
        return i-leaves;
    }
    void replace(size_t j, Vec p) {
        size_t i=leaves+j;
        tree[i]=weight(p);
        while (i>1) { i/=2; tree[i]=tree[2*i]+tree[2*i+1]; }
    }
};

double minkowski(Vec a, Vec b, double ea, double eb) {
    const double ab=dot(a,b);
    if (ab<=0) return ea*eb-ab;
    const Vec normal=cross(a,b);
    return (1.0+dot(a,a)+dot(b,b)+dot(normal,normal))/(ea*eb+ab);
}

Vec random_direction(Random& rng, double cone_degrees, Vec incoming) {
    const double pi=std::acos(-1.0);
    const double cosine=1.0-rng.draw()*(1.0-std::cos(cone_degrees*pi/180.0));
    const double transverse=std::sqrt((1.0-cosine)*(1.0+cosine));
    const double phase=2.0*pi*rng.draw();
    Vec direction{transverse*std::cos(phase),transverse*std::sin(phase),cosine};
    if (cone_degrees==180) return direction;
    Vec axis=unit(incoming);
    Vec first=unit(cross(std::abs(axis.z)<0.75 ? Vec{0,0,1} : Vec{1,0,0},axis));
    return direction.x*first + direction.y*cross(axis,first) + direction.z*axis;
}

int main(int argc, char** argv) {
    if (argc!=7) throw std::runtime_error("Expected momenta, output-prefix, same/mixed/opposite, degrees, seed, checkpoints.");
    std::ifstream input(argv[1],std::ios::binary|std::ios::ate);
    auto bytes=input.tellg();
    if (bytes<=0 || bytes%sizeof(Vec)!=0) throw std::runtime_error("Bad input");
    std::vector<Vec> p(static_cast<size_t>(bytes)/sizeof(Vec));
    input.seekg(0); input.read(reinterpret_cast<char*>(p.data()),bytes);
    const std::string prefix=argv[2], kind=argv[3];
    const double cone=std::stod(argv[4]);
    if (p.size()%4 || !(cone>0 && cone<=180)) throw std::runtime_error("Invalid trial");
    Random rng(std::stoull(argv[5]));
    std::ifstream checkpoints(argv[6]);
    std::vector<Group> groups;
    const bool opposite_only=kind=="opposite";
    if (kind=="mixed" || opposite_only) groups.emplace_back(p,0);
    else if (kind=="same") { groups.emplace_back(p,1); groups.emplace_back(p,-1); }
    else throw std::runtime_error("Invalid channels");
    double exposure=0, max_pair_error=0, max_majorant=0;
    uint64_t accepted=0, proposed=0, mixed=0;
    int target;
    while (checkpoints>>target) {
        while (accepted<static_cast<uint64_t>(target)*p.size()/32) {
            const double first=groups[0].proposal_rate(p.size());
            const double group_rate=first+(groups.size()==2 ? groups[1].proposal_rate(p.size()):0.0);
            // Each unlike pair draws at rate 2/(N-1), using v_M <= 2.
            const double total=opposite_only ? 2.0*(p.size()/2)*(p.size()/2)/(p.size()-1) : group_rate;
            exposure-=std::log(rng.draw())/total;
            Group& group=groups[(groups.size()==1 || rng.draw()*total<first)?0:1];
            size_t i, j;
            double bound;
            if (opposite_only) {
                i=2*rng.index(p.size()/2);
                j=2*rng.index(p.size()/2)+1;
                bound=2.0;
            } else {
                i=group.weighted_index(rng);
                j=rng.index(group.particles.size()-1);
                if (j>=i) ++j;
                bound=2.0*(group.eta(i)+group.eta(j));
            }
            const size_t ia=group.particles[i], ib=group.particles[j];
            Vec a=p[ia], b=p[ib];
            const double ea=energy(a), eb=energy(b);
            const double product=std::max(1.0,minkowski(a,b,ea,eb));
            const double flux=std::sqrt((product-1.0)*(product+1.0))/(ea*eb);
            const double probability=flux/bound;
            max_majorant=std::max(max_majorant,probability);
            if (probability>1.0+5e-9) throw std::runtime_error("Invalid pair majorant");
            ++proposed;
            if (rng.draw()>=probability) continue;
            const Vec total_p=a+b;
            const double total_e=ea+eb, root_s=std::sqrt(2.0*(1.0+product));
            const double k=std::sqrt(0.5*(product-1.0));
            const Vec in_cm=a-((ea+0.5*root_s)/(total_e+root_s))*total_p;
            const Vec outgoing=random_direction(rng,cone,in_cm);
            const double parallel=dot(outgoing,total_p)/(root_s*(total_e+root_s));
            p[ia]=0.5*total_p+k*(outgoing+parallel*total_p);
            p[ib]=total_p-p[ia];
            max_pair_error=std::max(max_pair_error,std::abs(energy(p[ia])+energy(p[ib])-total_e)/total_e);
            group.replace(i,p[ia]); group.replace(j,p[ib]);
            ++accepted;
            if (ia%2!=ib%2) ++mixed;
        }
        std::ofstream out(prefix+std::to_string(target), std::ios::binary);
        out.write(reinterpret_cast<const char*>(p.data()),bytes);
        if (!out) throw std::runtime_error("Cannot write snapshot");
        std::cout<<std::setprecision(17)<<target<<' '<<exposure<<' '<<accepted<<' '<<proposed<<' '<<mixed
                 <<' '<<max_pair_error<<' '<<max_majorant<<'\n';
    }
    if (accepted==0) throw std::runtime_error("No positive checkpoint");
}
