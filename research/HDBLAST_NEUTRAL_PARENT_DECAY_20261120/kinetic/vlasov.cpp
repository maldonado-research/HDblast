// Collisionless charged daughter shells on an expanding homogeneous gauge triplet.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Node { double t,a,f,R,scalar; };
struct Shell { double weight, potential, momentum; };
struct Moments { double current=0, gamma=0, inverse=0; };

Moments moments(const std::vector<Shell>& shells, double a, double A, double Q) {
    Moments sum;
    for(const auto& c:shells) {
        const double x=Q*(A-c.potential)/a,s=c.momentum/a;
        const double plus=std::sqrt(1+(x+s)*(x+s)),minus=std::sqrt(1+(x-s)*(x-s));
        const double scale=1/(plus+minus),up=plus*scale,down=minus*scale;
        const double vel=(2./3)*((x+s)/plus*(1-down*down)+(x-s)/minus*(1-up*up));
        sum.current+=c.weight*vel;
        sum.gamma+=c.weight*(2./3)*(plus+minus-plus*down);
        sum.inverse+=c.weight*2*scale;
    }
    return sum;
}

int main(int argc,char** argv) {
    if(argc!=4) throw std::runtime_error("Usage: solver input snapshots shells");
    std::ifstream input(argv[1]);
    int bins,substeps,interval;
    double Q,K,D,A,P;
    if(!(input>>bins>>substeps>>interval>>Q>>K>>D>>A>>P) || substeps%2!=0)
        throw std::runtime_error("Invalid parameters or odd substeps");
    std::vector<Node> grid(bins*substeps+1);
    for(auto& n:grid)
        if(!(input>>n.t>>n.a>>n.f>>n.R>>n.scalar))
            throw std::runtime_error("Truncated background");
    std::ofstream trace(argv[2]),packets(argv[3]);
    trace<<std::setprecision(16);
    packets<<std::setprecision(16);
    trace<<"t a A P f field daughter parent scalar current fluxwork volume_work net_residual field_residual birthkin forward_work returned_work\n";
    std::vector<Shell> shells;
    const double E0=.5*K*(P*P+4*A*A);
    const double total0=E0+D*grid[0].a*grid[0].R;
    double f=0,volume_work=0,fluxwork=0,birthkin=0,forward=0,returned=0;
    const int last=grid.size()-1;
    for(int i=0;i<=last;++i) {
        const auto& n=grid[i];
        if(i%substeps==substeps/2) {
            const auto& left=grid[i-substeps/2];
            const auto& right=grid[i+substeps/2];
            const double df=right.f-left.f;
            if(df < -1.e-11) throw std::runtime_error("Negative birth weight");
            if(df>0) {
                const double p=std::sqrt(std::max(0.,n.R*n.R-1));
                shells.push_back({df,A,n.a*p});
                f+=df;
                birthkin+=D*df*(std::sqrt(1.+p*p)-1);
            }
        }
        if(i%interval==0 || i==last) {
            const auto m=moments(shells,n.a,A,Q);
            const double field=.5*K*(P*P+4*A*A), daughter=D*n.a*m.gamma;
            const double parent=D*n.a*(1-f)*n.R;
            trace<<n.t<<' '<<n.a<<' '<<A<<' '<<P<<' '<<f<<' '<<field<<' '<<daughter<<' '<<parent<<' '
                 <<std::pow(n.a,4)*n.scalar<<' '<<m.current<<' '<<fluxwork<<' '<<volume_work<<' '
                 <<field+daughter+parent-total0-volume_work<<' '<<field+fluxwork-E0<<' '<<birthkin<<' '<<forward<<' '<<returned<<'\n';
        }
        if(i==last) break;
        const auto& z=grid[i+1];
        const double a=(n.a+z.a)/2., h=(z.t-n.t)/a;
        const double mid=A+h*P/2.;
        const auto m=moments(shells,a,mid,Q);
        const double nextP=P-h*(4*mid+(D*Q/K)*m.current);
        const double nextA=mid+h*nextP/2.;
        const double work=D*Q*m.current*(nextA-A);
        fluxwork+=work;
        forward+=std::max(0.,work);
        returned+=std::max(0.,-work);
        volume_work+=D*((m.inverse+(1-f)*(z.R+n.R)/2.)*(z.a-n.a)+a*(1-f)*(z.R-n.R));
        A=nextA;P=nextP;
    }
    for(const auto& s:shells)
        packets<<s.weight<<' '<<s.potential<<' '<<s.momentum<<'\n';
}
