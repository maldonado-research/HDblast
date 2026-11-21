// Linear k=2 Maxwell propagation with prescribed shell E and absorbing outer tube.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using std::array;
using std::vector;

constexpr int P=4;
const array<double,P> x={-1.,-std::sqrt(.2),std::sqrt(.2),1.};
const array<double,P> w={1./6.,5./6.,5./6.,1./6.};

struct Geometry {double r,rx,u,n,nx;};
struct Wave {double i=0,o=0,v=0;};
using State=vector<Wave>;
struct Budget {double wall=0,tube=0,geometry=0,dissipation=0;};

vector<double> load(const char *path) {
    std::ifstream f(path,std::ios::binary|std::ios::ate);
    if(!f) throw std::runtime_error("Missing input");
    const size_t n=f.tellg();
    if(n%sizeof(double)) throw std::runtime_error("Truncated input");
    vector<double> a(n/sizeof(double));
    f.seekg(0); f.read(reinterpret_cast<char*>(a.data()),n);
    if(!f) throw std::runtime_error("Failed to read");
    return a;
}

struct Response {
    const vector<double> geom,src;
    int base,cells,ng,ns,gr=0,sr=0,width;
    double dx,duration_step;
    array<array<double,P>,P> D{};
    vector<Geometry> geo;
    vector<array<double,P>> basis;
    State y,k1,k2,k3,k4,tmp;
    vector<Budget> totals;

    Response(const char* gf,const char* sf,int n,int refine,double depth,double cfl):
        geom(load(gf)),src(load(sf)),base(n),cells(n*refine),width(1+5*n*P),
        dx(depth/cells),duration_step(cfl*dx),geo(cells*P),basis(cells*P),
        y(3*cells*P),k1(y.size()),k2(y.size()),k3(y.size()),k4(y.size()),tmp(y.size()),totals(3) {
        if(cfl<=0||cfl>.15||refine<1||geom.size()%width||src.size()%4)
            throw std::runtime_error("Unexpected integration data");
        ng=geom.size()/width;ns=src.size()/4;
        if(ng<2||ns<2||src[0]<geom[0]||src[4*(ns-1)]>geom[width*(ng-1)]+1e-10)
            throw std::runtime_error("Uncovered source times");
        for(int k=0;k<P;++k)for(int j=0;j<P;++j){
            if(k==j){for(int a=0;a<P;++a)if(a!=k)D[k][j]+=1./(x[k]-x[a]);}
            else{D[k][j]=1./(x[j]-x[k]);for(int a=0;a<P;++a)if(a!=k&&a!=j)D[k][j]*=(x[k]-x[a])/(x[j]-x[a]);}
        }
        for(int z=0;z<cells*P;++z) {
            const double point=(2.*((z/P)%refine)+1.+x[z%P])/refine-1.;
            for(int j=0;j<P;++j) {
                double b=1.;
                for(int k=0;k<P;++k)if(j!=k)b*=(point-x[k])/(x[j]-x[k]);
                basis[z][j]=b;
            }
        }
    }

    void interpolate(double t, array<double,3>& source) {
        while(gr+2<ng && geom[(gr+1)*width]<t)++gr;
        while(sr+2<ns && src[(sr+1)*4]<t)++sr;
        // RK switching from k4 to the next k1 always keeps time nondecreasing.
        const double f=(t-geom[gr*width])/(geom[(gr+1)*width]-geom[gr*width]);
        const double h=(t-src[sr*4])/(src[(sr+1)*4]-src[sr*4]);
        for(int m=0;m<3;++m)source[m]=(t<src[0])?0.:src[sr*4+1+m]+h*(src[(sr+1)*4+1+m]-src[sr*4+1+m]);
        const int refine=cells/base;
        for(int z=0;z<cells*P;++z) {
            array<double,5> q{};
            for(int j=0;j<P;++j) {
                const int k=gr*width+1+5*((z/P)/refine*P+j);
                for(int a=0;a<5;++a) q[a]+=basis[z][j]*(geom[k+a]+f*(geom[k+width+a]-geom[k+a]));
            }
            geo[z]={q[0],q[1],q[2],q[3],q[4]};
        }
    }

    void rhs(double t,const State& b,State& d,array<Budget,3>& flux){
        array<double,3> ewall;interpolate(t,ewall);
        for(int z=0;z<cells*P;++z){
            const int k=z%P;const auto& g=geo[z];
            const double r=g.r, n=g.n, r3=r*r*r;
            const double expansion=n*g.u/r, grad=n*g.rx/r;
            const double measure=.5*dx*w[k];
            for(int m=0;m<3;++m) {
                const int q=3*z+m;const auto value=b[q];auto& energy=flux[m];
                const double e=value.i+value.o, h=value.i-value.o;
                double di=0.,doo=0.;
                for(int j=0;j<P;++j){
                    const auto other=b[q+3*(j-k)];
                    di+=D[k][j]*other.i;doo+=D[k][j]*other.o;
                }
                if(k==P-1) {
                    const double incoming=z+1==cells*P?0.:b[q+3].i;
                    di+=(incoming-value.i)/w[k];
                    energy.dissipation+=n*r3*std::pow(incoming-value.i,2);
                }
                if(k==0) {
                    const double outgoing=z==0?ewall[m]-value.i:b[q-3].o;
                    doo-=(outgoing-value.o)/w[k];
                    energy.dissipation+=n*r3*std::pow(outgoing-value.o,2);
                }
                const double se=(g.nx+2.*grad)*h-2.*expansion*e-4.*n*value.v/r3;
                const double sh=(g.nx+grad)*e-expansion*h;
                d[q]={2.*n*di/dx+.5*(se+sh),-2.*n*doo/dx+.5*(se-sh),n*r*e};
                energy.geometry+=measure*(r3*g.nx*e*h+.5*n*g.u*r*r*(h*h-e*e)-2.*n*g.u*value.v*value.v/(r*r));
                if(z==0)energy.wall=n*r3*(std::pow(ewall[m]-value.i,2)-value.i*value.i);
                if(z+1==cells*P)energy.tube=n*r3*value.o*value.o;
            }
        }
    }

    void stage(State& out,const State& d,double h){
        for(size_t z=0;z<y.size();++z)out[z]={y[z].i+h*d[z].i,y[z].o+h*d[z].o,y[z].v+h*d[z].v};
    }
    void advance(double t,double h) {
        array<Budget,3> f1{},f2{},f3{},f4{};
        rhs(t,y,k1,f1);stage(tmp,k1,.5*h);
        rhs(t+.5*h,tmp,k2,f2);stage(tmp,k2,.5*h);
        rhs(t+.5*h,tmp,k3,f3);stage(tmp,k3,h);
        rhs(t+h,tmp,k4,f4);
        for(size_t j=0;j<y.size();++j){
            auto& q=y[j];const auto a=k1[j],b=k2[j],c=k3[j],d=k4[j];
            q.i+=h/6*(a.i+2*b.i+2*c.i+d.i);
            q.o+=h/6*(a.o+2*b.o+2*c.o+d.o);
            q.v+=h/6*(a.v+2*b.v+2*c.v+d.v);
        }
        for(int m=0;m<3;++m) {
            auto accumulate=[&](double Budget::*member) {
                totals[m].*member+=h/6*(f1[m].*member+2*(f2[m].*member)+2*(f3[m].*member)+f4[m].*member);
            };
            accumulate(&Budget::wall);accumulate(&Budget::tube);accumulate(&Budget::geometry);accumulate(&Budget::dissipation);
        }
    }
    void run(const char* output){
        std::ofstream f(output,std::ios::binary);
        double t=src[0];
        for(int row=0;row<ns;++row){
            const double end=src[row*4];int steps=std::max(0,int(std::ceil((end-t)/duration_step)));
            const double start=t;
            for(int j=0;j<steps;++j){
                const double next=start+(end-start)*((j+1.)/steps);
                advance(t,next-t);t=next;
            }
            array<double,3> ewall;interpolate(t,ewall);
            for(int m=0;m<3;++m){
                double E=0.;
                for(int z=0;z<cells*P;++z){
                    const auto u=y[3*z+m];const double r=geo[z].r;
                    E+=.5*dx*w[z%P]*(r*r*r*(u.i*u.i+u.o*u.o)+2.*u.v*u.v/r);
                }
                const auto& b=totals[m];
                const array<double,7> result={y[m].i,y[(cells*P-1)*3+m].o,E,b.wall,b.tube,b.geometry,b.dissipation};
                f.write(reinterpret_cast<const char*>(result.data()),sizeof(result));
            }
        }
        if(!f)throw std::runtime_error("Failed to write response");
    }
};

int main(int argc,char** argv) {
    if(argc!=8)throw std::runtime_error("geometry, drive, output, cells, subdivision, depth, CFL");
    Response problem(argv[1],argv[2],std::stoi(argv[4]),std::stoi(argv[5]),std::stod(argv[6]),std::stod(argv[7]));
    problem.run(argv[3]);
}
