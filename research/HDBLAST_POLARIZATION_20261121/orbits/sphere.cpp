// Test particles on S^3 in two coexact electromagnetic harmonics on a prescribed mean wall.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using Complex=std::complex<double>;
using Point=std::array<Complex,2>;
using Row=std::array<double,8>;
struct Vector { double x,y,z; };
Vector operator+(Vector a,Vector b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vector operator*(double c,Vector a){return {c*a.x,c*a.y,c*a.z};}
Vector cross(Vector a,Vector b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double square(Vector a){return a.x*a.x+a.y*a.y+a.z*a.z;}

std::vector<Row> read(const std::string& path) {
    std::ifstream f(path,std::ios::binary);
    if(!f)throw std::runtime_error("open input");
    std::vector<Row> rows; Row r;
    while(f.read(reinterpret_cast<char*>(r.data()),sizeof(r))) rows.push_back(r);
    if(!f.eof() || f.gcount() || rows.empty())throw std::runtime_error("truncated input");
    return rows;
}

Complex harmonic(Complex z,int k){
    Complex power(1.,0.);
    for(int i=2;i<k;++i)power*=z;
    return std::sqrt(double(k-1))*std::conj(power);
}

struct Polynomial {
    int degree;
    std::vector<Complex> coefficient;
    explicit Polynomial(int k):degree(k-2),coefficient(k-1) {
        // |coefficient_j|² = choose(n,j) gives unit Haar RMS for the whole sum.
        const double pi=std::acos(-1.);
        for(int j=0;j<=degree;++j) {
            double mag=std::exp(.5*(std::lgamma(degree+1.)-std::lgamma(j+1.)-std::lgamma(degree-j+1.)));
            coefficient[j]=std::polar(mag,pi*std::sqrt(2.)*j*j);
        }
    }
    Complex eval(const Point& z) const {
        bool reverse=std::abs(z[0])>std::abs(z[1]);
        Complex lead=reverse?z[0]:z[1],r=(reverse?z[1]:z[0])/lead;
        Complex h=coefficient[reverse?0:degree];
        for(int j=1;j<=degree;++j) h=h*r+coefficient[reverse?j:degree-j];
        return std::conj(h*std::pow(lead,degree));
    }
};

void drift(Point& q,Vector u,double a,double duration){
    const double magnitude=std::sqrt(square(u)), energy=std::hypot(a,magnitude);
    const double angle=duration*magnitude/energy;
    const double scale=(magnitude==0.)?duration/energy:std::sin(angle)/magnitude;
    const Complex horizontal(u.x,u.y), vertical(0.,-u.z);
    const Point tangent={-horizontal*std::conj(q[1])+vertical*q[0],horizontal*std::conj(q[0])+vertical*q[1]};
    const double c=std::cos(angle);
    for(int j=0;j<2;++j)q[j]=c*q[j]+scale*tangent[j];
}

Vector displacement(const Point& q,const Point& to){
    const Complex horiz=q[0]*to[1]-q[1]*to[0];
    const Complex ip=std::conj(q[0])*to[0]+std::conj(q[1])*to[1];
    const Vector chord={horiz.real(),horiz.imag(),-ip.imag()};
    const double sn=std::sqrt(square(chord));
    const double ratio=(sn>0.)?std::atan2(sn,ip.real())/sn:1.;
    return ratio*chord;
}

Point shift(Point q,Vector delta){
    const double length=std::sqrt(square(delta));
    if(length>0.)drift(q,delta,0.,length);
    return q;
}

struct Model {
    const std::vector<Row>& fields;
    int k1,k2;
    double a1,a2;
    Complex phase;
    bool distributed;
    Polynomial poly;

    std::pair<Point,Vector> evolve(Point q,Vector u,double born,int sign)const{
        const int steps=int(fields.size()-1);
        const double t0=fields.front()[0],step=(fields.back()[0]-t0)/steps;
        const int first=std::clamp(int((born-t0)/step),0,steps-1);
        for(int j=first;j<steps;++j){
            const auto& left=fields[j];const auto& right=fields[j+1];
            const double lo=std::max(left[0],born),dt=right[0]-lo;
            if(dt<=0.)continue;
            const double s=(lo-left[0]+.5*dt)/(right[0]-left[0]);
            auto sample=[&](int col){return left[col]+s*(right[col]-left[col]);};
            const double a=sample(1),deta=dt/a,half=.5*deta;
            drift(q,u,a,half);
            const Complex y1=a1*harmonic(q[0],k1);
            const Complex y2=a2*phase*(distributed?poly.eval(q):harmonic(q[0],k2));
            const Complex electric=sample(4)*y1+sample(6)*y2;
            const Complex magnetic=sample(5)*y1+sample(7)*y2;
            const Vector E={electric.real(),electric.imag(),sample(2)};
            const Vector B={magnetic.real(),magnetic.imag(),sample(3)};
            const double h=.5e7*sign*a*a*deta;
            const Vector p=u+h*E;
            const Vector rot=(h/std::sqrt(a*a+square(p)))*B;
            const Vector w=p+cross(p,rot);
            u=p+cross(w,(2./(1.+square(rot)))*rot)+h*E;
            drift(q,u,a,half);
        }

        return {q,u};
    }
};

int main(int argc,char** argv){
    if(argc!=10 && argc!=11)throw std::runtime_error("fields seeds k1 k2 amp1 amp2 phase2 shape2 output [conditional]");
    const auto fields=read(argv[1]),seeds=read(argv[2]);
    const int k1=std::stoi(argv[3]),k2=std::stoi(argv[4]);
    if(k1<2||k2<2||fields.size()<2)throw std::runtime_error("input settings");
    const Model model{fields,k1,k2,std::stod(argv[5]),std::stod(argv[6]),std::polar(1.,std::stod(argv[7])),
                      std::string(argv[8])=="distributed",Polynomial(k2)};
    const bool conditional=argc==11;
    std::ofstream out(argv[9],std::ios::binary);
    if(!out)throw std::runtime_error("open output");
    for(const auto& birth:seeds)for(int sign:{1,-1}){
        const Point location={Complex(birth[1],birth[2]),Complex(birth[3],birth[4])};
        const Vector initial=double(sign)*Vector{birth[5],birth[6],birth[7]};
        Point start=location;
        auto result=model.evolve(start,initial,birth[0],sign);
        double residual=0.,det=1.,jacobian_change=0.;
        if(conditional){
            // Invert the map from uniform birth volume to final observation volume.
            for(int iter=0;iter<18;++iter){
                const Vector err=displacement(result.first,location);
                residual=std::sqrt(square(err));
                if(residual<2.e-12)break;
                start=shift(start,err);
                result=model.evolve(start,initial,birth[0],sign);
            }
            if(residual>2.e-10)throw std::runtime_error("conditional position failed to invert");
            // First-order geodesic differencing in three orthonormal tangent directions.
            const double eps=2.e-6;
            std::array<Vector,3> j;
            for(int axis=0;axis<3;++axis){
                Vector e={0.,0.,0.};
                if(axis==0)e.x=1.;
                if(axis==1)e.y=1.;
                if(axis==2)e.z=1.;
                Point plus=model.evolve(shift(start,eps*e),initial,birth[0],sign).first;
                j[axis]=(1./eps)*displacement(result.first,plus);
                jacobian_change+=square(j[axis]+(-1.)*e);
            }
            const auto normal=cross(j[1],j[2]);
            det=j[0].x*normal.x+j[0].y*normal.y+j[0].z*normal.z;
            if(!(det>.2 && det<5.))throw std::runtime_error("conditional spatial Jacobian");
        }
        const auto& q=result.first;
        const Vector u=(1./fields.back()[1])*result.second;
        const Row row{double(sign),q[0].real(),q[0].imag(),q[1].real(),q[1].imag(),u.x,u.y,u.z};
        out.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
        if(conditional){
            const std::array<double,3> extra{det,residual,std::sqrt(jacobian_change)};
            out.write(reinterpret_cast<const char*>(extra.data()),sizeof(extra));
        }
    }
    if(!out)throw std::runtime_error("write failure");
}
