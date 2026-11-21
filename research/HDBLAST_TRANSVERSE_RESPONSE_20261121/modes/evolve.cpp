// Linearized transverse Maxwell response of continuously injected cold charge-conjugate cohorts.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using Matrix = std::vector<double>;
constexpr double Q=1.e7, K=2.54;
struct Point {double a, guide, impulse;};
struct Cohort {double w, impulse_birth, source_scale; int step, group;};
struct Moment {double x=0.,y=0.,z=0.,d=0.;};

double operator_norm(const Matrix& b,const Matrix& e) {
    double u=0.,v=0.,m=0.;
    for (int j=0;j<4;++j) {u+=b[j]*b[j];v+=e[j]*e[j];m+=b[j]*e[j];}
    return std::sqrt(.5*(u+v+std::hypot(u-v,2*m)));
}
void rotate(double& b,double& e,double c,double s) {
    const double old=b;b=c*b+s*e;e=c*e-s*old;
}
int main(int argc,char** argv) {
    if(argc!=4) throw std::runtime_error("input k_comoving guide_scale");
    const double k=std::stod(argv[2]),guide_scale=std::stod(argv[3]);
    std::ifstream f(argv[1],std::ios::binary);
    int n,steps,groups; double dt;
    auto read=[&](auto& v){f.read(reinterpret_cast<char*>(&v),sizeof(v));};
    read(n);read(steps);read(groups);read(dt);
    std::vector<Point> bg(2*steps+1);
    for(auto& z:bg)read(z);
    std::vector<Cohort> births(n);
    for(auto& b:births){read(b.w);read(b.impulse_birth);read(b.source_scale);read(b.step);read(b.group);}
    if(!f || groups<0 || n<1 || steps<1)throw std::runtime_error("Bad input");
    const int cols=5+groups;
    // Four independent electromagnetic seeds, coherent birth modulation, then covariance factors.
    enum{BY,EZ,BZ,EY};
    using Block=std::array<Matrix,4>;
    Block fields,half_fields,mid_fields,vacuum;
    for(int r=0;r<4;r++){fields[r]=Matrix(cols);fields[r][r]=1.;}
    double comoving_phase=0., peak_perp_field=0.;
    Matrix old_jy(cols),old_jz(cols),jy(cols),jz(cols);
    std::vector<Moment> p(n*cols);
    using Seg=std::array<double,5>;
    Seg accum{}; // time integral of redshift-weighted force norm per initial field unit.
    std::vector<Seg> at_birth(n);
    int active=0;
    const double a0=bg[0].a;
    for(int l=0;l<steps;++l) {
        const auto& L=bg[2*l]; const auto& M=bg[2*l+1];const auto& R=bg[2*l+2];
        while(active<n && births[active].step==l) {
            const auto& b=births[active];at_birth[active]=accum;
            // d is fractional number modulation divided by Q. P_i=a delta p_i/(a0 Q).
            p[active*cols+4].d=1./Q;
            if(groups) p[active*cols+5+b.group].d=b.source_scale/Q;
            ++active;
        }
        const double ar=a0/M.a, physical_k=k/M.a;
        const double angle=.5*dt*physical_k,c=std::cos(angle),s=std::sin(angle);
        half_fields=fields;
        for(int j=0;j<cols;j++){
            rotate(half_fields[BY][j],half_fields[EZ][j],c,s);
            rotate(half_fields[BZ][j],half_fields[EY][j],c,s);
        }
        mid_fields=half_fields;
        for(int j=0;j<cols;j++){
            mid_fields[EY][j]-=.5*dt*old_jy[j];
            mid_fields[EZ][j]-=.5*dt*old_jz[j];
        }
        std::fill(jy.begin(),jy.end(),0.);std::fill(jz.begin(),jz.end(),0.);
        for(int i=0;i<active;i++) {
            const auto& b=births[i];
            const double zl=(L.impulse-b.impulse_birth)/L.a,zr=(R.impulse-b.impulse_birth)/R.a;
            const double gl=std::hypot(1.,zl),gr=std::hypot(1.,zr),vl=zl/gl,vr=zr/gr;
            const double dz=zr-zl;
            double invg,avg3,first3;
            if(std::abs(dz)>1.e-7*std::max(1.,std::abs(zl))) {
                invg=(std::asinh(zr)-std::asinh(zl))/dz;
                // Compute mean gamma^-3 from a rationalized difference of velocities.
                if(zr*zl>0.) avg3=(zr+zl)/((zr*gl+zl*gr)*gl*gr);
                else avg3=(vr-vl)/dz;
                first3=((zr+zl)/(gr*gl*(gl+gr))-zl*avg3)/dz;
            } else {
                const double gm=std::hypot(1.,.5*(zl+zr));
                invg=1./gm;avg3=invg*invg*invg;first3=.5*avg3;
            }
            const double v=(zr+zl)/(gr+gl);
            const double hangle=.5*dt*Q*guide_scale*M.guide*invg;
            const double gc=std::cos(hangle),gs=std::sin(hangle);
            const double halfdt=.5*dt;
            const double gsin=std::abs(hangle)>1.e-6 ? halfdt*gs/hangle : halfdt*(1.-hangle*hangle/6.);
            const double gcos=std::abs(hangle)>1.e-6 ? halfdt*2.*std::pow(std::sin(.5*hangle),2)/hangle : halfdt*.5*hangle;
            for(int j=0;j<cols;j++){
                auto& z=p[i*cols+j];
                const double fx=-ar*v*mid_fields[BY][j],fy=ar*mid_fields[EY][j],fz=ar*mid_fields[EZ][j];
                const double addx=gsin*fx+gcos*fy,addy=gsin*fy-gcos*fx;
                const double xm=gc*z.x+gs*z.y+addx,ym=gc*z.y-gs*z.x+addy;
                const double dd=dt*physical_k*ar*invg*xm;
                jy[j]+=b.w*(ar*invg*ym);
                jz[j]+=b.w*(ar*(z.z*avg3+dt*fz*first3)+v*(z.d+.5*dd));
                z.x=gc*xm+gs*ym+addx;z.y=gc*ym-gs*xm+addy;
                z.z+=dt*fz;z.d+=dd;
            }
        }
        fields=half_fields;
        const double source=Q*Q/(K*M.a*a0*a0);
        for(int j=0;j<cols;j++){
            old_jy[j]=source*jy[j];old_jz[j]=source*jz[j];
            fields[EY][j]-=dt*old_jy[j];fields[EZ][j]-=dt*old_jz[j];
            rotate(fields[BY][j],fields[EZ][j],c,s);rotate(fields[BZ][j],fields[EY][j],c,s);
        }
        const double midphase=comoving_phase+angle; comoving_phase+=2*angle;
        const double peak=operator_norm(mid_fields[BY],mid_fields[EY]);
        peak_perp_field=std::max(peak_perp_field,ar*ar*peak);
        Matrix induced_b=mid_fields[BY],induced_e=mid_fields[EY];
        induced_b[0]-=std::cos(midphase); induced_b[1]-=std::sin(midphase);
        induced_e[2]+=std::sin(midphase); induced_e[3]-=std::cos(midphase);
        double noise=0.;
        for(int j=5;j<cols;j++)noise+=mid_fields[BY][j]*mid_fields[BY][j]+mid_fields[EY][j]*mid_fields[EY][j];
        Seg integrand{peak,operator_norm(induced_b,induced_e),
            std::hypot(mid_fields[BY][0],mid_fields[EY][0]),
            std::hypot(mid_fields[BY][4],mid_fields[EY][4])/std::sqrt(2.),std::sqrt(noise)};
        for(size_t j=0;j<accum.size();++j)accum[j]+=dt*M.a*ar*ar*integrand[j];
    }
    const double af=bg.back().a, ar=a0/af;
    double total=0.,shotnet=0.,netcoherent=0.; Seg bound{};double gram[4][4]{};
    for(int i=0;i<active;++i){
        const double w=births[i].w;total+=w;
        for(size_t j=0;j<bound.size();j++)bound[j]+=w*std::pow((accum[j]-at_birth[i][j])*Q/af,2);
        for(int j=0;j<cols;++j){
            auto const& u=p[i*cols+j];
            if(j<4)for(int h=0;h<4;++h){auto const& v=p[i*cols+h];gram[j][h]+=w*(u.x*v.x+u.y*v.y);}
            else if(j==4)netcoherent+=w*(u.x*u.x+u.y*u.y)/2.;
            else shotnet+=w*(u.x*u.x+u.y*u.y);
        }
    }
    if(active!=n || total<=0.)throw std::runtime_error("Birth schedule incomplete");
    std::cout<<std::setprecision(17);
    for(auto x:bound)std::cout<<std::sqrt(x/total)<<' ';
    std::cout<<Q*ar*std::sqrt(netcoherent/total)<<' '<<Q*ar*std::sqrt(shotnet/total)<<' '<<peak_perp_field<<' '<<comoving_phase;
    for(int i=0;i<4;i++)for(int j=0;j<4;j++)std::cout<<' '<<gram[i][j]*Q*Q*ar*ar/total;
    for(int i=0;i<4;i++)for(int j=0;j<4;j++)std::cout<<' '<<fields[i][j]*ar*ar;
    std::cout<<'\n';
}
