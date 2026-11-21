#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace hdboost {
struct Grid {
    int nt = 0, ny = 0;
    double t0, dt, dy;
    std::vector<double> data, scales;
};
inline Grid grid;
inline void load(const char* path) {
    std::ifstream input(path);
    Grid g;
    if (!(input >> g.nt >> g.ny >> g.t0 >> g.dt >> g.dy) ||
        g.nt < 4 || g.ny < 4 || g.dt <= 0.0 || g.dy <= 0.0)
        throw std::runtime_error("Invalid boosted drag table header");
    g.data.resize(g.nt * g.ny);
    g.scales.resize(g.nt);
    for (int i = 0; i < g.nt; ++i) {
        const double theta = std::pow(10., g.t0 + g.dt * i);
        g.scales[i] = std::sqrt(theta / (1. + theta));
        for (int j = 0; j < g.ny; ++j)
            if (!(input >> g.data[i * g.ny + j]) || !std::isfinite(g.data[i * g.ny + j]))
                throw std::runtime_error("Invalid boosted drag table entry");
    }
    grid = std::move(g);
}
inline double cubic(const double* arr, double x) {
    return arr[0] * (-(x-1.)*(x-2.)*(x-3.)/6.) + arr[1] * (x*(x-2.)*(x-3.)/2.)
         + arr[2] * (-x*(x-1.)*(x-3.)/2.) + arr[3] * (x*(x-1.)*(x-2.)/6.);
}
inline double cold_cross_section(double u) {
    const double e = std::hypot(1.,u), s=4*e*e, g=.30*(s+2.), d=s-8300.*8300.;
    const double b=g*u/e;
    return 3.*std::acos(-1.)*(g/e)*(g/e)/(d*d+b*b);
}
inline double cold_rate(double theta, double u) {
    if (theta == 0.) return 2.*u*cold_cross_section(u);
    const double thermal = std::sqrt(theta), z=u/thermal, root_pi=std::sqrt(std::acos(-1.));
    if (z < .08) {
        const double x=z*z;
        return 2.*thermal*cold_cross_section(u)*(8./3.+x*(8./15.+x*(-4./105.+x*(4./945.))))/root_pi;
    }
    if (z > 30.) return 2.*cold_cross_section(u)*(u+thermal/z-.25*thermal/(z*z*z));
    const double shape=(z+1./z-.25/(z*z*z))*std::erf(z)
                       +(1.+.5/(z*z))*std::exp(-z*z)/root_pi;
    return 2.*thermal*cold_cross_section(u)*shape;
}
inline double rate(double theta, double u) {
    u = std::abs(u);
    if (grid.nt == 0 || theta < 0. || !std::isfinite(theta) || !std::isfinite(u))
        throw std::runtime_error("Invalid boosted drag query or table not loaded");
    if (theta == 0.) return cold_rate(theta,u);
    const double it=(std::log10(theta)-grid.t0)/grid.dt;
    if (it < 0.) return cold_rate(theta,u);
    if (it > grid.nt-1+1.e-8) throw std::out_of_range("Boosted drag temperature outside table");
    const int i0=std::clamp(static_cast<int>(it)-1,0,grid.nt-4);
    double rows[4];
    for (int k=0;k<4;++k) {
        const int i=i0+k;
        const double iy=std::asinh(u/grid.scales[i])/grid.dy;
        if (iy > grid.ny-1+1.e-8) throw std::out_of_range("Boosted drag rapidity outside table");
        const int j=std::clamp(static_cast<int>(iy)-1,0,grid.ny-4);
        rows[k]=cubic(&grid.data[i*grid.ny+j],iy-j);
    }
    return std::exp(cubic(rows,it-i0));
}
}
inline void load_thermal_table(const char* path) { hdboost::load(path); }
inline double boosted_rate(double theta, double proper_velocity) {return hdboost::rate(theta,proper_velocity);}
