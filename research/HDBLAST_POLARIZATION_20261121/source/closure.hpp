// Angular zero-mode projection: the axial mode and each added source are orthonormal.
struct Extra {
    double e2=0., b2=0., i2=0., o2=0.;
};

Extra extra_wall(const vector<double>& y, double a, double gauge) {
    Extra v;
    if (!backreact) return v;
    for (int m=0;m<unforced;m++) {
        const double e=y[PROBE_MOM+m]/gauge, b=probe_lambda[m]*y[PROBE_C+m]/a;
        const double inc=y[COUNT+3+m], out=e-inc;
        v.e2+=e*e; v.b2+=b*b; v.i2+=inc*inc; v.o2+=out*out;
    }
    return v;
}

// For each Maxwell mode, <e²>, <z²>, <b²> and <ez> enter the spherically projected stress.
Geo radial(Geo g, Fields w) {
    double e=w.i+w.o, z=w.i-w.o, b=2*w.a/(g.r*g.r);
    double e2=e*e,z2=z*z,b2=b*b,ez=e*z;
    if (backreact) for (int m=0;m<unforced;m++) {
        e=w.more[m][0]+w.more[m][1]; z=w.more[m][0]-w.more[m][1];
        b=probe_lambda[m]*w.more[m][2]/(g.r*g.r);
        e2+=e*e; z2+=z*z; b2+=b*b; ez+=e*z;
    }
    const double mr=g.r*(1+g.r*g.r+g.u*g.u-g.rx*g.rx);
    return {g.rx, g.r+mr/(g.r*g.r)-g.r*(e2+z2+b2), -2*g.r*ez, g.nx,
            g.n*(1-3*mr/(g.r*g.r*g.r)+5*b2-z2+e2)};
}
