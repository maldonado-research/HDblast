// Included inside the archived Solver to reuse its stage geometry and nodal fluxes.
void probe_cell(const vector<double>& y, vector<double>& dy, int cell,
                const Geo (&geom)[P], double emitted, double gauge) {
    for (int k=0; k<P; ++k) {
        const int node = cell*P+k;
        const Geo &g = geom[k];
        const double invr = 1./g.r;
        const double in_speed = 2.*g.n/dx;
        const double se_z = g.nx + 2.*g.n*g.rx*invr;
        const double sz_e = g.nx + g.n*g.rx*invr;
        const double expansion = g.n*g.u*invr;
        const double potential = g.n*invr*invr*invr;
        for (int m=0; m<probe_count; ++m) {
            const int ii = 3+m, oo = ii+probe_count, aa = oo+probe_count;
            const double inc=y[index(node,ii)], out=y[index(node,oo)];
            double di=0., dout=0.;
            for (int j=0; j<P; ++j) {
                di += D[k][j]*y[index(cell*P+j,ii)];
                dout += D[k][j]*y[index(cell*P+j,oo)];
            }
            if (k==P-1) {
                const double right = cell+1==cells ? (m==unforced ? 0. : emitted) : y[index(node+1,ii)];
                di += (right-inc)/weight[k];
            }
            if (k==0) {
                const double left = cell==0 ? y[PROBE_MOM+m]/gauge-inc : y[index(node-1,oo)];
                dout -= (left-out)/weight[k];
            }
            const double e=inc+out, z=inc-out;
            const double se=se_z*z-2.*expansion*e-probe_lambda[m]*probe_lambda[m]*potential*y[index(node,aa)];
            const double sz=sz_e*e-expansion*z;
            dy[index(node,ii)]=in_speed*di+.5*(se+sz);
            dy[index(node,oo)]=-in_speed*dout+.5*(se-sz);
            dy[index(node,aa)]=g.n*g.r*e;
        }
    }
}

void probe_wall(const vector<double>& y, vector<double>& dy, const Values& w, double a, double h) {
    for (int m=0; m<probe_count; ++m) {
        const double inc=y[index(0,3+m)], mom=y[PROBE_MOM+m], e=mom/w.k;
        const double lambda=probe_lambda[m], c=y[PROBE_C+m], b=lambda*c/a;
        dy[PROBE_C+m]=e-h*c;
        dy[PROBE_MOM+m]=-2.*h*mom-lambda*lambda*w.k*c/(a*a)-2.*(e-2.*inc);
        const double a4=a*a*a*a;
        dy[PROBE_NET+m]=2.*a4*(inc*inc-(e-inc)*(e-inc));
        dy[PROBE_DRIVER+m]=.5*a4*w.dk*y[PI]*(b*b-e*e);
    }
}
