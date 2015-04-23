#if 0
/*
Time-stamp: <sphdensity.cc on Thursday, 9 April, 2015 at 10:28:55 MST (pinto)>

  Example of sph density determination
     using O(N^2) loops and Newton-Raphson iteration

 Version which works for 1, 2, and 3 dimensions

 */

#include <cassert>
#include "pprint.cc"
#include "smallvec.cc"
#include "mykern.cc"

#include "particleset.cc"
#include "getparticles.cc"
#endif

// class to find SPH density and smoothing lengths of a set of particles
template <class MyKernel, int D>
class SPHDensity : public MyKernel {
public:
    SPHDensity(double eta) : eta(eta) {}

    void rhosums(dvec sink_r, dvec *source_r, double *source_m, int nsource, double h,
                 double &rho, double &omega, int &nn);
    void funcd(dvec sink_r, double sink_m, dvec *source_r, double *source_m, int nsource, double h,
               double &f, double &df, double &rho, double &omega, int &nn);
    double newton(dvec sink_r, double sink_m, dvec *source_r, double *source_m, int nsource, double h,
                  double &rho, double &omega, int &nn);
    void getrho(dvec sink_r, double sink_m, dvec *source_r, double *source_m, int nsource,
                double &h, double &rho, double &omega, int &nn);
    double rtsafe(const dvec sink_r, const double sink_m, dvec *source_r, double *source_m, int nsource, const double hlower, const double hupper,
                  double &rho, double &omega, int &nn);

    double eta;
};

template <class MyKernel, int D>
void SPHDensity<MyKernel, D>::getrho(dvec sink_r, double sink_m, dvec *source_r, double *source_m, int nsource,
                               double &h, double &rho, double &omega, int &nn) {

    assert(h>0);
    double newh = newton(sink_r, sink_m, source_r, source_m, nsource, h, rho, omega, nn);
    if(newh<0) {
        //        fpprint(std::cerr,"trying rtsafe\n");
        double l = 0.5 * h;
        double u = 2   * h;
        int it = 0;
        while(it<10) {
            newh = rtsafe(sink_r, sink_m, source_r, source_m, nsource, l, u, rho, omega, nn);
            if(newh > 0)  break;
            l *= 0.5;
            u *= 2;
            it++;
        }
        if(it==10) {
            fpprint(std::cerr,"finding h failed:\n");
            fpprint(std::cerr,"% e % e\n", l, u);
            exit(1);
        }
        //            pprint("finished rtsafe\n");
    }

    assert(newh > 0);
    h = newh;
}

// get SPH density, omega, and number of neighbours for particle i
template <class MyKernel, int D>
void SPHDensity<MyKernel, D>::rhosums(const dvec sink_r, dvec *source_r, double *source_m, int nsource, double h,
                                double &rho, double &omega, int &nn) {
    assert( h > 0);
    rho = 0;
    omega = 0;
    nn = 0;
    for(int j=0; j<nsource; j++) {
        double rij = (sink_r - source_r[j]).norm();
        rho += source_m[j] * this->kernel(rij, h);
        omega += source_m[j] * this->dKerneldh(rij, h);
        if(rij/h < 1.0) {
            nn++;
        }
    }
    omega = 1.0 + h/(D*rho)*omega;
}

template <class MyKernel, int D>
void SPHDensity<MyKernel, D>::funcd(const dvec sink_r, const double sink_m, dvec *source_r, double *source_m, int nsource, double h,
                              double &f, double &df, double &rho, double &omega, int &nn) {
    assert(h>0);
    rhosums(sink_r, source_r, source_m, nsource, h, rho, omega, nn);
    f = sink_m * pow(eta/h, D) - rho;
    df = -(D/h)*rho*omega;
}

template <class MyKernel, int D>
double SPHDensity<MyKernel, D>::newton(const dvec sink_r, const double sink_m, dvec *source_r, double *source_m, int nsource, double h,
                                 double &rho, double &omega, int &nn) {
    assert(h>0);
    int it = 0;
    double dh = 1e30;
    while(fabs(dh)/h > 1e-8) {
        double f, df;
        assert(h>0);
        funcd(sink_r, sink_m, source_r, source_m, nsource, h, f, df, rho, omega, nn);
        dh = -f/df;
        h += dh;
        it++;
        if(it>100 || !std::isfinite(h) || h <= 0 || h > 100) return -1;
    }
    return h;
}

template <class MyKernel, int D>
double SPHDensity<MyKernel, D>::rtsafe(const dvec sink_r, const double sink_m,
                                       dvec *source_r, double *source_m, int nsource,
                                       const double hlower, const double hupper,
                                       double &rho, double &omega, int &nn) {
    const double htol = 1.0e-08;
    const int maxhiterations = 100;

    double fl, fu, hl, hu, f, df;
    funcd(sink_r, sink_m, source_r, source_m, nsource, hlower, fl, df, rho, omega, nn);
    funcd(sink_r, sink_m, source_r, source_m, nsource, hupper, fu, df, rho, omega, nn);
    if ((fl > 0.0 && fu > 0.0) || (fl < 0.0 && fu < 0.0)) {
        fpprint(std::cerr,"rtsafe: root must be bracketed by input\n");
        return -1;
    }

    if (fl == 0.0) return hlower;
    if (fu == 0.0) return hupper;
    if (fl < 0.0) {
        hl=hlower;  hu=hupper;
    } else {
        hu=hlower;  hl=hupper;
    }

    double hm = 0.5*(hlower+hupper);
    double dhold = fabs(hupper-hlower);
    double dh = dhold;
    funcd(sink_r, sink_m, source_r, source_m, nsource, hm, f, df, rho, omega, nn);

    for (int j=0; j< maxhiterations; j++) {
        if ( ( ((hm-hu)*df-f)*((hm-hl)*df-f) > 0.0 ) || ( fabs(2.0*f) > fabs(dhold*df) ) ) {
            dhold = dh;
            dh = 0.5*(hu-hl);
            hm = hl + dh;
            if ( hl == hm ) return hm;
        }
        else {
            dhold = dh;
            dh = f/df;
            double temp = hm;
            hm -= dh;
            if ( temp == hm ) return hm;
        }
        if ( fabs(dh) < htol ) return hm;

        funcd(sink_r, sink_m, source_r, source_m, nsource, hm, f, df, rho, omega, nn);
        if ( f < 0.0 )
            hl=hm;
        else
            hu=hm;
    }

    fpprint(std::cerr,"rtsafe: maximum number of iterations exceeded\n");
    exit(1);
}
