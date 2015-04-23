/*
Time-stamp: <sphdensityND.cc on Monday, 9 March, 2015 at 23:54:59 MST (pinto)>

  Example of sph density determination
     using O(N^2) loops and Newton-Raphson iteration

 Version which works for 1, 2, and 3 dimensions

 */

#include <cmath>
#include <cassert>
//#include "pprint.cc"
//#include "smallvec.cc"

// class implementing Gaussian kernel
template <int D>
class Gaussian {
public:
    // Gaussian kernel normalized to unit integral over space
    double W(double rij, double h) {
        double q = rij/h;
        return exp(-q*q)/(pow(M_PI,0.5*D) * pow(h,D));
    }

    // and its h-derivative
    double dWdh(double rij, double h) {
        return W(rij, h)*(D*h*h-2*rij*rij)/(h*h*h);
    }

};

// class implementing M3 cubic spline kernel with support q = r/h on [0,1)
template <int D>
class CubicSpline {
public:
    // cubic spline kernel
    double W(double rij, double h) {

        double q = rij/h;
        double w;
        if(q<=0.5)
            w =  0.5 - 3*q*q + 3*q*q*q;
        else if(q<=1.0)
            w =  1.0 - 3*q + 3*q*q - q*q*q;
        else
            w = 0;

        return knorm[D]*pow(h, -D) * w;
    }

    // and its h-derivative
    double dWdh(double rij, double h) {
        double q = rij/h;
        double dW;
        if(q<=0.5)
            dW =  -1.5 + 15*q*q - 18*q*q*q;
        else if(q<=1.0)
            dW =  -3 + 12*q - 15*q*q + 6*q*q*q;
        else
            dW = 0;

        return knorm[D]*pow(h, -(D+1)) * dW;
    }

    // normalization constants for all three dimensions (first element is zero since we have D=1,2,or 3)
    const double knorm[4] = { 0, 8.0/3.0, 16.0*5.0/(7.0*M_PI), 16.0/M_PI };
};

// a struct to keep particle data:
template <int D>
struct Particle {
    SmallVec<double, D> r;
    double m;
    double rho, omega;
    double h;
    int id;
    int nn;
};

// fill a unit D-cube (-0.5,0.5)^D with n1d^D particles on a grid
template <int D>
Particle<D> *unifill(int n1d, int &N) {

    double eps = 1e-6;

    N = pow(n1d, D);
    double mpp = 1.0/(double)N;
    double dx = 1.0/(double)(n1d);

    Particle<D> *p = new Particle<D>[N];

    SmallVec<int, D> idx(0);

    int m = 0;
    while(1) {
        for(int i=0; i<D; i++) {
            double q = -0.5 + 0.5*dx + idx[i]*dx;
            q += eps*(1.0-0.5*drand48());
            p[m].r[i] = q;
        }
        p[m].m = mpp;
        p[m].h = dx;
        m++;

        int i;
        for(i=D-1; i>=0; i--) {
            idx[i]++;
            if(idx[i] < n1d) break;
            idx[i] = 0;
        }
        if( i < 0 ) break;
    }
    assert( m == N );

    return p;
}

// class to find SPH density and smoothing lengths of a set of particles
template <class MyKernel, int D>
class Doit : public MyKernel {
public:
    Doit(double eta, Particle<D> *p, int N) : eta(eta), p(p), N(N) {}
    double rhosmooth(double h);
    void rhosums(int i, double h, double &rho, double &omega, int &nn);
    void funcd(int i, double h, double &f, double &df, double &rho, double &omega, int &nn);
    void newton(int i, double &h, double &rho, double &omega, int &nn);
    void getrho(int i);

    int N;
    Particle<D> *p;
    double eta;
};

template <class MyKernel, int D>
double Doit<MyKernel, D>::rhosmooth(double h) {
    return pow(eta/h, D);
}

// get SPH density, omega, and number of neighbours for particle i
template <class MyKernel, int D>
void Doit<MyKernel, D>::rhosums(int i, double h, double &rho, double &omega, int &nn) {
    rho = 0;
    omega = 0;
    nn = 0;
    for(int j=0; j<N; j++) {
        double rij = (p[i].r - p[j].r).norm();
        rho += p[j].m * this->W(rij, h);
        omega += p[j].m * this->dWdh(rij, h);
        if(rij/h < 1.0) {
            nn++;
        }
    }
    omega = 1.0 + h/(D*p[i].m*rhosmooth(h))*omega;
}

template <class MyKernel, int D>
void Doit<MyKernel, D>::funcd(int i, double h, double &f, double &df, double &rho, double &omega, int &nn) {
    rhosums(i, h, rho, omega, nn);
    f = p[i].m * rhosmooth(h) - rho;
    df = -D*omega*p[i].m*rhosmooth(h)/h;
}

template <class MyKernel, int D>
void Doit<MyKernel, D>::newton(int i, double &h, double &rho, double &omega, int &nn) {
    int it = 0;
    double dh = 1e30;

    while(fabs(dh)/h > 1e-8) {
        double f, df;
        funcd(i, h, f, df, rho, omega, nn);
        dh = -f/df;
        h += dh;
        it++;
    }
}

template <class MyKernel, int D>
void Doit<MyKernel, D>::getrho(int i) {
    double h, rho, omega;
    int nn;
    h = p[i].h;
    newton(i, h, rho, omega, nn);
    p[i].h = h;
    p[i].rho = rho;
    p[i].omega = omega;
    p[i].nn = nn;
}

int main() {

    const int D = 3;
    int n1d = pow(1<<14, 1.0/(double)D); // make about 16K particles, no matter what dimension D
    int N;
    Particle<D> *p = unifill<D>(n1d, N);

    double eta = 2.5; // should give about 80 neighbours in 3D

    // get instance of Doit using desired kernel function
    //    Doit<Gaussian<D>, D> doit(eta, p, N);
    Doit<CubicSpline<D>, D> doit(eta, p, N);

#pragma omp parallel for
    for(int i=0; i<N; i++) {
        doit.getrho(i);
#pragma omp critical
        if(p[i].r.norm() < 0.3) pprint("%e  %e  %d\n", p[i].r.norm(), p[i].rho, p[i].nn);
    }
}
