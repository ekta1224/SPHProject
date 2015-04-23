/*
Time-stamp: <foo.cc on Monday, 9 March, 2015 at 23:57:36 MST (pinto)>

  Example of sph density determination
     using O(N^2) loops and Newton-Raphson iteration

 */

#include <cmath>
#include <cassert>
//#include "pprint.cc"
//#include "smallvec.cc"
#include <cstdlib>

// class implementing Gaussian kernel
class Gaussian {
public:
    // Gaussian kernel normalized to unit integral over space
    double W(double rij, double h) {
        double q = rij/h;
        return exp(-q*q)/(pow(M_PI,1.5) * pow(h,3));
    }

    // and its h-derivative
    double dWdh(double rij, double h) {
        return W(rij, h)*(3*h*h-2*rij*rij)/(h*h*h);
    }

};

// class implementing M3 cubic spline kernel with support q = r/h on [0,1)
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

        return (16.0/M_PI)*pow(h, -3) * w;
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

        return (16.0/M_PI)*pow(h, -4) * dW;
    }

};

// a struct to keep particle data:
struct Particle {
    SmallVec<double, 3> r;
    double m;
    double rho, omega;
    double h;
    int id;
    int nn;
};

// fill a unit cube (-0.5,0.5)^3 with n1d^3 particles on a (slightly perturbed) grid
Particle *unifill(int n1d, int &N) {

    double eps = 1e-6;

    N = pow(n1d, 3);
    double mpp = 1.0/(double)N;
    double dx = 1.0/(double)(n1d);

    Particle *p = new Particle[N];

    int m = 0;
    for(int i=0; i<n1d; i++) {
        double x = -0.5 + 0.5*dx + i*dx;
        for(int j=0; j<n1d; j++) {
            double y = -0.5 + 0.5*dx + j*dx;
            for(int k=0; k<n1d; k++) {
                double z = -0.5 + 0.5*dx + k*dx;
                x += eps*(1.0-0.5*drand48());
                y += eps*(1.0-0.5*drand48());
                z += eps*(1.0-0.5*drand48());
                p[m].r = {x, y, z};
                p[m].m = mpp;
                p[m].h = dx;
                m++;
            }
        }
    }

    assert( m == N );

    return p;
}

// class to find SPH density and smoothing lengths of a set of particles
template <class MyKernel>
class Doit : public MyKernel {
public:
    Doit(double eta, Particle *p, int N) : eta(eta), p(p), N(N) {}
    double rhosmooth(double h);
    void rhosums(int i, double h, double &rho, double &omega, int &nn);
    void funcd(int i, double h, double &f, double &df, double &rho, double &omega, int &nn);
    void newton(int i, double &h, double &rho, double &omega, int &nn);
    void getrho(int i);

    int N;
    Particle *p;
    double eta;
};

template <class MyKernel>
double Doit<MyKernel>::rhosmooth(double h) {
    return pow(eta/h, 3);
}

// get SPH density, omega, and number of neighbours for particle i
template <class MyKernel>
void Doit<MyKernel>::rhosums(int i, double h, double &rho, double &omega, int &nn) {
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
    omega = 1.0 + h/(3*p[i].m*rhosmooth(h))*omega;
}

template <class MyKernel>
void Doit<MyKernel>::funcd(int i, double h, double &f, double &df, double &rho, double &omega, int &nn) {
    rhosums(i, h, rho, omega, nn);
    f = p[i].m * rhosmooth(h) - rho;
    df = -3*omega*p[i].m*rhosmooth(h)/h;
}

template <class MyKernel>
void Doit<MyKernel>::newton(int i, double &h, double &rho, double &omega, int &nn) {
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

template <class MyKernel>
void Doit<MyKernel>::getrho(int i) {
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

    int n1d = pow(1<<14, 1.0/3.0); // make about 16K particles
    int N;
    Particle *p = unifill(n1d, N);

    double eta = 2.5; // should give about 80 neighbours in 3D

    // get instance of Doit using desired kernel function
    //    Doit<Gaussian> doit(eta, p, N);
    Doit<CubicSpline> doit(eta, p, N);

#pragma omp parallel for
    for(int i=0; i<N; i++) {
        doit.getrho(i);
#pragma omp critical
        if(p[i].r.norm() < 0.3) pprint("%e  %e  %d\n", p[i].r.norm(), p[i].rho, p[i].nn);
    }
}
