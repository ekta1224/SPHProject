#include <string>
#include <algorithm>

#include "pprint.cc"
#include "smallvec.cc"

/*
  Time-stamp: <mykern.cc on Tuesday, 10 March, 2015 at 00:11:10 MST (pinto)>

  provides constants and basic kernel functions for B-spline kernels of order 3, 4, and 5
   and Wendland C2, C4, and C6 kernels

  NB: Radius of support is unity for all kernels.

  Notes on gravity kernels:
     The potential kernel, phi, has a constant, a 1/q term, and then terms starting with q^2 up to q^(O+2)
     The force kernel, phip, has no constant, a 1/q^2 term, and then terms starting with q up to q^(O+1)
     The dphidh kernel, dphidh, has a constant, and then terms starting with q^2 up to q^(O+2)
     pc -- potential coefficients of q starting at pc[][2], pcconst[] constant, and pcinv[] the 1/q coefficient
     fc -- force coefficients of q starting at pc[][1], fcinv[] the 1/q^2 coefficient
     dc -- dphidh coefficients of q starting at pc[][2], dcconst[] the constant

*/

// kernel base class
// (actual kernels are defined below in derived classes)
template <int D, int NS, int O>
class kernels {
public:

    typedef SmallVec<double, D> dvec;

    kernels();
    void setdata();
    void kernelInfo();

    double geteta(int Ndesired);

    double kernel(const double &r, const double &h);
    double checkkernel(const double &r, const double &h);
    double gradKernel(const double &r, const double &h);
    void checkGradKernel(double r, double h);
    double dKerneldh(const double &r, const double &h);
    void checkdhkernel(double r, double h);

    double phikernel(const double &r, const double &h);
    double forcekernel(const double &r, const double &h);
    double dphikerneldh(const double &r, const double &h);

    std::string kernelname;

    static constexpr int ns = NS;
    static constexpr int nt = O+1;

    double cr[ns][nt], co[ns][nt], cd[ns][nt];
    double hs[ns];

    double pc[ns+1][nt+2], fc[ns+1][nt+2];
    double dc[ns+1][nt+2];
    double pcinv[ns+1], fcinv[ns+1];
    double pcconst[ns], dcconst[ns];

    double kernconst;
    double invD;

    //    static constexpr double KD[4] = {0, 2.0, M_PI, 4*M_PI/3.0};
    double KD[4];

};

template <int D, int NS, int O>
kernels<D,NS,O>::kernels() {
    invD = 1.0/(double)D;
    KD[0] = 0;
    KD[1] = 2.0;
    KD[2] = M_PI;
    KD[3] = 4*M_PI/3.0;
}

template <int D, int NS, int O>
double kernels<D,NS,O>::geteta(int Ndesired) {
    assert( D>0 && D<4 );
    double eta;

    eta = pow( (double)Ndesired/KD[D], 1.0/((double)D) );

    pprint("Ndesired: %d  KD: %e  eta = %e\n", Ndesired, KD[D], eta);

    return eta;
}

template <int D, int NS, int O>
void kernels<D,NS,O>::setdata() {

    // multiply by kernconst here so we don't have to in every kernel evaluation
    for(int s=0; s<ns; s++)
        for(int t=0; t<nt; t++) cr[s][t] *= kernconst;

    // manufacture the d(kernel)/dh and grad kernel coefficients
    for(int s=0; s<ns; s++) {
        for(int t=0; t<nt; t++) co[s][t] = -cr[s][t]*(t+D);
        for(int t=1; t<nt; t++) cd[s][t-1] = cr[s][t]*t;
    }

    // gravity only meaningful w/ D=3:
    if( D == 3 ) {

        // gravitational force kernel phi':
        for(int k=0; k<ns; k++) {
            fc[k][0] = 0;
            fc[k][nt+1] = 0;                        // there is no fc[k][nt+1], but we sum over it in gravity_sum
            for(int i=0; i<nt; i++) {
                fc[k][i+1] = 4*M_PI*cr[k][i]/(i+3); // coefficient of q^(i+1)
            }

            fcinv[k] = 0;                           // coefficient of 1/q^2 stored

            for(int j=0; j<k; j++) {
                for(int i=0; i<nt; i++) {
                    fcinv[k] += 4*M_PI*pow(hs[j], i+3)/(i+3) * (cr[j][i]-cr[j+1][i]);
                }
            }
        }

        // for k == ns, we are outside the kernel; these coefficients had better sum to 1!
        fc[ns][0] = 0;
        for(int j=0; j<ns-1; j++) {
            for(int i=0; i<nt; i++) {
                fcinv[ns] += 4*M_PI*pow(hs[j], i+3)/(i+3) * (cr[j][i]-cr[j+1][i]);
            }
        }
        for(int i=0; i<nt; i++) {
            fcinv[ns] += 4*M_PI*pow(hs[ns-1], i+3)/(i+3) * cr[ns-1][i];
        }
        pprint("fabs(1-fcinv[ns]) = %e\n", fabs(1-fcinv[ns]));
        assert( fabs(1-fcinv[ns]) < 1e-13 );

        // gravitational potential kernel phi:
        for(int k=0; k<ns; k++) {
            pc[k][0] = pc[k][1] = 0;
            for(int i=0; i<nt; i++) {
                pc[k][i+2] = 4*M_PI*cr[k][i]/((i+2)*(i+3)); // coefficient of q^(i+2)
            }

            pcconst[k] = 0;                                 // constant in polynomial
            for(int j=k; j<ns-1; j++) {
                for(int i=0; i<nt; i++) {
                    pcconst[k] += - 4*M_PI*pow(hs[j],i+2)/(i+2) * (cr[j][i] - cr[j+1][i]);
                }
            }
            for(int i=0; i<nt; i++) {
                pcconst[k] += - 4*M_PI*pow(hs[ns-1], i+2)/(i+2) * cr[ns-1][i];
            }

            pcinv[k] = 0;                                  // coefficient of 1/q
            for(int j=0; j<k; j++) {
                for(int i=0; i<nt; i++) {
                    pcinv[k] += - 4*M_PI*pow(hs[j],i+3)/(i+3) * ( cr[j][i] - cr[j+1][i] );
                }
            }
        }

        // for k == ns, we are outside the kernel; these coefficients had better sum to -1!
        pcinv[ns] = 0;
        for(int j=0; j<ns-1; j++) {
            for(int i=0; i<nt; i++) {
                pcinv[ns] += - 4*M_PI*pow(hs[j], i+3)/(i+3) * ( cr[j][i] - cr[j+1][i] );
            }
        }
        for(int i=0; i<nt; i++) {
            pcinv[ns] += - 4*M_PI*pow(hs[ns-1], i+3)/(i+3) * cr[ns-1][i];
        }
        assert( fabs(1+pcinv[ns]) < 1e-13 );


        // gravitational d(phi)/dh kernel = -1/h^2 (h phi - q h^2 phi')
        for(int k=0; k<ns; k++) {
            dc[k][0] = dc[k][1] = 0;
            for(int i=1; i<=nt+1; i++) {
                dc[k][i] = -( pc[k][i] + fc[k][i-1] ); // coefficient of q^i
            }

            dcconst[k] = -pcconst[k];                  // constant is -(constant in phi)
            assert( fabs(pcinv[k]+fcinv[k]) < 1e-13 ); // there is no 1/q^n term...
        }

    }
}

// report on kernel setup
template <int D, int NS, int O>
void kernels<D,NS,O>::kernelInfo() {
    pprint("--------------------------------------------------------------------------------\n");
    pprint("--------------------------------------------------------------------------------\n");
    pprint("\nKernel: %s\n\n", kernelname);

    pprint("      dimension: %d\n", D);
    pprint("   kernel order: %d\n", O);
    pprint("kernel segments: %d\n", ns);
    pprint("  kernel radius: %3.1f\n", 1.0);
    pprint("kernel normalization: %e\n", kernconst);
    pprint("\n");

    pprint("kernel polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("cr[%d] = {",i);
        for(int j=0; j<nt-1; j++) {
            pprint("% 14.8e, ", cr[i][j]+0);
        }
        pprint("% 14.8e}\n", cr[i][nt-1]+0);
    }
    pprint("\n");

    pprint("d(kernel)/dh polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("co[%d] = {",i);
        for(int j=0; j<nt-1; j++) {
            pprint("% 14.8e, ", co[i][j]+0);
        }
        pprint("% 14.8e}\n", co[i][nt-1]+0);
    }
    pprint("\n");

    pprint("grad kernel polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("cd[%d] = {",i);
        for(int j=0; j<nt-2; j++) {
            pprint("% 14.8e, ", cd[i][j]+0);
        }
        pprint("% 14.8e}\n", cd[i][nt-2]+0);
    }
    pprint("\n");

    pprint("gravitational potential kernel phi polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("pc[%d] = {",i);
        pprint("% 14.8e, ", pcconst[i]);
        for(int j=1; j<=nt; j++) {
            pprint("% 14.8e, ", pc[i][j]+0);
        }
        pprint("% 14.8e}", pc[i][nt+1]+0);
        pprint("   % 14.8e q^%d\n",pcinv[i], -1);
    }
    for(int j=0; j<=(10+17*(nt+2)); j++) pprint(" ");
    pprint("% 14.8e q^%d\n",pcinv[ns], -1);
    pprint("\n");

    pprint("gravitational force kernel phi' polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("fc[%d] = {",i);
        for(int j=0; j<nt; j++) {
            pprint("% 14.8e, ", fc[i][j]+0);
        }
        pprint("% 14.8e}", fc[i][nt]);
        pprint("   % 14.8e q^%d\n",fcinv[i]+0, -2);
    }
    for(int j=0; j<=(10+17*(nt+1)); j++) pprint(" ");
    pprint("% 14.8e q^%d\n",fcinv[ns], -2);
    pprint("\n");

    pprint("gravitational d(phi')/dh kernel polynomial coefficients:\n");
    for(int i=0; i<ns; i++) {
        pprint("dc[%d] = {",i);
        pprint("% 14.8e, ", dcconst[i]);
        for(int j=1; j<=nt; j++) {
            pprint("% 14.8e, ", dc[i][j]+0); // takes care of -0 (!?!)
        }
        pprint("% 14.8e}\n", dc[i][nt+1]+0);
    }
    pprint("--------------------------------------------------------------------------------\n");
    pprint("--------------------------------------------------------------------------------\n");
}

// the smoothing kernel
template <int D, int NS, int O>
inline double kernels<D,NS,O>::kernel(const double &r, const double &h) {
    assert(r>=0); assert(h>0);

    double kern = 0;
    double q = r/h;

    for(int s=0; s<ns; s++) {
        if(q<hs[s]) {
            double fac = 1;
            for(int t=0; t<nt; t++) {
                kern += cr[s][t]*fac;
                fac *= q;
            }
            break;
        }
    }

    if(kern<0) kern = 0;

    return kern * pow(h,-D);
}

// gradkernel returns the magnitude of grad W such that the vector quantity
//     Nabla W(|r|,h) = rhat F(|r|,h)
template <int D, int NS, int O>
inline double kernels<D,NS,O>::gradKernel(const double &r, const double &h) {
    assert(r>=0); assert(h>0);

    double gradkern = 0;
    double q = r/h;

    for(int s=0; s<ns; s++) {
        if(q<hs[s]) {
            double fac = 1;
            for(int t=0; t<nt-1; t++) {
                gradkern += cd[s][t]*fac;
                fac *= q;
            }
            break;
        }
    }

    return gradkern * pow(h,-D-1);
}

// derivative of smoothing kernel w/r to smoothing length
template <int D, int NS, int O>
inline double kernels<D,NS,O>::dKerneldh(const double &r, const double &h) {
    assert(r>=0); assert(h>0);

    double dkerndh = 0;
    double q = r/h;

    for(int s=0; s<ns; s++) {
        if(q < hs[s]) {
            double fac = 1;
            for(int t=0; t<nt; t++) {
                dkerndh += co[s][t]*fac;
                fac *= q;
            }
            break;
        }
    }
    return dkerndh * pow(h,-D-1);
}

// the gravitational potential kernel phi (only defined for D=3)
template <int D, int NS, int O>
inline double kernels<D,NS,O>::phikernel(const double &r, const double &h) {
    assert( D==3 ); assert(r>=0); assert(h>0);

    double phikern = 0;
    double hinv = 1.0/h;
    double q = r * hinv;

    for(int s=0; s<ns; s++) {
        if(q<hs[s]) {
            double fac = q*q;
            for(int t=2; t<=nt+1; t++) {
                phikern += pc[s][t]*fac;
                fac *= q;
            }
            phikern += pcinv[s]/q;
            phikern += pcconst[s];
            break;
        }
    }
    if(q>hs[ns-1]) phikern = -1./q;

    return phikern * hinv;
}

// the gravitational force kernel phi' (only defined for D=3)
template <int D, int NS, int O>
inline double kernels<D,NS,O>::forcekernel(const double &r, const double &h) {
    assert( D==3 ); assert(r>=0); assert(h>0);

    double forcekern = 0;
    double hinv = 1.0/h;
    double q = r * hinv;

    for(int s=0; s<ns; s++) {
        if(q<hs[s]) {
            double fac = q;
            for(int t=1; t<nt+1; t++) {
                forcekern += fc[s][t]*fac;
                fac *= q;
            }
            forcekern += fcinv[s]/(q*q);
            break;
        }
    }
    if(q>hs[ns-1]) forcekern = 1./(q*q);

    return forcekern * hinv * hinv;
}

// the gravitational d(phi)/dh kernel (only defined for D=3)
template <int D, int NS, int O>
inline double kernels<D,NS,O>::dphikerneldh(const double &r, const double &h) {
    assert( D==3 ); assert(r>=0); assert(h>0);

    double dphikern = 0;
    double hinv = 1.0/h;
    double q = r * hinv;

    for(int s=0; s<ns; s++) {
        if(q<hs[s]) {
            double fac = q*q;
            for(int t=2; t<=nt+1; t++) {
                dphikern += dc[s][t]*fac;
                fac *= q;
            }
            dphikern += dcconst[s];
            break;
        }
    }

    return dphikern * hinv * hinv;
}

// --------------------------------------------------------------------------------

// Individual kernels:

// --------------------------------------------------------------------------------

// cubic B-spline kernel
template <int D>
class CubicSpline : public kernels<D,2,3> {
public:
    using kernels<D,2,3>::cr;
    using kernels<D,2,3>::kernconst;
    using kernels<D,2,3>::hs;
    using kernels<D,2,3>::setdata;
    using kernels<D,2,3>::kernelname;

    CubicSpline() {

        kernelname = "Cubic Spline";

        if(D==1)
            kernconst = 8./3.;
        else if(D==2)
            kernconst = 80./(7.*M_PI);
        else if(D==3)
            kernconst = 16./M_PI;
        else {
            fpprint(std::cerr,"Cubic B-spline constructor: bad dimension: D=%d\n", D);
            exit(1);
        }

        hs[0] = 1.0/2.0;
        hs[1] = 1.0;

        cr[0][0] = 0.5;  cr[0][1] =  0.;  cr[0][2] = -3.;  cr[0][3] =  3.;
        cr[1][0] = 1.0;  cr[1][1] = -3.;  cr[1][2] =  3.;  cr[1][3] = -1.;

        setdata();
    }

    double analytic_kernel(double r, double h) {
        assert( r >=0 ); assert( h > 0 );
        double q = r/h;
        return kernconst*(-4*pow(std::max(0.0,0.5-q),3) + pow(std::max(0.0,1.0-q),3))*pow(h,-D);
    }
    double analytic_dkerneldh(double r, double h) {
        assert( r >=0 ); assert( h > 0 );
        double q = r/h;
        if(q<0.5) {
            return kernconst * (-1.5*(h*h*h - 10*h*r*r + 12*r*r*r)*pow(h,-7));
        }
        else if(q<1.0) {
            return kernconst * (-3*(h-2*r)*(h-r)*(h-r)*pow(h,-7));
        }
        return 0;
    }
};

// cubic B-spline kernel with radius of support r/h = 2
template <int D>
class CubicSplineOn2 : public kernels<D,2,3> {
public:
    using kernels<D,2,3>::cr;
    using kernels<D,2,3>::kernconst;
    using kernels<D,2,3>::hs;
    using kernels<D,2,3>::setdata;
    using kernels<D,2,3>::kernelname;

    CubicSplineOn2() {

        kernelname = "Cubic Spline";

        if(D==1)
            kernconst = 2./3.;
        else if(D==2)
            kernconst = 10./(7.*M_PI);
        else if(D==3)
            kernconst = 1./M_PI;
        else {
            fpprint(std::cerr,"Cubic B-spline constructor: bad dimension: D=%d\n", D);
            exit(1);
        }

        hs[0] = 1.0;
        hs[1] = 2.0;

        cr[0][0] = 1.;  cr[0][1] =  0.;  cr[0][2] = -1.5;  cr[0][3] =  0.75;
        cr[1][0] = 2.;  cr[1][1] = -3.;  cr[1][2] =  1.5;  cr[1][3] = -0.25;

        setdata();
    }
};



// quartic B-spline kernel
template <int D>
class QuarticSpline : public kernels<D,3,4> {
public:
    using kernels<D,3,4>::cr;
    using kernels<D,3,4>::kernconst;
    using kernels<D,3,4>::hs;
    using kernels<D,3,4>::setdata;
    using kernels<D,3,4>::kernelname;

    QuarticSpline() {

        kernelname = "Quartic Spline";

        if(D==1)
            kernconst = 3125./768.;
        else if(D==2)
            kernconst = 46875./(2398.*M_PI);
        else if(D==3)
            kernconst = 15625./(512.*M_PI);
        else {
            fpprint(std::cerr,"Quartic B-spline constructor: bad dimension: D=%d\n", D);
            exit(1);
        }

        hs[0] = 1.0/5.0;
        hs[1] = 3.0/5.0;
        hs[2] = 1.0;

        cr[0][0] =  46./125.;  cr[0][1] =      0.;  cr[0][2] = -12./5.;  cr[0][3] =  0.;  cr[0][4] =  6.;
        cr[1][0] =  44./125.;  cr[1][1] =  8./25.;  cr[1][2] = -24./5.;  cr[1][3] =  8.;  cr[1][4] = -4.;
        cr[2][0] =        1.;  cr[2][1] =     -4.;  cr[2][2] =      6.;  cr[2][3] = -4.;  cr[2][4] =  1.;

        setdata();
    }

};

// quintic B-spline kernel
template <int D>
class QuinticSpline : public kernels<D,3,5> {
public:
    using kernels<D,3,5>::cr;
    using kernels<D,3,5>::kernconst;
    using kernels<D,3,5>::hs;
    using kernels<D,3,5>::setdata;
    using kernels<D,3,5>::kernelname;

    QuinticSpline() {

        kernelname = "Quintic Spline";

        if(D==1)
            kernconst = 243./40.;
        else if(D==2)
            kernconst = 15309./(478.*M_PI);
        else if(D==3)
            kernconst = 2187./(40.*M_PI);
        else {
            fpprint(std::cerr,"Quintic B-spline constructor: bad dimension: D=%d\n", D);
            exit(1);
        }

        hs[0] = 1.0/3.0;
        hs[1] = 2.0/3.0;
        hs[2] = 1.0;

        cr[0][0] = 22./81.;  cr[0][1] =      0.;  cr[0][2] = -20./9.;  cr[0][3] =     0.;  cr[0][4] =   10.;  cr[0][5] =  -10.;
        cr[1][0] = 17./81.;  cr[1][1] = 25./27.;  cr[1][2] = -70./9.;  cr[1][3] = 50./3.;  cr[1][4] =  -15.;  cr[1][5] =    5.;
        cr[2][0] =      1.;  cr[2][1] =     -5.;  cr[2][2] =     10.;  cr[2][3] =   -10.;  cr[2][4] =    5.;  cr[2][5] =   -1.;

        setdata();
    }
};


// C2 Wendland kernel for D>1
template <int D>
class C2 : public kernels<D,1,5> {
public:
    using kernels<D,1,5>::cr;
    using kernels<D,1,5>::kernconst;
    using kernels<D,1,5>::hs;
    using kernels<D,1,5>::setdata;
    using kernels<D,1,5>::kernelname;

    C2() {

        kernelname = "C2 Wendland D>1";

        if(D==2) {
            kernconst = 7./M_PI;
        }
        else if(D==3) {
            kernconst = 21.0/(2.*M_PI);
        }
        else {
            fpprint(std::cerr,"Wendland C2 constructor:  bad dimension: D=%d\n", D);
        }

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -10.; cr[0][3] = 20.; cr[0][4] = -15.; cr[0][5] = 4.;

        setdata();
    }
};

// C2 Wendland kernel for D=1
template <>
class C2<1> : public kernels<1,1,4> {
public:
    using kernels<1,1,4>::cr;
    using kernels<1,1,4>::kernconst;
    using kernels<1,1,4>::hs;
    using kernels<1,1,4>::setdata;
    using kernels<1,1,4>::kernelname;

    C2() {

        kernelname = "C2 Wendland D=1";

        kernconst = 5./4.;

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -6.; cr[0][3] = 8.; cr[0][4] = -3.;

        setdata();
    }
};


// C4 Wendland kernel for D>1
template <int D>
class C4 : public kernels<D,1,8> {
public:
    using kernels<D,1,8>::cr;
    using kernels<D,1,8>::kernconst;
    using kernels<D,1,8>::hs;
    using kernels<D,1,8>::setdata;
    using kernels<D,1,8>::kernelname;

    C4() {

        kernelname = "C4 Wendland D>1";

        if(D==2) {
            kernconst = 9./M_PI;
        }
        else if(D==3) {
            kernconst = 495./(32.*M_PI);
        }
        else {
            fpprint(std::cerr,"Wendland C4 constructor:  bad dimension: D=%d\n", D);
        }

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -28./3.; cr[0][3] = 0.; cr[0][4] = 70.; cr[0][5] = -448./3.;
        cr[0][6] = 140.; cr[0][7] = -64.; cr[0][8] = 35./3.;

        setdata();
    }
};


// C4 Wendland kernel for D=1
template <>
class C4<1> : public kernels<1,1,7> {
public:
    using kernels<1,1,7>::cr;
    using kernels<1,1,7>::kernconst;
    using kernels<1,1,7>::hs;
    using kernels<1,1,7>::setdata;
    using kernels<1,1,7>::kernelname;

    C4() {

        kernelname = "C4 Wendland D=1";

        kernconst = 3./2.;

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -7.; cr[0][3] = 0.; cr[0][4] = 35.; cr[0][5] = -56;
        cr[0][6] = 35; cr[0][7] = -8.;

        setdata();
    }
};


// C6 Wendland kernel for D>1
template <int D>
class C6 : public kernels<D,1,11> {
public:
    using kernels<D,1,11>::cr;
    using kernels<D,1,11>::kernconst;
    using kernels<D,1,11>::hs;
    using kernels<D,1,11>::setdata;
    using kernels<D,1,11>::kernelname;

    C6() {

        kernelname = "C6 Wendland D>1";

        if(D==2) {
            kernconst = 78./(7.*M_PI);
        }
        else if(D==3) {
            kernconst = 1365./(64.*M_PI);
        }
        else {
            fpprint(std::cerr,"Wendland C6 constructor:  bad dimension: D=%d\n", D);
        }

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -11.; cr[0][3] = 0.; cr[0][4] = 66.; cr[0][5] = 0.;
        cr[0][6] = -462.; cr[0][7] = 1056.; cr[0][8] = -1155.; cr[0][9] = 704.; cr[0][10] = -231.; cr[0][11] = 32.;

        setdata();
    }
};


// C6 Wendland kernel for D=1
template <>
class C6<1> : public kernels<1,1,10> {
public:
    using kernels<1,1,10>::cr;
    using kernels<1,1,10>::kernconst;
    using kernels<1,1,10>::hs;
    using kernels<1,1,10>::setdata;
    using kernels<1,1,10>::kernelname;

    C6() {

        kernelname = "C6 Wendland D=1";

        kernconst = 55./32.;

        hs[0] = 1.0;

        cr[0][0] = 1.; cr[0][1] = 0.; cr[0][2] = -9.; cr[0][3] = 0.; cr[0][4] = 42.; cr[0][5] = 0.;
        cr[0][6] = -210.; cr[0][7] = 348.; cr[0][8] = -315.; cr[0][9] = 128.; cr[0][10] = -21.;

        setdata();
    }
};

// --------------------------------------------------------------------------------
