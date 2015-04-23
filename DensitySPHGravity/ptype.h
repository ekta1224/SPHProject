#ifndef __PTYPE_INCLUDE__
#define __PTYPE_INCLUDE__

template <int D>
struct ptype {
    SmallVec<double, D> r, v;
    SmallVec<double, D> acc, acc0;
    double pot, pot0;
    double m;
    double rho, omega, h;
    int nn;
    int id;
};

#endif // __PTYPE_INCLUDE__
