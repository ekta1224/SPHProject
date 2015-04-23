/*
Time-stamp: <grav.cc on Thursday, 12 March, 2015 at 10:50:18 MST (pinto)>

g++ -DSORTED -std=c++11 -O3 -Ilib -fopenmp main.cc

 */


#include <cassert>
#include <omp.h>
#include <cstring>
#include <algorithm>
#include "pprint.cc"
#include "smallvec.cc"
#include "STimer.cc"
#include "minheap.cc"

template <int D>
struct ptype {
    SmallVec<double, D> r;
    SmallVec<double, D> acc;
    double m;
    int id;
};

#include "bh.cc"

template <int D>
void mkpoints(SmallVec<double, D> lowerleft, SmallVec<double, D> upperright, int N,
              ptype<D> *p) {
    pprint("making %d uniform particles on %e to %e\n", N, lowerleft, upperright);

    double mpp = 1.0/(double)N;

    SmallVec<double, D> delta = upperright - lowerleft;
    for(int i=0; i<N; i++) {
        for(int d=0; d<D; d++) {
            double coord;
            do {
                coord = lowerleft[d] + drand48()*delta[d];
            } while( !(coord>=lowerleft[d] && coord<upperright[d]) );
            p[i].r[d] = coord;
        }
        p[i].m = mpp;
        p[i].id = i;
    }
}


template <int D>
SmallVec<double, D> analyticAcc(int sink, double theta, double eps, ptype<D> *p, int N) {
    double eps2 = eps*eps;
    SmallVec<double,D> acc(0);

    for(int j=0; j<N; j++) {
        if( sink!=j ) {
            SmallVec<double, D> dr = p[sink].r - p[j].r;
            double nr = sqrt(dr.norm2() + eps2);
            double ir3 = p[j].m/(nr*nr*nr);
            acc -= ir3*dr;
        }
    }
    return acc;
}

int main(int argc, char **argv) {

    if(argc!=3) { pprint("usage %s <number of points> <maxleafsize>\n", argv[0]); exit(1); }
    int N = atoi(argv[1]);
    int maxleafsize = atoi(argv[2]);

    const int D = 3;
    typedef SmallVec<double, D> dvec;

    dvec center(0);
    double halfwidth = 1;
    ptype<D> *p = new ptype<D>[N];
    mkpoints(center-dvec(halfwidth), center+dvec(halfwidth), N, p);

    Octree<D> tree;
    tree.initTree(center, halfwidth);
    STimer foo; foo.Clear(); foo.Start();
    tree.buildTree(p, N, maxleafsize);
    foo.Stop(); pprint("time to build tree: %e  rate: %e\n", foo.Elapsed(), (double)N/foo.Elapsed());
    tree.checkTree();
    pprint("nnodes: %d\n", tree.n_nodes);
    pprint("nleaves: %d\n", tree.n_leafSet);

    tree.propagateCOM(0);
    tree.propagateBmax(0);

    double eps = 0.01;
    double theta = 0.3;
    foo.Clear(); foo.Start();
    tree.accAll(theta, eps);
    foo.Stop(); pprint("Barnes-Hut time: %e  rate: %e\n", foo.Elapsed(), (double)N/foo.Elapsed() );

    foo.Clear(); foo.Start();
    SmallVec<double, D> *anacc = new SmallVec<double, D>[N];
    for(int i=0; i<N; i++) anacc[i] = analyticAcc(i, theta, eps, p, N);
    foo.Stop(); pprint("      full time: %e  rate: %e\n", foo.Elapsed(), (double)N/foo.Elapsed() );

    double relerr[N];
    for(int i=0; i<N; i++) {
        dvec err = (anacc[i] - p[i].acc)/anacc[i].norm();
        relerr[i] = 0;
        for(int d=0; d<D; d++) relerr[i] = std::max(relerr[i], fabs(err[d]));
    }
    std::sort( &(relerr[0]), &(relerr[N]) );
    pprint("   max error: %e\n", relerr[N-1]);
    int n = 0.95*N;
    pprint("  0.95 error: %e\n", relerr[n]);
    n = 0.5*N;
    pprint("median error: %e\n", relerr[n]);


}
