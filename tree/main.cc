/*
Time-stamp: <main.cc on Thursday, 12 March, 2015 at 10:15:31 MST (pinto)>

g++ -DSORTED -std=c++11 -O3 -Ilib -fopenmp main.cc

 */


#include <cassert>
//#include <omp.h>
#include <cstring>
#include <cmath>
#include "pprint.cc"
#include "smallvec.cc"
#include "STimer.cc"
#include "minheap.cc"

template <int D>
struct ptype {
    SmallVec<double, D> r;
    int id;
};

#include "octree.cc"

template <int D>
void mkpoints(SmallVec<double, D> lowerleft, SmallVec<double, D> upperright, int N,
              ptype<D> *p) {
    pprint("making %d uniform particles on %e to %e\n", N, lowerleft, upperright);

    SmallVec<double, D> delta = upperright - lowerleft;
    for(int i=0; i<N; i++) {
        for(int d=0; d<D; d++) {
            double coord;
            do {
                coord = lowerleft[d] + drand48()*delta[d];
            } while( !(coord>=lowerleft[d] && coord<upperright[d]) );
            p[i].r[d] = coord;
        }
    }
}


int main() {

    const int D = 3;
    typedef SmallVec<double, D> dvec;

    int KNN = 32;

    dvec center(0);
    double halfwidth = 1;
    int N = 1<<22;
    ptype<D> *p = new ptype<D>[N];
    mkpoints(center-dvec(halfwidth), center+dvec(halfwidth), N, p);

    Octree<D> tree;
    tree.initTree(center, halfwidth);
    STimer foo; foo.Clear(); foo.Start();
    tree.buildTree(p, N, 2*KNN);
    foo.Stop(); pprint("time to build tree: %e  rate: %e\n", foo.Elapsed(), (double)N/foo.Elapsed());
    tree.checkTree();
    pprint("nnodes: %d\n", tree.n_nodes);
    pprint("nleaves: %d\n", tree.n_leafSet);
}
