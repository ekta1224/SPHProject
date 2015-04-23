#include <omp.h>
#include "smallvec.cc"
#include "pprint.cc"
#include "util.cc"

#include "ptype.h"
#include "particleset.cc"
#include "getparticleshalogen.cc"

#include "minheap.cc"
#include <cstring>
#include "octree.cc"

void ensureKNN(const int knn, int *ids, const SmallVec<double, 3> center,
               const double r_max, Particleset<3> &part) {
    for(int i=0; i<knn; i++) {
        assert( (part[ids[i]].r - center).norm() <= r_max );
    }
}

template <int D>
void testknn(Octree<D> &OT, const int knn, Particleset<D> &part) {
    int nprocs = omp_get_num_procs();
    int ids[nprocs][knn];
#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<1000; i++) {
        int g = omp_get_thread_num();
        int j = drand48() * part.np;
        double r_max;
        OT.knn(part[j].r, knn, r_max, ids[g]);
        assert( r_max > 0 );
        ensureKNN( knn, ids[g], part[j].r, r_max, part);
    }
}

template <int D>
void testBallSearch(Octree<D> &OT, Particleset<D> &part, const double rmax) {

    int nprocs = omp_get_num_procs();
    int ids[nprocs*part.np];
#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<1000; i++) {
        int g = omp_get_thread_num();
        double radius = drand48() * rmax;
        int j = drand48() * part.np;
        int nids;
        OT.ballSearch(part[j].r, radius, &(ids[g*part.np]), nids);
    }
}

template <int D>
void testTree(Octree<D> &OT, Particleset<D> &part) {

    int nn = 0, nl = 0;
    int sum;
    for(int i=0; i<OT.n_nodes; i++) {
        if(OT.isLeaf(i)) nl++; else nn++;
        sum += OT.nodeSize(i);
        SmallVec<double, D> cen = OT.nodeCentre(i);
    }

    for(int i=0; i<part.np; i++) {
        int j = OT.whereAmI(part[i].r);
        assert( OT.isIn(part[i].r, j) );
    }


}

int main() {
    const int D = 3;
    int nprocs = omp_get_num_procs();

    int maxleafsize = 1;


    Particleset<D> part;
    double halfwidth;
    SmallVec<double, D> center;
    getParticlesHalogen("testparticles", part, center, halfwidth);
    center.zero();

    Octree<D> OT(omp_get_num_procs());
    OT.buildTree(center, halfwidth, part, maxleafsize);
    OT.checkTree();

    testTree(OT, part);
    pprint("done testTree\n");

    // doing knn searches
    testknn(OT, 20, part);
    testknn(OT, 40, part);
    testknn(OT, 60, part);
    pprint("done testknn\n");

    // doing ball searches
    testBallSearch(OT, part, halfwidth);
    pprint("done testBallSearche\n");

}
