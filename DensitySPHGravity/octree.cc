/*
Time-stamp: <octree.cc on Thursday, 9 April, 2015 at 09:55:41 MST (pinto)>

  basic "octree" in D = {1,2,3} dimensions

  with current version of MortonKey class, can do up to 32 levels
  i.e. each particle coordinate in the Morton key is represented by an unsigned 32 bit integer
  (about 9 digits of precision in scaled coordinates)

*/

#ifndef __OCTREE_INCLUDE__
#define __OCTREE_INCLUDE__

#include "morton.cc"

// source list data type
template <int D>
struct xyzm {
    SmallVec<double, D> pos; // 24
    double m;                // 32
};

#define FORALLDAUGHTERS(D,node)   for(int D=tree[node].fd; D<tree[node].fd+tree[node].dcnt; D++)
#define FORALLPOINTS(I,node)      for(int I=tree[node].begin; I<=tree[node].end; I++)

template <int D>
class Octree : public MortonKey<D> {
public:
private:
    typedef SmallVec<int, D> ivec;
    typedef SmallVec<double, D> dvec;

    // internal point data type
    struct otpoint {
        SmallVec<double, D> xyz;  // 24
        double m;                 // 32
        int32_t id;                // 36 (probably want size_t?)
    };

    struct treenode {
        dvec center;     // center coordinate of cell
        double hw;       // half-width of the cell
        int begin, end;  // beginning and ending particle number (inclusive) in node
        int parent;      // node of parent
        uint32 fd;       // node of first daughter
        int16_t wd;      // daughter direction from parent
        int8_t dcnt;     // how many daughters (0 to 8)
        int8_t level;    // level in tree
    };

public:

    Octree(int const maxthreads);
    ~Octree();
    void reset();
    void buildTree(dvec const _center, double const _halfwidth, Particleset<D> const &pset, int const _maxleafsize);
    void checkTree();
    void knn(const dvec pos, const int knn, double &r_knn, int *ids);
    void ballSearch(const dvec center, const double radius, int *ids, int &nids);

//private:

    int isLeaf(int const node)      { assert(node < n_nodes); return (tree[node].dcnt == 0); }
    int nodeSize(int const node)    { assert(node < n_nodes); return tree[node].end - tree[node].begin + 1; }
    dvec nodeCentre(int const node) { assert(node < n_nodes); return tree[node].center; }
    int isIn(dvec const pos, int const node) { return within(pos, tree[node].center, tree[node].hw); }

    void sortPoints();
    void IOT1(int const level, int const b, int const e, int const parent);
    void IOT2(int const level, int const b, int const e, int const parent, dvec const center, double const hw);
    int within(dvec const pos, dvec const nodecenter, double const hw);
    void reportWithin(dvec const pos, dvec const nodecenter, double const hw);
    void traverseCheckTree(int const node, int const level, dvec const nodecenter, double const hw);
    void nodeContents(int const node, xyzm<D> const *points);
    void nodeIds(int const node, int const *ids);
    int internalWhereAmI(typename MortonKey<D>::mortonkey const mymorton, int const node, int const level);
    int whereAmI(dvec const pos);
    int sphereNodeIntersect(dvec const center, double const r2, int const node);
    void traverseKNN(const dvec pos, int node, const double dist2, MinHeap<double> &heap);
    void internalBallSearch(const dvec pos, int node, const double r2, int *ids, int &nids);

    static const int MAXLEVEL = 32;
    typename MortonKey<D>::mortonkey *morton;

    int n_leafSet, *leafSet; // list of leaf nodes
    int n_nodes, *node2leaf; // map from leafSet index into tree node

    int np;
    otpoint *ps;
    treenode *tree;

    int maxthreads;
    MinHeap<double> *heaps;

    int maxleafsize;
    int maxdaughter; // maximum number of daughters (2, 4, and 8) for D of (1, 2, 3)
    dvec ROOTcenter;
    double halfwidth;
    int ROOT, ROOTlevel, nodeptr;
    int levelcount[MAXLEVEL], levelptr[MAXLEVEL];
};


template <int D>
Octree<D>::Octree(int const maxthreads)
    : maxthreads(maxthreads), tree(NULL), morton(NULL), leafSet(NULL), node2leaf(NULL), ps(NULL) {
    maxdaughter = (1<<D) - 1;
    ROOT = 0;    ROOTlevel = 1;
    heaps = new MinHeap<double>[maxthreads];
}

template <int D>
Octree<D>::~Octree() {
    reset();
    delete[] heaps;
}

template <int D>
void Octree<D>::reset() {
    np = 0;
    ps = NULL; tree = NULL; morton = NULL; leafSet = NULL; node2leaf = NULL;
    if(      tree != NULL ) delete[] tree;       tree = NULL;
    if(    morton != NULL ) delete[] morton;     morton = NULL;
    if(   leafSet != NULL ) delete[] leafSet;    leafSet = NULL;
    if( node2leaf != NULL ) delete[] node2leaf;  node2leaf = NULL;
    if(        ps != NULL ) delete[] ps;         ps = NULL;
}

template <int D>
void Octree<D>::buildTree(dvec const _center, double const _halfwidth, Particleset<D> const &pset, int const _maxleafsize) {

    reset();

    halfwidth = _halfwidth;
    ROOTcenter = _center;
    maxleafsize = _maxleafsize;

    np = pset.np;
    ps = new otpoint[np];
    for(int i=0; i<np; i++) {
        ps[i].xyz = pset[i].r;
        ps[i].m = pset[i].m;
        ps[i].id = i;        // original index of particle (save before sorting ps in morton order)
    }

    // compute Morton keys and sort ps in Morton order
    this->initMortonKey(_center-dvec(halfwidth), _center+dvec(halfwidth));
    morton = new typename MortonKey<D>::mortonkey[np];
    for(int p=0; p<np; p++) morton[p] = this->Morton(ps[p].xyz, p);
    std::sort(&(morton[0]), &(morton[np]), AscendingMorton());
    for(int p=0; p<np-1; p++) assert( AM(morton[p], morton[p+1]) );
    sortPoints();

    n_leafSet = 0;
    levelcount[0] = 1;
    for(int l=1; l<MAXLEVEL; l++) levelcount[l] = 0;

    // run through the data once to determine space required for tree
    n_nodes = 1; // count the root!
    IOT1(ROOTlevel, 0, np-1, ROOT);

    int nc = 0;
    for(int l=0; l<MAXLEVEL; l++) nc += levelcount[l];
    assert( nc == n_nodes );

    // allocate space for leafSet and node2leaf index mapping
    node2leaf = new int[n_nodes];
    for(int i=0;i<n_nodes;i++) node2leaf[i] = -1;
    leafSet = new int[n_leafSet];
    // allocate tree
    tree = new treenode[n_nodes];

    levelptr[0] = 0;
    for(int l=1; l<MAXLEVEL; l++) levelptr[l] = levelptr[l-1] + levelcount[l-1];

    nodeptr = 0;
    tree[ROOT].fd = 0;            tree[ROOT].wd = 0;               tree[ROOT].dcnt = 0;
    tree[ROOT].level = ROOTlevel; tree[ROOT].center = ROOTcenter;  tree[ROOT].hw = halfwidth;
    tree[ROOT].begin = 0;         tree[ROOT].end = np-1;
    tree[ROOT].parent = -1;

    // run through the data a second time to build the tree
    n_leafSet = 0;
    IOT2(ROOTlevel, 0, np-1, ROOT, ROOTcenter, 0.5*halfwidth);
    assert( nodeptr+1 == n_nodes );

    delete[] morton;
}

template <int D>
void Octree<D>::sortPoints() {
    otpoint *tmp = new otpoint[np];
    for(int i=0;i<np;i++) tmp[i] = ps[ this->pindex(morton[i]) ];
    std::memcpy( ps, tmp, sizeof(otpoint)*np );
    delete[] tmp;
}

// first pass to count cells and leaves
template <int D>
void Octree<D>::IOT1(int const level, int b, int const e, int const parent) {

    int direction = this->key(morton[b],level);
    while( (direction<=maxdaughter) && (b<=e) ) {

        int count = 0;
        while( b <= e ) {
            if( this->key(morton[b],level) == direction ) { b++; count++; }
            else break;
        }
        assert(count>0);

        levelcount[level]++;
        n_nodes++;

        if( count <= maxleafsize )
            n_leafSet++;
        else
            IOT1(level+1, b-count, b-1, n_nodes);

        if(b<=e) direction = this->key(morton[b],level);
    }
}

// second pass to fill in tree data
template <int D>
void Octree<D>::IOT2(int const level, int b, int const e, int const parent, dvec const center, double const hw) {

    assert( level < MAXLEVEL ); assert( e>b );
    assert( tree[parent].fd == 0 ); assert( tree[parent].dcnt == 0 );  // shouldn't be here yet if parent
                                                                       // has any daughters

    int direction = this->key(morton[b],level);
    while( b<=e ) {
        assert( direction <= maxdaughter );

        // count number of particles in this octant
        int count = 0;
        while( b <= e ) {
            if( this->key(morton[b],level) == direction ) { b++; count++; }
            else break;
        }
        assert( count > 0 );

        // get daughter node number in tree
        int daughter = levelptr[level];
        levelptr[level]++;

        if( tree[parent].fd == 0 ) {
            // first daughter
            assert( tree[parent].dcnt == 0 );
            tree[parent].fd=daughter;
            tree[parent].dcnt = 1;
        }
        else {
            // subsequent daughters
            tree[parent].dcnt++;
            assert( tree[parent].dcnt <= maxdaughter+1 );
        }

        tree[daughter].level = level+1;
        tree[daughter].parent = parent;
        tree[daughter].begin = b-count;
        tree[daughter].end = b-1;
        tree[daughter].hw = hw;
        for(int d=0; d<D; d++) tree[daughter].center[d] = center[d] + hw*Mask<D>::mask[direction][d];
        tree[daughter].wd = direction;
        tree[daughter].fd = 0;
        tree[daughter].dcnt = 0;
        nodeptr++;
        assert( nodeptr < n_nodes );

        if( count <= maxleafsize ) {
            // node has <= maxleafsize particles -- make it a leaf
            leafSet[n_leafSet] = daughter;
            node2leaf[daughter] = n_leafSet;
            n_leafSet++;
        }
        else {
            // node has too many particles -- recurse to divide into daughters
            IOT2(tree[daughter].level, b-count, b-1, daughter, tree[daughter].center, 0.5*hw);
        }

        // more points left, get next octant to process
        if(b<=e) direction = this->key(morton[b],level);
    }
}

//---------------------------------------------------------------------------------------------------------------------
// Tree checking code
//---------------------------------------------------------------------------------------------------------------------

// is a particle within the node
template <int D>
int Octree<D>::within(dvec const pos, dvec const nodecenter, double const hw) {
    int success = 1;
    for(int d=0; d<D; d++)
        if( !(pos[d] >= nodecenter[d] - hw && pos[d] <= nodecenter[d] + hw) ) success = 0;
    return success;
}

// format complaint
template <int D>
void Octree<D>::reportWithin(dvec const pos, dvec const nodecenter, double const hw) {
    std::string name[3] = {"x", "y", "z"};
    pprint("nodecenter: % .15e  hw: % .15e\n", nodecenter, hw);
    for(int d=0; d<D; d++)
        pprint("%s:  % .15e <= % .15e <= % .15e\n",
               name[d], nodecenter[d]-hw, pos[d], nodecenter[d]+hw);
}

// traverse the tree checking that each point is within its purported node and
// the levels, widths, and centers are correct
template <int D>
void Octree<D>::traverseCheckTree(int const node, int const level, dvec const nodecenter, double const hw) {

    assert( tree[node].level == level );
    assert( tree[node].hw == hw );

    FORALLPOINTS(p, node) {
        if( !within( ps[p].xyz, nodecenter, hw) ) reportWithin(ps[p].xyz, nodecenter, hw);
        assert( within( ps[p].xyz, nodecenter, hw) );
    }

    FORALLDAUGHTERS(daughter, node) {
        dvec dcenter = nodecenter;
        for(int d=0; d<D; d++) dcenter[d] += 0.5*hw*Mask<D>::mask[tree[daughter].wd][d];
        dvec dc = tree[daughter].center;
        assert(dc == dcenter );
        traverseCheckTree(daughter, level+1, dcenter, 0.5*hw);
    }
}

// public interface for checking tree
template <int D>
void Octree<D>::checkTree() {
    traverseCheckTree(ROOT, ROOTlevel, ROOTcenter, halfwidth);
    fpprint(std::cerr,"checkTree passed\n");
}

//---------------------------------------------------------------------------------------------------------------------
// Accessing points
//---------------------------------------------------------------------------------------------------------------------

// copy the particle data of node into array xyzm<D> *P
//    assumed to be allocated with sufficient space, so that length is already know from
//    a call to nodeSize
template <int D>
void Octree<D>::nodeContents(int const node, xyzm<D> const *points) {
    assert(node < n_nodes );
    int m = 0;
    FORALLPOINTS(i, node) {
        points[i].pos = ps[i].xyz;
        points[i].m = ps[i].m;
    }
}

// copy the particle data of node into array xyzm<D> *P
//    assumed to be allocated with sufficient space, so that length is already know from
//    a call to nodeSize
template <int D>
void Octree<D>::nodeIds(int const node, int const *ids) {
    assert(node < n_nodes );
    int m = 0;
    FORALLPOINTS(i, node) {
        ids[i] = ps[i].id;
    }
}

//---------------------------------------------------------------------------------------------------------------------
// utility geometry functions
//---------------------------------------------------------------------------------------------------------------------

// given a MortonKey, determine the node number of the leaf cell where this key belongs
template <int D>
int Octree<D>::internalWhereAmI(typename MortonKey<D>::mortonkey const mymorton, int const node, int const level) {

    // if a leaf, terminate with my node
    if(isLeaf(node)) return node;

    // else recurse into the correct direction
    int direction = this->key(mymorton,level);
    int daughter = -1;
    FORALLDAUGHTERS(d, node) {
        if( tree[d].wd == direction ) {
            daughter = internalWhereAmI(mymorton, d, level+1);
            break;
        }
    }
    if( daughter == -1 ) {
        std::cerr << "internalWhereAmI: leaf cell doesn't exist in tree\n";
        assert( daughter >= 0 );
    }

    return daughter;
}

// return the index of the node containing pos
//    giving an error if no leaf node contains this pos
template <int D>
int Octree<D>::whereAmI(dvec const pos) {
    typename MortonKey<D>::mortonkey mymorton = this->Morton(pos, 0);
    return internalWhereAmI(mymorton, ROOT, ROOTlevel);
}

// sphere-cube intersection
//   returns true if any part of node is within the sphere
template <int D>
int Octree<D>::sphereNodeIntersect(dvec const center, double const r2, int const node) {
    assert(node < n_nodes);
    dvec nodecenter = tree[node].center;

    // fast rejection: sphere-sphere intersection
    double r = sqrt(r2);
    if( (center-nodecenter).norm2() > sqr( sqrt(3)*tree[node].hw + r ) ) return 0;

    // they might intersect -- do the more expensive exact test
    double hw = tree[node].hw;
    dvec nmin = nodecenter - dvec(hw);
    dvec nmax = nodecenter + dvec(hw);

    double mindist2 = 0;
    for(int d=0; d<D; d++) {
        if(      center[d] < nmin[d] ) mindist2 += sqr( center[d] - nmin[d] );
        else if( center[d] > nmax[d] ) mindist2 += sqr( center[d] - nmax[d] );
    }
    return mindist2 <= r2;
}

//--------------------------------------------------------------------------------------------------------------------
// K-nearest neighbours
//--------------------------------------------------------------------------------------------------------------------

template <int D>
void Octree<D>::traverseKNN(const dvec pos, int node, const double dist2, MinHeap<double> &heap) {

    if( sphereNodeIntersect(pos, dist2, node) ) {
        if(isLeaf(node))
            FORALLPOINTS(i, node) {
                double dr2 = (pos - ps[i].xyz).norm2();
                heap.add(i, dr2); // adds the id OF THE MORTON-ORDER PARTICLE
            }
        else
            FORALLDAUGHTERS(d, node) traverseKNN(pos, d, dist2, heap);
    }
}

// find knn nearest neighbours (including self) to the given postion
//    return the original particle ids in ids
//    (which must be preallocated of length >= knn)
template <int D>
void Octree<D>::knn(const dvec pos, const int knn, double &r_knn, int *ids) {
    assert( knn <= np );

    int g = omp_get_thread_num();
    heaps[g].alloc(knn);

    // get an approximation to knn distance from local node
    int node = whereAmI(pos);
    if( nodeSize(node) == 1 ) node = tree[node].parent;
    double maxdr2 = 0;
    FORALLPOINTS(i, node) maxdr2 = std::max(maxdr2, (pos - ps[i].xyz).norm2() );

    // try with this distance; if heap cannot be filled, increase distance
    // until success
    do {
        heaps[g].clear();
        traverseKNN(pos, ROOT, maxdr2, heaps[g]);
        maxdr2 *= 2;
    } while(!heaps[g].isfull() );// && maxdr2<3*halfwidth*halfwidth);

    // once more, with a will...
    maxdr2 = heaps[g].maxval();
    heaps[g].clear();
    traverseKNN(pos, ROOT, maxdr2, heaps[g]);

    // convert morton-order ids to original ids
    for(int i=0; i<heaps[g].length(); i++) {
        int j = heaps[g][i].id;
        ids[i] = ps[j].id;
    }
    r_knn = sqrt(heaps[g].maxval());
}

template <int D>
void Octree<D>::internalBallSearch(const dvec pos, int node, const double r2, int *ids, int &nids) {
    if( sphereNodeIntersect(pos, r2, node) ) {
        if(isLeaf(node)) {
            FORALLPOINTS(i, node)
                if( (pos - ps[i].xyz).norm2() <= r2 ) ids[nids++] = i;
        }
        else
            FORALLDAUGHTERS(d, node) internalBallSearch(pos, d, r2, ids, nids);
    }
}

// return a list of the original particle ids within radius of center
template <int D>
void Octree<D>::ballSearch(const dvec center, const double radius, int *ids, int &nids) {
    double r2 = radius*radius;
    assert( r2 > 0 );
    nids = 0;
    internalBallSearch(center, ROOT, r2, ids, nids);
    for(int i=0; i<nids; i++) {
        int j = ids[i];
        ids[i] = ps[j].id;
    }
}

#undef FORALLDAUGHTERS
#undef FORALLPOINTS

#endif // __OCTREE_INCLUDE__
