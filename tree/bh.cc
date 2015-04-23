/*
Time-stamp: <bh.cc on Thursday, 12 March, 2015 at 10:40:43 MST (pinto)>

  basic "octree" in D = {1,2,3} dimensions

  with current version of MortonKey class, can do up to 32 levels
  i.e. each particle coordinate in the Morton key is represented by an unsigned 32 bit integer
  (about 9 digits of precision in scaled coordinates)

*/

#ifndef __OCTREE_INCLUDE__
#define __OCTREE_INCLUDE__

#include "morton.cc"

template <class T>
inline T SQR(T x) {return x*x;}

#define FORALLDAUGHTERS(D,node)   for(int D=tree[node].fd; D<tree[node].fd+tree[node].dcnt; D++)
#define FORALLPOINTS(I,node)      for(int I=tree[node].begin; I<=tree[node].end; I++)

#ifdef SORTED
#define pm(i) ps[i]
#else
#define pm(l) ps[this->pindex(morton[l])]
#endif

template <int D>
class Octree : public MortonKey<D> {
public:
private:
    typedef SmallVec<int, D> ivec;
    typedef SmallVec<double, D> dvec;

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

    struct bhdata {
        dvec com;       // center of mass
        double mass;   // mass of node
        double Bmax;   // maximum distance of particle from com
    };

public:
    Octree();
    ~Octree();
    void initTree(const dvec center, const double halfwidth);
    void buildTree(ptype<D> *const ps, int np, int mls);
    void checkTree();

    int isLeaf(const int node)      { assert(node < n_nodes); return (tree[node].dcnt == 0); }
    int nodeSize(const int node)    { assert(node < n_nodes); return tree[node].end - tree[node].begin + 1; }
    dvec nodeCenter(const int node) { assert(node < n_nodes); return tree[node].center; }

    int whereAmI(dvec pos);
    void nodeContents(int node, ptype<D> *const P);
    void replaceNodeContents(int node, ptype<D> *const P);

    void propagateCOM(int node);
    void propagateBmax(int node);
    void direct(int sink, dvec source, double srcmass);
    void acc1(int p, double theta, double eps);
    void accAll(double theta, double eps);
    void bhSingle(int p, double theta, int node);


    treenode *tree;
    int n_leafSet, *leafSet; // list of leaf nodes
    int n_nodes, *node2leaf; // map from leafSet index into tree node

    bhdata *BH;
    double eps2;

private:
    void MakeMorton();
    void IOT1(int level, int b, int e, int parent);
    void IOT2(int level, int b, int e, int parent, dvec center, double hw);
    void sortPoints();

    int within(const dvec p, dvec nodecenter, double hw);
    void reportWithin(dvec pos, dvec nodecenter, double hw);
    void traverseCheckTree(int node, int level, dvec nodecenter, double hw);

    int internalWhereAmI(typename MortonKey<D>::mortonkey mymorton, int node, int level);
    int sphereNodeIntersect(dvec sinkcenter, double r2, int node);

    static const int MAXLEVEL = 32;
    typename MortonKey<D>::mortonkey *morton;

    int np;
    ptype<D> *ps;

    int maxleafsize;
    int maxdaughter; // maximum number of daughters (2, 4, and 8) for D of (1, 2, 3)
    dvec ROOTcenter;
    double halfwidth;
    int ROOT, ROOTlevel, nodeptr;
    int levelcount[MAXLEVEL], levelptr[MAXLEVEL];
};

template <int D>
Octree<D>::Octree() {
    maxdaughter = (1<<D) - 1;
    ROOT = 0;    ROOTlevel = 1;
    tree = NULL; morton = NULL; leafSet = NULL; node2leaf = NULL;
    BH = NULL;
}

template <int D>
void Octree<D>::initTree(dvec _center, double _halfwidth) {
    halfwidth = _halfwidth;
    ROOTcenter = _center;
    this->initMortonKey(_center-dvec(halfwidth), _center+dvec(halfwidth));
}

template <int D>
Octree<D>::~Octree() {
    if(      tree != NULL ) delete[] tree;       tree = NULL;
    if(        BH != NULL ) delete[] BH;         BH = NULL;
    if(    morton != NULL ) delete[] morton;     morton = NULL;
    if(   leafSet != NULL ) delete[] leafSet;    leafSet = NULL;
    if( node2leaf != NULL ) delete[] node2leaf;  node2leaf = NULL;
}

// sort points into morton order
template <int D>
void Octree<D>::sortPoints() {
    ptype<D> *tmp = new ptype<D>[np];
    for(int i=0;i<np;i++) tmp[i] = ps[ this->pindex(morton[i]) ];
    std::memcpy( ps, tmp, sizeof(ptype<D>)*np );
    delete[] tmp;
}

template <int D>
void Octree<D>::buildTree(ptype<D> *const _ps, int _np, int _maxleafsize) {

    ps = _ps;
    np = _np;
    maxleafsize = _maxleafsize;

    morton = new typename MortonKey<D>::mortonkey[np];
    for(int p=0; p<np; p++) morton[p] = this->Morton(ps[p].r, p);
    std::sort(&(morton[0]), &(morton[np]), AscendingMorton());
    for(int p=0; p<np-1; p++) assert( AM(morton[p], morton[p+1]) );

#ifdef SORTED
    sortPoints();
#endif

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
    BH = new bhdata[n_nodes];

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
}

// first pass to count cells and leaves
template <int D>
void Octree<D>::IOT1(int level, int b, int e, int parent) {

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
void Octree<D>::IOT2(int level, int b, int e, int parent, dvec center, double hw) {

    assert( level < MAXLEVEL ); assert( e>b );
    // shouldn't be here if parent has any daughters
    assert( tree[parent].fd == 0 ); assert( tree[parent].dcnt == 0 );

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

//------------------------------------------------------------------------------------------------------------------------
// baseTree checking code
//------------------------------------------------------------------------------------------------------------------------

// is a particle within the node
template <int D>
int Octree<D>::within(dvec pos, dvec nodecenter, double hw) {
    int success = 1;
    for(int d=0; d<D; d++)
        if( !(pos[d] >= nodecenter[d] - hw && pos[d] <= nodecenter[d] + hw) ) success = 0;
    return success;
}

// format complaint
template <int D>
void Octree<D>::reportWithin(dvec pos, dvec nodecenter, double hw) {
    std::string name[3] = {"x", "y", "z"};
    pprint("nodecenter: % .15e  hw: % .15e\n", nodecenter, hw);
    for(int d=0; d<D; d++)
        pprint("%s:  % .15e <= % .15e <= % .15e\n",
               name[d], nodecenter[d]-hw, pos[d], nodecenter[d]+hw);
}

// traverse the tree checking that each point is within its purported node and
// the levels, widths, and centers are correct
template <int D>
void Octree<D>::traverseCheckTree(int node, int level, dvec nodecenter, double hw) {

    assert( tree[node].level == level );
    assert( tree[node].hw == hw );

    FORALLPOINTS(p, node) {
        if( !within( pm(p).r, nodecenter, hw) ) reportWithin(pm(p).r, nodecenter, hw);
        assert( within( pm(p).r, nodecenter, hw) );
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
    pprint("checkTree passed\n");
}

//------------------------------------------------------------------------------------------------------------------------
// Accessing points
//------------------------------------------------------------------------------------------------------------------------

// copy the particle data of node into array P (assumed to be allocated with sufficient space)
template <int D>
void Octree<D>::nodeContents(int node, ptype<D> *const P) {
    assert(node < n_nodes );
#ifdef SORTED
    std::memcpy( P,  &(ps[tree[node].begin]), sizeof(ptype<D>) * nodeSize(node) );
#else
    FORALLPOINTS(i, node)   P[i-BEGIN(node)] = pm(i);
#endif
}

// copy the particle data P into the particles of node; P is assumed to contain the correct number of points
template <int D>
void Octree<D>::replaceNodeContents(int node, ptype<D> *const P) {
    assert( node < n_nodes );
#ifdef SORTED
    std::memcpy(&(ps[tree[node].begin]), P, sizeof(ptype<D>)*nodeSize(node) );
#else
    int j = 0;
    FORALLPOINTS(i, node) { pm(i) = P[j]; j++; }
#endif
}

//------------------------------------------------------------------------------------------------------------------------
// utility geometry functions
//------------------------------------------------------------------------------------------------------------------------

// given a MortonKey, determine the node number of the leaf cell where this key belongs
template <int D>
int Octree<D>::internalWhereAmI(typename MortonKey<D>::mortonkey mymorton, int node, int level) {

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

// return the index of the node containing pos,  giving an error if no leaf node contains this pos
template <int D>
int Octree<D>::whereAmI(dvec pos) {
    typename MortonKey<D>::mortonkey mymorton = this->Morton(pos, 0);
    return internalWhereAmI(mymorton, ROOT, ROOTlevel);
}

// sphere-cube intersection; returns true if any part of node is within the sphere
template <int D>
int Octree<D>::sphereNodeIntersect(dvec center, double r2, int node) {
    assert(node < n_nodes);
    dvec nodecenter = tree[node].center;

    // fast rejection: sphere-sphere intersection
    double r = sqrt(r2);
    if( (center-nodecenter).norm2() > SQR( sqrt(3)*tree[node].hw + r ) ) return 0;

    // they might intersect -- do the more expensive exact test
    double hw = tree[node].hw;
    dvec nmin = nodecenter - dvec(hw);
    dvec nmax = nodecenter + dvec(hw);

    double mindist2 = 0;
    for(int d=0; d<D; d++) {
        if(      center[d] < nmin[d] ) mindist2 += SQR( center[d] - nmin[d] );
        else if( center[d] > nmax[d] ) mindist2 += SQR( center[d] - nmax[d] );
    }
    return mindist2 <= r2;
}

//------------------------------------------------------------------------------------------------------------------------
// Barnes-Hut acceleration
//------------------------------------------------------------------------------------------------------------------------

template <int D>
void Octree<D>::propagateCOM(int node) {
    assert(!isLeaf(node));

    FORALLDAUGHTERS(d,node) {
        if(isLeaf(d)) {
            FORALLPOINTS(i, d) {
                BH[node].com += pm(i).m * pm(i).r;
                BH[node].mass += pm(i).m;
            }
        }
        else {
            propagateCOM(d);
            BH[node].com  += BH[d].mass * BH[d].com;
            BH[node].mass += BH[d].mass;
        }
    }
    BH[node].com /= BH[node].mass;
}

template <int D>
void Octree<D>::propagateBmax(int node) {
    assert(!isLeaf(node));

    double bmax2 = 0;
    BH[node].Bmax = 0;
    FORALLPOINTS(i, node) bmax2 = std::max( bmax2, (pm(i).r - BH[node].com).norm2() );
    BH[node].Bmax = sqrt(bmax2);

    FORALLDAUGHTERS(d,node) if(!isLeaf(d)) propagateBmax(d);
}

template <int D>
inline void Octree<D>::direct(int sink, dvec source, double srcmass) {
    dvec dr = pm(sink).r - source;
    double nr = sqrt(dr.norm2() + eps2);
    double ir3 = srcmass/(nr*nr*nr);
    pm(sink).acc -= ir3*dr;
}

template <int D>
void Octree<D>::acc1(int p, double theta, double eps) {
    eps2 = eps*eps;
    bhSingle(p, theta, ROOT);
}

template <int D>
void Octree<D>::accAll(double theta, double eps) {
    eps2 = eps*eps;
    for(int p=0; p<np; p++) bhSingle(p, theta, ROOT);
}

template <int D>
void Octree<D>::bhSingle(int p, double theta, int node) {
    if(isLeaf(node)) {
        FORALLPOINTS(i, node) if(i != p) direct(p, pm(i).r, pm(i).m);
    }
    else {
        double dr = (pm(p).r - BH[node].com).norm();
        if( BH[node].Bmax < theta * dr )
            direct(p, BH[node].com, BH[node].mass);
        else
            FORALLDAUGHTERS(d,node) bhSingle(p, theta, d);
    }
}

#undef pm
#undef FORALLDAUGHTERS
#undef FORALLPOINTS

#endif // __OCTREE_INCLUDE__
