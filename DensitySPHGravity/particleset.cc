#ifndef __PARTICLESET_INCLUDE__
#define __PARTICLESET_INCLUDE__

#include "smallvec.cc"
#include "ptype.h"

template <int D>
struct Particleset {

    Particleset() : psize( sizeof(ptype<D>) ), ps(NULL) {}

    Particleset(size_t N) : psize( sizeof(ptype<D>) ), np(N) {
        ps = new ptype<D>[np];
    }

    ~Particleset() {
        if( ps!=NULL ) delete[] ps;
    }

    void alloc(size_t N) {
        np = N;
        if( ps != NULL) delete[] ps;
        ps = new ptype<D>[np];
    }

    ptype<D> operator[] (const size_t i) const {
        assert(i < np);
        return *(ps + i);
    }

    ptype<D>& operator[] (const size_t i) {
        assert(i < np);
        return *(ps + i);
    }

    size_t psize;
    size_t np;
    ptype<D> *ps;
};

#endif // __PARTICLESET_INCLUDE__

#if 0
int main() {

    Particleset<3> pset(10);
    ptype<3> foo = pset[4];

}
#endif
