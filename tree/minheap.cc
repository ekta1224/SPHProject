#ifndef __MINHEAP_INCLUDE__
#define __MINHEAP_INCLUDE__

#include <algorithm>
#include <limits>

template <class T>
class PointDist {
public:
    PointDist() {}
    PointDist(int id, T dist) : id(id), dist(dist) {}
    inline T val(void) { return dist; }

    // needed for std::sort
    bool operator() (const PointDist &a, const PointDist &b) { return a.dist < b.dist; }

    int id;
    T dist;
};

template <class T>
class MinHeap {
public:
    MinHeap() : maxsize(0) {}
    ~MinHeap() { if(maxsize!=0) delete[] pq; }

    void alloc(int _maxsize) {
        if(_maxsize!=maxsize) {
            if(maxsize != 0) delete[] pq;
            pq = new PointDist<T>[_maxsize+1];
            maxsize = _maxsize;
        }
    }

    inline void clear(void) { maxv = -std::numeric_limits<T>::max(); heaplength = 0; }

    void removemax(void)  {
        exch(1, heaplength);
        sink(1, heaplength-1);
        heaplength--;
    }

    void insert(PointDist<T> v) {
        if(heaplength>=maxsize) {
            if(v.dist<maxv) {
                removemax();
                pq[++heaplength] = v;
                swim(heaplength);
                maxv = pq[1].dist;
            }
        }
        else {
            if(v.dist>maxv) maxv = v.dist;
            pq[++heaplength] = v;
            swim(heaplength);
        }
    }

    int length(void) { return heaplength; }
    int isfull(void) { return heaplength == maxsize; }
    T maxval(void) { return maxv; }

    PointDist<T> *pq;

private:
    PointDist<T> tempelement;

    inline void exch(int i, int j) { tempelement = pq[i]; pq[i] = pq[j]; pq[j] = tempelement; }

    void sink(int k, int N)  {
        int j;
        while (2*k <= N) {
            j = 2*k;
            if (j < N && pq[j].val() < pq[j+1].val() ) j++;
            if ( !(pq[k].val() < pq[j].val() ) ) break;
            exch(k, j);
            k = j;
        }
    }

    void swim(int k) {
        while (k > 1 && ( pq[k/2].val() < pq[k].val() ) )  {
            exch(k, k/2);
            k = k/2;
        }
    }

    int maxsize;
    T maxv;
    int heaplength;
};

#endif // __MINHEAP_INCLUDE

#if 0
#include <cstdioh>
#include <cstdlib>
#include <cmath>

int main(void) {

    MinHeap<double> heap;

    heap.alloc(5);
    heap.clear();
    for(int i=1; i<100; i++) {
        double dist = i;
        heap.insert(PointDist<double>(i, dist) );
    }


    printf("%d\n", heap.length());
    printf("%e\n", heap.maxval());

#if 1
    for(int i=1; i<=heap.length(); i++) {
        printf("%d  %f\n", i-1, heap.pq[i].dist);
    }
#endif
}
#endif
