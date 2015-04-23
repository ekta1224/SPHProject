/*
Time-stamp: <minheap.cc on Tuesday, 7 April, 2015 at 19:42:57 MST (pinto)>

  Minheap:
     maintains a list of the N smallest values added to the heap

  Usage:

  minheap<double> heap;
  heap.alloc(heap_length);
  heap.clear();

  many times: heap.add(index, value)

  heap.maxval() is largest of heap_length minimum values
  heap[i].id is the array (zero-offset) of heap_length indices
  heap[i].val is the array of values in the heap
  heap.removemax() pops the largest value off of the heap

*/

#ifndef __MINHEAP_INCLUDE__
#define __MINHEAP_INCLUDE__

#include <utility>
#include <limits>

namespace heapspace {
    template <class T>
    struct Pair {
        Pair() {}
        Pair(int id, T value) : id(id), val(value) {}
        int id;
        T val;
    };
};

template <class T>
class MinHeapBase {
public:
    MinHeapBase() : maxsize(0) {}
    ~MinHeapBase() { if(maxsize!=0) delete[] list; }

    void alloc(int _maxsize) {
        if(_maxsize!=maxsize) {
            if(maxsize != 0) delete[] list;
            list = new heapspace::Pair<T>[_maxsize+1];
            maxsize = _maxsize;
        }
    }

    inline void clear(void) { maxv = -std::numeric_limits<T>::max(); heaplength = 0; }

    T removemax(void)  {
        // remove top element and return it
        double oldmaxv = maxv;
        // replace top element with end of heap
        // then call sink to bubble down the value into its correct place in the tree
        std::swap(list[1], list[heaplength]);
        sink(1, heaplength-1);
        heaplength--;
        maxv = list[1].val;
        return oldmaxv;
    }

    void print() {
        for(int i=1; i<=heaplength; i++) pprint("(%d, %2.2e) ",list[i].id, list[i].val);
        pprint("\n");
    }

    void insert(heapspace::Pair<T> v) {
        // place new element in first free position in array
        //    if heap is full, add only if smaller than maximum value
        // then call swim to bubble up the value to the correct place in tree
        if(heaplength>=maxsize) {
            if(v.val<maxv) {
                removemax();
                list[++heaplength] = v;
                swim(heaplength);
                maxv = list[1].val;
            }
        }
        else {
            if(v.val>maxv) maxv = v.val;
            list[++heaplength] = v;
            swim(heaplength);
        }
    }

    int length(void) { return heaplength; }
    int isfull(void) { return heaplength == maxsize; }
    T maxval(void) { return maxv; }

    // give access to heap with 0-offset indexing (list is on [1,heaplength])
    inline heapspace::Pair<T> operator[] (int i) {
        return list[i+1];
    }

    heapspace::Pair<T> *list;

private:
    // minheap invarients:
    //   value in node <= value in node's parents
    //   binary tree is complete (every level is completely full but the last)
    //
    // addressing: for a node i:
    //      left child is node 2i
    //      right child is node 2i+1
    //      parent is node i/2

    void sink(int k, int N)  {
        int j;
        // while the moved element is greater than one of its children,
        // swap it with the smaller-valued child
        while (2*k <= N) {
            j = 2*k;
            if ( (j < N) && (list[j].val < list[j+1].val) ) j++;
            if ( !(list[k].val < list[j].val ) ) break;
            std::swap(list[k], list[j]);
            k = j;
        }
    }

    void swim(int k) {
        // while new element (k) is less than its parent (k/2), swap with parent
        while (k > 1 && ( list[k/2].val < list[k].val ) )  {
            std::swap(list[k], list[k/2]);
            k = k/2;
        }
    }

    int maxsize;
    T maxv;
    int heaplength;
};

template <class T>
class MinHeap : public MinHeapBase<T> {
public:
    inline void add(int i, T d) { this->insert( heapspace::Pair<T>(i, d) ); }
};

#endif // __MINHEAP_INCLUDE

#ifdef TEST_MINHEAP
// unit test
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(void) {

    // create large number of values to add to the heap
    int N = 1000;
    double values[N];
    srand48(12093123);
    for(int i=0; i<N; i++) values[i] = drand48();

    // create heap
    MinHeap<double> heap;
    int Nmin = 5;
    heap.alloc(Nmin);
    heap.clear();

    // add values to heap
    for(int i=0; i<N; i++) heap.add(i, values[i]);

    // extract results
    printf("heap.length() = %d\n", heap.length());
    printf("heap.maxval() = %e\n", heap.maxval());
    printf("Heap contents:\n");
    for(int i=0; i<heap.length(); i++) {
        printf("  %4d  % e\n", heap[i].id, heap[i].val);
    }

    // see if these are the same as the smallest values added
    while(heap.length() > 0) {
        // sort values
        std::sort( values, &(values[N]) );
        // get current maxval in heap
        double maxval = heap.maxval();
        // find this in the minimum heap.length() values
        int flag = 0;
        for(int i=0; i<heap.length(); i++) {
            if(values[i] == maxval) {
                flag = 1;
                // if found, "remove" it from the values
                values[i] = 1e308;
                break;
            }
        }
        // pop max value from heap and ensure it was the correct value
        assert( maxval == heap.removemax() );
        // ensure we found the value
        assert(flag==1);
    }
    printf("heap PASSED\n");

}
#endif
