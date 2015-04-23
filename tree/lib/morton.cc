/*
Time-stamp: <morton.cc on Thursday, 12 March, 2015 at 08:47:25 MST (pinto)>

Morton Key code for D=1, 2, 3
 (to simplify, Morton keys for D<3 simply have zeros in the appropriate places of the 3D key)

  A MortonKey is a 128-bit (16 byte) number whose leftmost 32 bits are the particle index; the
    remaining 96 bits are the three integer coordinates interleaved into a Morton key.

 */


#ifndef __MORTON__
#define __MORTON__

#define uint128 __uint128_t
#define uint64 __uint64_t
#define uint32 __uint32_t
#define int64 __int64_t
#define int32 __int32_t

// I/O convenience functions
#define BITL(x, bit) (( (x)>>(bit) )&1UL)

void out128(uint128 r) {
    uint64 l,h;
    h = (r>>64);
    l = r;
    printf("%016llX%016llX\n",h,l);
}

void out64Binary(uint64 x) {
    for(int l=63; l>=0; l--) {
        printf("%llu", BITL(x,l) );
    }
}

void out32Binary(uint32 x) {
    for(int l=31; l>=0; l--) {
        printf("%lu", BITL(x,l) );
    }
}

void out96Binary(uint128 x) {
    uint64 l,h;
    h = (x>>64);
    l = x;
    uint32 ll;
    ll = h;
    out32Binary(ll);
    out64Binary(l);
}

void out128Binary(uint128 x) {
    uint64 l,h;
    h = (x>>64);
    l = x;
    out64Binary(h);
    out64Binary(l);
}
#undef BITL


template <int D>
struct Mask {
    static const SmallVec<int, D> mask[1<<D];
};

// direction masks: NB: x coordinate changes most rapidly
template <>
const SmallVec<int, 3> Mask<3>::mask[8] = { {-1,-1,-1},
                                            { 1,-1,-1},
                                            {-1, 1,-1},
                                            { 1, 1,-1},
                                            {-1,-1, 1},
                                            { 1,-1, 1},
                                            {-1, 1, 1},
                                            { 1, 1, 1} };

template <>
const SmallVec<int, 2> Mask<2>::mask[4] = { {-1,-1},
                                            { 1,-1},
                                            {-1, 1},
                                            { 1, 1} };

template <>
const SmallVec<int, 1> Mask<1>::mask[2] = { {-1},
                                            { 1} };

// base class; uses 128 bit key for all dimensions
class BaseMortonKey {
public:

    BaseMortonKey();
    ~BaseMortonKey();
    int key(uint128 m, int l);
    int pindex(uint128 m);

    void outKey(uint128 r);

    typedef uint128 mortonkey;

    void InitializeMortonConstants(void);
    uint128 dilate3_32(int t);
    uint32 double2int(double x);

    uint128 m1, m2, c1, c2, c3, c4, c5, c6, upper32mask0;

    // Magic numbers for double2int: note that the double 0x1p+0 cannot be converted this way
    // so that the range of numbers is, strictly,  [0, 1)
    static constexpr double MAGIC = 6755399441055744.0; // 2^52 + 2^51
    static constexpr double MAXIMUMINTEGER = 4294967294.0; // 2^32 - 2
};

bool AM(uint128 x, uint128 y) {
    uint128 xx = (x << 32);
    uint128 yy = (y << 32);
    return ( xx < yy );
}

struct AscendingMorton {
    bool operator() (uint128 x, uint128 y) {
        return AM(x,y);
    }
};

BaseMortonKey::BaseMortonKey() {
    // Check for portably defined types
    assert(sizeof(uint32)  ==  4);
    assert(sizeof(int32)   ==  4);
    assert(sizeof(uint64)  ==  8);
    assert(sizeof(int64)   ==  8);
    assert(sizeof(uint128) == 16);

    InitializeMortonConstants();
}

BaseMortonKey::~BaseMortonKey() {
}

// return the particle index from the MortonKey
inline int BaseMortonKey::pindex(uint128 m) {
    return (m>>96);
}


#define ONE ((uint128) 1)
void BaseMortonKey::InitializeMortonConstants(void) {

    m1 = (ONE<<64)+1;
    m2 = (ONE<<64) + (ONE<<32) + 1;

    uint128 x;

    x = 0x00000000ffffffffUL;
    c1 = (x<<96)+x;

    x = 0xffffUL;
    c2 = (x<<96) + (x<<48) + x;

    x = 0xffUL;
    c3 = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x;

    x = 0xf00fUL;
    c4 = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x;

    x = 0xc30c3UL;
    c5  = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x;

    x = 0x249249UL;
    c6 =  (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x;

    x = 0xffffffffUL;
    upper32mask0 = ~(x<<96);

}
#undef ONE

// spread the bits of t 3 apart: i.e. t=1011 becomes r=001 000 001 001
uint128 BaseMortonKey::dilate3_32(int t) {
    uint128 r = t;

    r = (r*m1) & c1;
    r = (r*m2) & c2;
    r = (r*0x100010001UL) & c3;
    r = (r*0x10101UL)     & c4;
    r = (r*0x00111UL)     & c5;
    r = (r*0x00015UL)     & c6;
    r = r & upper32mask0;
    return r;
}

// convert a double on [0, 1) to an unsigned 32 bit integer
uint32 BaseMortonKey::double2int(double d) {
    double _t = d * MAXIMUMINTEGER + MAGIC;
    return *( (uint32 *) &(_t) );
}

// extract byte b from the 96 bit key
inline int BaseMortonKey::key(uint128 m, int b) {
    int shr = 93-3*(b-1);
    return (m>>shr) & 7UL;
}

void BaseMortonKey::outKey(uint128 r) {
    // top 32 bits are the particle index stored as an unsigned int at the top
    unsigned int k = (r>>96);
    printf("index: %10u morton: ", k);
    // next 96 bits is the morton key, printed as octal digits
    for(int i=1; i<=32; i++) {
        printf("%d", key(r,i) );
    }
    printf("\n");
}


// ====================================================================================================

template <int D>
class MortonKey : public BaseMortonKey {
public:

    typedef SmallVec<int, D> ivec;
    typedef SmallVec<float, D> fvec;
    typedef SmallVec<double, D> dvec;

    void initMortonKey(dvec boxmin, dvec boxmax);
    template <class U>
    uint128 Morton(SmallVec<U, D> pos, int n);
    uint128 Morton(SmallVec<uint32, D> pos, int n);

    ivec mask[1<<D];
    dvec scale, boxmin, boxmax;
};

// boxmin and boxmax define the bounding box in user coordinates to be mapped to [0,1)^D
template <int D>
void MortonKey<D>::initMortonKey(dvec _boxmin, dvec _boxmax) {
    boxmin = _boxmin;
    boxmax = _boxmax;
    for(int d=0; d<D; d++) scale[d] = 1.0/(boxmax[d]-boxmin[d]);
}

// convert a position vector pos and particle index n into a 128-bit MortonKey
template <int D>
template <class U>
uint128 MortonKey<D>::Morton(SmallVec<U, D> pos, int n) {
    dvec pos_scaled = pos - boxmin;
    for(int d=0; d<D; d++) pos_scaled[d] *= scale[d];

    SmallVec<uint32, D> ii;
    for(int d=0; d<D; d++) ii[d] = double2int(pos_scaled[d]);
    return Morton(ii, n);
}

// convert a scaled integer position and index into a 128-bit MortonKey
template <int D>
uint128 MortonKey<D>::Morton(SmallVec<uint32, D> pos, int n) {
    uint128 result = ((uint128)n<<96);
    for(int d=0; d<D; d++) result |= (dilate3_32(pos[d])<<(d));
    return result;
}

#endif // __MORTON__
