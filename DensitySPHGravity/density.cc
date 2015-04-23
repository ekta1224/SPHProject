
template <int D>
void getDensity1st(Particleset<D> &p, Octree<D> &OT, int Nneighbours, double eta) {

    SPHDensity<KERNEL<D>, D> doit(eta);

    int nprocs = omp_get_num_procs();
    int nsource = 3*Nneighbours;
    dvec source_r[nprocs][nsource];
    double source_m[nprocs][nsource];
    int ids[nprocs][nsource];

#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<p.np; i++) {
        int g = omp_get_thread_num();
        double r_max;
        OT.knn(p[i].r, nsource, r_max, ids[g]);
        p[i].h = r_max;
        assert(p[i].h > 0);
        for(int k=0; k<nsource; k++) {
            assert( ids[g][k] < p.np );
            source_r[g][k] = p[ids[g][k]].r;
            source_m[g][k] = p[ids[g][k]].m;
        }
        doit.getrho(p[i].r, p[i].m, source_r[g], source_m[g], nsource, p[i].h, p[i].rho, p[i].omega, p[i].nn);
   }
}

template <int D>
void getDensity(Particleset<D> &p, Octree<D> &OT, int Nneighbours, double eta) {

    SPHDensity<KERNEL<D>, D> doit(eta);

    int nprocs = omp_get_num_procs();
    int nsource = p.np;
    dvec source_r[nprocs][nsource];
    double source_m[nprocs][nsource];
    int ids[nprocs][nsource];

#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<p.np; i++) {
        int g = omp_get_thread_num();
        int nids;
        OT.ballSearch(p[i].r, 1.2*p[i].h, ids[g], nids);
        for(int k=0; k<nids; k++) {
            assert( ids[g][k] < p.np );
            source_r[g][k] = p[ids[g][k]].r;
            source_m[g][k] = p[ids[g][k]].m;
        }
        doit.getrho(p[i].r, p[i].m, source_r[g], source_m[g], nids, p[i].h, p[i].rho, p[i].omega, p[i].nn);
   }
}
