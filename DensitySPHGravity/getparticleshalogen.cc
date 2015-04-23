/*
Time-stamp: <getparticleshalogen.cc on Thursday, 9 April, 2015 at 10:29:18 MST (pinto)>

 Read a Halogen for MUSE IC.ascii file with mass, positions, and velocities
 in to a Particleset<D>
 */


int fileLines(std::string name) {
    std::ifstream file(name);
    std::size_t lines_count =0;
    std::string line;
    while (std::getline(file , line)) ++lines_count;
    file.close();
    return lines_count;
}

template <int D>
void getParticlesHalogen(const char *fname, Particleset<D> &pset, SmallVec<double, D> &center, double &halfwidth) {

    int N = fileLines(fname);
    pset.alloc(N);

    FILE *fp = fopen(fname,"r");
    SmallVec<double, D> com(0), vcom(0);
    double mtot = 0;
    for(int i=0; i<pset.np; i++) {
        int q;
        double m;
        SmallVec<double, 3> pos, vel;
        int dummy = fscanf(fp,"%d %le %le %le %le %le %le %le",
                           &q, &m, &(pos[0]), &(pos[1]), &(pos[2]), &(vel[0]), &(vel[1]), &(vel[2]) );
        pset[i].id = q;
        pset[i].m = m;
        for(int d=0; d<D; d++) {
            pset[i].r[d] = pos[d];
            pset[i].v[d] = vel[d];
        }

        com += m*pset[i].r;
        mtot += m;
    }
    fclose(fp);

    // subtract off to put COM at centre
    com /= mtot;
    for(int i=0; i<pset.np; i++) {
        pset[i].r -= com;
    }

    // find halfwidth which defines bounding cube
    SmallVec<double, D> maxpos(-1e308), minpos(1e308);
    for(int p=0; p<pset.np; p++) {
        for(int d=0; d<3; d++) {
            minpos[d] = std::min(minpos[d], pset[p].r[d]);
            maxpos[d] = std::max(maxpos[d], pset[p].r[d]);
        }
    }
    //    printf("x on [% e, % e]\n", minpos[0], maxpos[0]);
    //    printf("y on [% e, % e]\n", minpos[1], maxpos[1]);
    //    printf("z on [% e, % e]\n", minpos[2], maxpos[2]);

    halfwidth = 0;
    for(int d=0; d<3; d++) {
        halfwidth = std::max(halfwidth, fabs(minpos[d]));
        halfwidth = std::max(halfwidth, fabs(maxpos[d]));
    }
    center.zero();
    fpprint(std::cerr,"got %d particles\n", pset.np);
    fpprint(std::cerr,"center: % e   halfwidth: %e\n", center, halfwidth);

}
