// equivalent to doing a tree all on all
#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<p.np; i++) {
        for(int j=0; j<p.np; j++) {
            if(i!=j) {
                dvec dr = p[i].r - p[j].r;
                double r = dr.norm();
                dvec Fij =  p[j].m * dr/(r*r*r);
                double Potij = - p[j].m / r;
                p[i].acc += Fij;
                p[i].pot += Potij;
            }
        }
    }
