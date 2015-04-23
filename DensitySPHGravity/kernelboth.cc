
    int srcids[p.np];

    for(int i=0; i<p.np; i++) {

        int nsrcids = 0;
        tree.ballSearch(p[i].r, p[i].h, srcids, nsrcids);
        assert( nsrcids > 0 );

        for(int k=0; k<nsrcids; k++) {
            int j = srcids[k];
            if(i!=j) {
                dvec dr = p[i].r - p[j].r;
                double r = dr.norm();
                if( r < p[i].h ) {
                    dvec a = 0.5 * K.forcekernel(r, p[i].h) * dr/r;
                    double b = 0.5 * K.phikernel(r, p[i].h);
                    dvec zippo = dr/(r*r*r);

                    p[i].acc += p[j].m * ( a - 0.5 * zippo );
                    p[i].pot += p[j].m * ( b + 0.5/r );
                    p[j].acc -= p[i].m * ( a - 0.5 * zippo );
                    p[j].pot += p[i].m * ( b + 0.5/r );
                }
            }
        }
    }
