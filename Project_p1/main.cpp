
#include <cmath>
#include <cassert>
//#include <cstdlib>
#include <random>
//#include <iostream>
//#include <vector>
//#include <sstream>
//#include <string>
#include <stdio.h>
#include <algorithm>

using namespace std;

class Gaussian {  // Stolen from the SPHFirst package
public:
  // Gaussian kernel normalized to unit integral over space
  double W(double rij, double h) {
    double q = rij/h;
    return exp(-q*q)/(pow(M_PI,1.5) * pow(h,3));
  }

  // and its h-derivative
  double dWdh(double rij, double h) {
    return W(rij, h)*(3*h*h-2*rij*rij)/(h*h*h);
  }

};




int main () {

  default_random_engine gen;
  uniform_real_distribution<double> d(-1.,1.);

  int nmax=10;
  double x[nmax];
  double y[nmax];

  for(int i=0;i<nmax;i++) {
    x[i]=d(gen);
    y[i]=d(gen);
  }

  int zz;
  for( int i=0; i < nmax; i++) {
    double myx=x[i];
    double myy=y[i];
    double dist[nmax-1];
    for( int j=0; j < i; j++) {
      // printf("%d\n",j);
      dist[j]=sqrt(pow(x[j]-myx,2.)+pow(y[j]-myy,2.));
    }
    for( int j=i+1; j< nmax; j++) {
      // printf("%d\n",j);
      dist[j-1]=sqrt(pow(x[j]-myx,2.)+pow(y[j]-myy,2.));

    }
    sort(dist,dist+nmax-1);


  }
  
  return 0;
}
