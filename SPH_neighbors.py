from pylab import *
import math
import numpy
from numpy import random
import matplotlib

nump=100; nclose=5
xs=uniform(-1.,1.,nump)
ys=uniform(-1.,1.,nump)
length=[]
coord=array(zip(xs,ys,range(nump)))

for i in range(len(coord)):
    me = coord[i]
    them=delete(coord[:],i,0)
    dist=sqrt(add(square(subtract(them[:,0],me[0])),square(subtract(them[:,1],me[1]))))
    near= them[argsort(dist)][:nclose]
    neard=dist[argsort(dist)][:nclose]
    length += [neard[-1]]

#print length

figure(figsize=(10,10))
errorbar(xs,ys,xerr=length,yerr=length,marker='o',color='k',ls='None')
show()
