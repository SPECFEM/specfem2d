#!/usr/bin/env python
from pylab import *
from os import *
from sys import *
import math

xmin,xmax,dx=0.,4000.,5.
ymin,ymax=0.,4000.

y0=2000.
# Amplitude de la sinusoide
a=50.
# Nombre de periodes
n=50
#
x = arange(xmin,xmax+dx,dx)
y = y0 + a*sin(2.*pi*n*(x-xmin)/(xmax-xmin))

x1,y1,x3,y3=xmin,y0,xmax,y0
x2,y2=(xmin+xmax)/2.,y0+a*3./2.

ma=(y2-y1)/(x2-x1)
mb=(y3-y2)/(x3-x2)

x1j=complex(x1,y1)
mx1j=abs(x1j)

xc=(ma*mb*(y1-y3)+mb*(x1+x2)-ma*(x2+x3))/2./(mb-ma)
yc=-(xc-(x1+x2)/2.)/ma+(y1+y2)/2.
xcj=complex(xc,yc)

R=abs(x1j-xcj)

ang1=math.atan2(y1-yc,x1-xc)
ang2=math.atan2(y3-yc,x3-xc)
angb=msort([ang1,ang2])
# Nombre de points pour la sinusoide
nsin=5000.
da = fabs((ang1-ang2)/nsin)
ang = arange(angb[0],angb[1]+da,da)

xi = xc + R*cos(ang)
yi = yc + R*sin(ang)


# Interface circulaire modulee

xim = xc + cos(ang)*(R + a*sin(2.*pi*n*(ang-ang[0])/(angb[1]-angb[0])))
yim = yc + sin(ang)*(R + a*sin(2.*pi*n*(ang-ang[0])/(angb[1]-angb[0])))

#
# Preparation des donnees en vue ecriture
#
DataS=zeros((len(x),2),Float)
DataC=zeros((len(xi),2),Float)
DataM=zeros((len(xim),2),Float)
DataS[:,0]=x
DataS[:,1]=y
DataC[:,0]=xi[::-1]
DataC[:,1]=yi[::-1]
DataM[:,0]=xim[::-1]
DataM[:,1]=yim[::-1]

DataM=asarray(DataM)
DataS=asarray(DataS)
DataC=asarray(DataC)

plot(x,y,'k',[xc],[yc],'ro')
plot(x,y,'k',xi,yi,'r')
plot(xim,yim,'b')
axis([xmin,xmax,ymin,ymax])
grid(True)
title('Interface sinusoidale')
#
# Ecriture du fichier interface
#
NbElemSpec=150
#
# Format d'ecriture en ascii des coordonnees
fmt='%f'
#
SEM_DIR='/home/cristini/SEM/SEM_2D_Dimitri'
#
chdir(SEM_DIR+'/DATA')
#
FileMod='InterfaceCercleMod'+str(int(a))+'_'+str(n)+'.dat'
print FileMod
#
f = file(FileMod,'w')
#
f.write("#\n# Number of interfaces\n#\n")
f.write("3\n")
f.write("#\n#\n#\n")
f.write("#\n# Interface 1 (Bottom of the mesh)\n#\n")
f.write('2\n%i %i\n%i %i\n'%(xmin,ymin,xmax,ymin))
f.write("#\n# Interface 2\n#\n")
f.write('%i\n'%len(xim))
for row in DataM:
   f.write(' '.join([fmt%val for val in row]) + '\n')
f.write("#\n# Interface 3\n#\n")
f.write('2\n%i %i\n%i %i\n'%(xmin,ymax,xmax,ymax))
f.write("#\n#Number of spectral elements\n#\n")
f.write("#\n#\n#\n")
f.write("%i\n"%NbElemSpec)
f.write("#\n#\n#\n")
f.write("%i\n"%NbElemSpec)
f.close()
#
show()

