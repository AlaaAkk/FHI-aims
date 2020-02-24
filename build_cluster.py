"""
This code builds o cluster benzene on copper from bulk of six layer
"""

from ase.build import molecule
from ase.build import fcc111, add_adsorbate
from ase.io import read, write
from ase import Atoms
from numpy import array, zeros, ones,ravel, float64, append, sqrt, arange,identity, newaxis, delete,linalg, sum
from numpy import genfromtxt
import numpy as np
from pylab import savetxt, transpose
import copy, os, shutil


################# Slab with metal and molecule #####################
slab = fcc111('Cu', size=(3,3,6), a=3.49331674)
slab.center(vacuum=50.0, axis=2)
atoms = molecule('C6H6')
atoms.write('molecule.in',format='aims')
structure=add_adsorbate(slab, atoms, 2.69,position=(3.7051999999999996,3.5344500000000005))
slab.write('system.in',format='aims')#geometry is written

g=atoms.get_center_of_mass(scaled=False)  #center of mass of the molecule
print('center of mass is found : '+ str(g))

#################### Slab ith metal only ##################################
slab2 = fcc111('Cu', size=(3,3,6), a=3.49331674)
slab2.center(vacuum=50.0, axis=2)
slab2.write('slab.in',format='aims')#geometry is written
geomerty = read('slab.in', 0, 'aims')




# ASE image 
write('system.png', slab * (3, 3, 1), rotation='10z,-80x')

#Reading the molecule from the slab
geomerty = read('system.in', 0, 'aims')
molecule=open('system.in','r')
lines=molecule.readlines()
geo_mol=[]
i=0
for line in lines:
    if line.rfind('atom')!=-1:
       if i> 43: 
      
          geo_mol.append(line)
    i=i+1      
geo_mol=array(geo_mol) #adsorbed molecule  
# xx is the choosen distance from the center of mass:
# d is the absolute distance between the atoms in the bulk 
# and the center of mass of the molecule:
xx=13

materials=slab2.get_chemical_symbols()
materials=materials[0]
i=0
geo_metal=[]
rr=slab2.get_positions()

for j in rr:
     
    dd=(j-g)
    dr=np.linalg.norm(dd)
    print(dr) 
    if i <9:
   
         #print(i)
         #print(rr[i,:])
         if dr < 0.75*xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
            
    elif i >= 9 and i<18:
         if dr < 0.6*xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
    elif i >= 18 and i<27:
         if dr < 0.7*xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
    elif i >= 27 and i<36:
         if dr < 0.8*xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
    elif i>= 36 and i<45:
         if dr < 0.9*xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
    elif i>= 45 and i < 54:
         if dr < xx:
            new='atom',j[0] , j[1], j[2], materials
            geo_metal.append(new)
    i=i+1 
geo_metal=array(geo_metal)
print i

############## Creating geomety for cluster ##########################
with open('geometry.in', 'w') as f:
    for item in geo_metal:
        #f.write("%s\n" % item)
        print >> f, ' '.join(item) 
    for item in geo_mol:
        f.write("%s\n" % item)





