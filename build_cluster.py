from ase.build import molecule
from ase.build import fcc111, add_adsorbate
from ase.io import read, write
from ase.build import make_supercell
from ase import Atoms
from numpy import array, zeros, ones,ravel, float64, append, sqrt, arange,identity, newaxis, delete,linalg, sum
from numpy import genfromtxt
import numpy as np
from pylab import savetxt, transpose
import copy, os, shutil


################# Slab with metal and molecule #####################
slab = fcc111('Cu', size=(3,3,6), a=3.49331674)
slab.center(vacuum=50.0, axis=2)
slab = make_supercell(slab, np.diag([3,3,1]))

atoms = molecule('C6H6')
atoms.write('molecule.in',format='aims')


gg=slab.get_center_of_mass(scaled=False)  #center of mass of the slab
print(gg)




structure=add_adsorbate(slab, atoms, 2.69,position=(gg[0],gg[1]+(1.3952479999999987))) #molecule on center
slab.write('system.in',format='aims')#geometry is written



g=atoms.get_center_of_mass(scaled=False)  #center of mass of the molecule
print('center of mass is found : '+ str(g))


#################### Slab with metal only ##################################
slab2 = fcc111('Cu', size=(3,3,6), a=3.49331674)

#slab2=slab2 * (3, 3, 1)
slab2.center(vacuum=50.0, axis=2)
slab2.write('slab.in',format='aims')#geometry is written
#atom = read("slab.in", format='aims')
slab2 = make_supercell(slab2, np.diag([3,3,1]))
write("slab.in", slab2, format="aims")
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
       if line.rfind('Cu')==-1:
           if i<12: 
              geo_mol.append(line)
           i=i+1      
geo_mol=array(geo_mol) #adsorbed molecule  geometry

# xx is the choosen distance from the center of mass:
xx= 10


materials=slab2.get_chemical_symbols()
materials=materials[0]
i=0
geo_metal=[]
rr=slab2.get_positions()  #position of metal in the slab
#print(slab2.get_tags())
for j in rr:
    dd=(j-g)
    dr=np.linalg.norm(dd)
#    print('distance from COM for atom '+str(i)+' is: '+str(dr))
    if dr < xx:
       new='atom',j[0] , j[1], j[2], materials
       geo_metal.append(new)
    i=i+1 
geo_metal=array(geo_metal)


############## Creating geomety for cluster ##########################
with open('geometry.in', 'w') as f: #saving it in geometry.in file
    for item in geo_metal:
        #f.write("%s\n" % item)
        print >> f, ' '.join(item) 
    for item in geo_mol:
        f.write("%s\n" % item)
###################Viewing the cluster using jmol#####################
os.system('jmol geometry.in')



