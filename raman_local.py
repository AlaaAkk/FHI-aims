"""
 !  Alaa Akkoush (Fritz Haber Institute)
 !  HISTORY
 !  February 2020 
 This script calculates the local polarizability of a cluster at different positions the user chooses, and 
 plots the tip enhanced raman images(TERS). You have to run this script after runing get_vibrtations.py 
 in mode (3).
 It works in three modes respectively:
 1. Mode (1). The folders for positive and negative displacement along normal modes are created.
   Example: python raman_local.py  H2O 1 -N 7 -f 0.003 -R  0 0 1 -n 5 5 -d 0.2 0.2 
   f: is the factor we multiply the normal modes with (-f #)
   N: The normal mode the user wants (-N #)
   R: The tip position in real space (-R Rx Ry Rz)
   n: The number of points in the grid (-n nx ny)
   d: the step size (-d dx dy)
   
 2. Mode(2). Runs aims local DFPT calculations in the created folder in Mode(1). 
   Example: python raman_local.py  -s srun -r $aims excutable path H2O 2 -N 7 -f 0.003 -R  0 0 1 -n 5 5 -d 0.2 0.2 -l 0 24
   l: The number of folders the users choose to run aims in in the order they are created, since there is a limit of 4000 
   srun in draco/cobra.
 3. Mode(3). Collect the polarizabilities and plots the TERS image.
    Example: python raman_local.py H2O 3 -N 7 -f 0.003 -R  0 0 1 -n 5 5 -d 0.2 0.2  -p
    p: for plotting.

"""
from scipy import constants, signal
from numpy import array, zeros, ones,ravel, float64, append, sqrt, arange,identity, newaxis, delete, linalg, sum, concatenate
from pylab import savetxt, transpose, eig, argsort, sort, sign, pi, dot, sum, linspace, argmin, r_, convolve
import copy, os, shutil
import numpy as np
import time
USAGE = """%prog [options] <name> <mode> 

<name> will be used as prefix for all output files.
<mode> selects preprocessing (and calculation if a path to FHI-aims binary is 
       is provides) and postprocessing - calcultion of vibrational modes from 
       forces obtained from DFT calculations.
       
"""
def split_line(lines):
  """Split input line"""
  line_array=array(lines.strip().split(' '))
  line_vals=line_array[line_array!='']
  return line_vals

def replace_submission(template_job, name, counter, filename):
    """
     Prepare submission script
     Only jobname and output file will be replaced
     Feel free to add more
    """
    template_job=template_job.replace('<jobname>',name+"_"+str(counter))
    template_job=template_job.replace('<outfile>',filename)
    job_out=open('job.sh','w')
    job_out.write(template_job)
    job_out.close()

def construct_grid(gridinfo):
    '''
    Returns a grid matrix of shape (ngrid, 3)
    '''
    orgx, orgy, orgz = gridinfo['org']
    nx, ny, nz = gridinfo['nstep']
    stepx, stepy, stepz = gridinfo['stepsize']
    x = np.linspace(orgx, orgx+(nx-1)*stepx, nx)
    y = np.linspace(orgy, orgy+(ny-1)*stepy, ny)
    z = np.linspace(orgz, orgz+(nz-1)*stepz, nz)
    gx, gy, gz = np.meshgrid(x, y, z, indexing='ij')
    gx = gx.flatten()
    gy = gy.flatten()
    gz = gz.flatten()
    return np.stack((gx, gy, gz)).transpose()
def check_if_output_exists(file_name, string_to_search):
    """ Check if any line in the file contains given string """
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if string_to_search in line:
                return True
    return False    

def main():
  import optparse
  from numpy import sum

  # Parse coxand line
  parser = optparse.OptionParser(usage=USAGE)
  parser.add_option("-i", "--info", action="store_true",
                      help="Set up/ Calculate vibrations & quit")
  parser.add_option("-p", "--plot", action="store_true",
                    help="Generate TERS image")

  parser.add_option("-s", "--suffix", action="store",
                    help="Call suffix for binary e.g. 'mpirun -np 12 '",
                    default='')
  parser.add_option("-t", "--submit", action="store",
                    help="""\
                    Path to submission script, string <jobname>
                    will be replaced by name + counter, string
                            <outfile> will be replaced by filename""")

  parser.add_option("-r", "--run", action="store",
                    help="path to FHI-aims binary",default='')
  parser.add_option("-N", "--number", action="store", type="int",
                    help="Mode number")
  parser.add_option("-f", "--fraction", action="store", type="float",
                    help="Finite difference fraction", default=0.01)
  parser.add_option("-R", action="store", type="float", nargs=3, dest="point")
  parser.add_option("-n", action="store", type="int", nargs=2, dest="step")
  parser.add_option("-d", action="store", type="float", nargs=2, dest="size")
  parser.add_option("-l", action="store", type="int", nargs=2, dest="bound")
  options, args = parser.parse_args()
  if options.info:
      print(__doc__)
      sys.exit(0)
  
  AIMS_CALL=options.suffix+' '+options.run
  name=args[0]
  mode=args[1]
  number=options.number
  fraction=options.fraction
  par=options.point
  n=options.step
  d=options.size
  l=options.bound
  p=options.plot
  if not options.number: 
      parser.error("Specify the vibrational mode you want by typing -N #")
  run_aims=False
  if options.run!='': run_aims=True
  
  submit_script = options.submit is not None
  
  # Asign all filenames
  inputgeomerty = 'geometry.in.'+name
  inputcontrol  = 'control.in.'+name
  
  irname   = name+'.data';

  

  f=open('control.in','r')                   # read control.in template
  template_control=f.read()
  f.close

  if submit_script:
    f=open(options.submit,'r')               # read submission script template
    template_job=f.read()
    f.close

  folder=''                                  # Duxy
 


  #constants
  c=constants.value('speed of light in vacuum')
  Ang=1.0e-10
  pi=constants.pi


# Construct scanning grid. Store coordinates in scanxyz
  gridinfo = {}
  gridinfo['org'] = [par[0], par[1],par[2]]
  gridinfo['stepsize'] = [d[0], d[1], 0]
  gridinfo['nstep'] = [n[0], n[1], 1]
  scanxyz = construct_grid(gridinfo)     
  #print('grid',scanxyz)
  #print(scanxyz.shape)
  # Checking for relaxed geometry 
  folder=name+'_relaxation' 
  if os.path.exists(folder+'/geometry.in.next_step'): 
    print("Relaxed geometry file found")
    geometry1=open(folder+'/geometry.in.next_step','r')
    geometry=open(folder+'/geometry.in.next_step','r')
    template_geometry=geometry.read()
    geometry.close
  else:
    print("Warning: directory name_relaxtion doesn't exit, non-relaxed geometry structure will be used")
    geometry1=open('geometry.in','r')
    geometry=open('geometry.in','r')
    template_geometry=geometry.read()
    geometry.close
  # Checking for Normal modes
  if os.path.exists('normalmodes.'+name+'.dat'): 
    print("Normal modes found")
    normalmodes=open('normalmodes.'+name+'.dat','r')
  else:
    print("Normal modes not found, run get_vibrations.py in mode 3")
   # Read input normal modes  
  num_line=1
  lines=normalmodes.readlines()
  norm_mode1=array
  for line in lines:
    if line.rfind('Atoms: ')!=-1:
         n_constrained= int(split_line(line)[1])  #3*N
    if num_line ==2*number:
        print('For mode '+ str(number)+', displacement along specified normal mode:\n'+str(line))
        norm_mode=line
       

    num_line=num_line+1
  if number > (num_line)/2:
      print('The mode you are requesting doesnt exist :)')
  # seperating displacement for each N
  q=zeros(shape=(n_constrained,3))
  x=zeros(shape=(n_constrained,3))
  coord=zeros(shape=(n_constrained,3))
  x_0=zeros(shape=(n_constrained,3))
  materials=''
  save=''
  norm_mode1=norm_mode
  norm_mode = norm_mode.split('[')[1].split(']')[0] # remove the bracket

 # Reading Mass weighted Normal modes
  j=0
  for i in range(n_constrained):

      q[i,:]= float64(split_line(norm_mode)[0+j:3+j])
      j=j+3
 # Normal mode in Cartesian Coordinates:

  temp=open('car_eig_vec.'+name+'.dat','r')
  lines_eig=temp.readlines()
  k=1
  
  for num in lines_eig:
     if k==number:
         qq=num

     k=k+1

  j=0
  for i in range(n_constrained):

      x[i,:]=float64(split_line(qq)[0+j:3+j])
      j=j+3


  #print('Normal mode in cartesian space: \n')+str(x)
 
  ##################################################
   # Read input geometry.in 
  n_line=0
  lines=geometry1.readlines()
  ii=0
  i=0
  lattice=[]
  slab=[]
  for line in lines:

      if line.rfind('constrain_relaxation .true.')!=-1 or line.rfind('Cu')!=-1:
           slab.append(line)
           i=i+1
      elif line.rfind('lattice_vector')!=-1:
           lattice.append(line)

      else:
           save=str(save+ line.split()[0]+'\n')
           x_0[ii,:]= float64(split_line(line)[1:4])
           materials = str(materials + line.split()[4]+'\n')
           ii=ii+1
  slab=array(slab)
  print(slab)
  lattice=array(lattice)
  n_atoms=ii
  print('Number of atoms:'+str(n_atoms))   
  print('Relaxed Geometry:\n'+str(x_0))

  #Displacement along normal mode:
  delta_x_pos=(x_0+(fraction*x))
  delta_x_neg=(x_0-(fraction*x))
  I=[]
  ############### 
  import numpy as np
  dq=fraction*q
  dx=fraction*x  
  dq_norm=np.linalg.norm(dq)
  dx_norm=np.linalg.norm(dx)
  
  print('||fraction*dx||=')+str(dx_norm)


  ###############################################################################################################################
  # creating geometry.in file along positive normal displacement
  i=0
  geo1 = []
  for lines in delta_x_pos:
      m1= save.splitlines()[i],lines[0], lines[1], lines[2],materials.splitlines()[i]
      geo1.append(m1)
      i=i+1
  geo1=array(geo1)
  # creating geometry.in file along negative normal displacement
  I=0
  geo2 = []
  for lines in delta_x_neg:
      m2= save.splitlines()[I],lines[0], lines[1], lines[2],materials.splitlines()[I]
      geo2.append(m2)
      I=I+1
  geo2=array(geo2)
  
# read M-1/2
  mass=open('mass_vec.'+name+'.dat','r')
  lines_mass=mass.readlines()
  k=0
  masses=[]

  for num in lines_mass:
      mass=num
      masses.append(mass.strip('\n').split(' '))            
      k=k+1
  masses=array(masses)
  masses=masses.astype(float)
 
  #print(masses)

  norm_mode1=array([float(i) for i in norm_mode.split()])
  #print(norm_mode1)

  value=dot((masses),transpose(norm_mode1))  
  #print('M^-1/2Li',linalg.norm(value))
  newline_ir='\n'
  if mode=='1': 
   # Set up / Read folders for displaced atoms
     for col in range(n[0]*n[1]):
        # for row in arange(3):
        filename='aims.out'
        folder1=name+'_positive_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])
        folder2=name+'_negative_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])


   ##### positive direction #######
        if not os.path.exists(folder1): os.mkdir(folder1)
        with open(folder1+'/geometry.in','w+') as f:
             for item in lattice:
                f.write("%s\n" % item)
             for item in slab:
                f.write("%s\n" % item)
                #print >> f, ' '.join(item)
             for item in geo1:
                print >> f, ' '.join(item)
        f.close    
        new_control=open(folder1+'/control.in','w')
        new_control.write(template_control+'DFPT local_polarizability '+ '{} {} {} '.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])+'1.3 1.3 5 1 1\n'+'DFPT_sc_accuracy_dm 0.001\n')
        new_control.close()
        os.chdir(folder1)                                   # Change directoy
        os.chdir('..')
 
 
  ###### Negative direction #####
        if not os.path.exists(folder2): os.mkdir(folder2)
        with open(folder2+'/geometry.in','w+') as f2:
             for item in lattice:
                f.write("%s\n" % item)
             for item in slab:
                f2.write("%s\n" % item)
               # print >> f2, ' '.join(item)
             for item in geo2:
          #sys.stdout.write(str(item))
                print >> f2, ' '.join(item)
          
  
        new_control=open(folder2+'/control.in','w')
        new_control.write(template_control+'DFPT local_polarizability '+ '{} {} {} '.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])+'1.3 1.3 5 1 1\n'+'DFPT_sc_accuracy_dm 0.001\n')
        new_control.close()
        os.chdir(folder2)                                   # Change directoy
        os.chdir('..')
#     print(n[0]*n[1]) 

  
  if mode=='2': 
   # Set up / Read folders for displaced atoms
     i=l[0]
     if  i>=l[0] and i<l[1]:
        for col in range(l[0],l[1]):
           # for row in arange(3):
           filename='aims.out'
           folder1=name+'_positive_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])
           folder2=name+'_negative_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])


   ##   ### positive direction #######
           if os.path.exists(folder1): 
              os.chdir(folder1)                                   # Change directoy
              if run_aims:
                   os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
              os.chdir('..')
           time.sleep(4)  
  ###   ### Negative direction #####
           if os.path.exists(folder2): 
              os.chdir(folder2)                                   # Change directoy
              if run_aims:
                   os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
              os.chdir('..')
 
           time.sleep(4)  
        i=i+1
  if mode=='3':
     for col in range(n[0]*n[1]):
        # for row in arange(3):
        filename='aims.out'
        folder1=name+'_positive_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])
        folder2=name+'_negative_disp_'+'{}_{}_{}'.format(scanxyz[col][0], scanxyz[col][1], scanxyz[col][2])
  # read polarizability for positive 
        os.chdir(folder1)                                   # Change directoy
        if os.path.exists(filename):
           data1=open(filename)
           if check_if_output_exists(filename, 'Have a nice day.'): 
              os.chdir('..')
           else:
                  if run_aims:
                    os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
                  os.chdir('..')
        else:
            if run_aims:
              os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
            os.chdir('..')
        data1=open(folder1+'/'+filename)
        lines= data1.readlines()
        for line in lines:
            if line.rfind('Polarizability')!=-1:
                alpha1 = float64(split_line(line)[4]) # alpha_zz 
                #print('pos',alpha1)
        os.chdir(folder2)                                   # Change directoy
        if os.path.exists(filename):
           data2=open(filename)
           if check_if_output_exists(filename, 'Have a nice day.'): 
              os.chdir('..')
           else:
                  if run_aims:
                    os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
                  os.chdir('..')
        else:
            if run_aims:
              os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
            os.chdir('..')
        data2=open(folder2+'/'+filename)
        lines= data2.readlines()
        for line in lines:
            if line.rfind('Polarizability')!=-1:
                alpha2 = float64(split_line(line)[4])
                #print('neg',alpha2)
        # Intensity 
        I=abs((alpha1-alpha2)/(2*(dq_norm)))**2 # Bohr^4/amu 
        # Saving Intensity
        newline_ir=newline_ir+'{} {} {}\n'.format(scanxyz[col][0], scanxyz[col][1],I)

     
        ir=open(irname,'w')
        ir.writelines('#x    y   intensity')
        ir.writelines(newline_ir)
     ir.close()
     if options.plot: 
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from matplotlib.colors import LogNorm
        import numpy as np
        import matplotlib.transforms as tr
        from matplotlib import rc
        from matplotlib import rcParams

        rc('text', usetex=True)
        SMALL_SIZE = 14
        MEDIUM_SIZE = 16
        BIGGER_SIZE = 20
        plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        #plt.rcParams['font.family'] = 'serif'
        #plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        params = {
                'axes.labelsize': 20,
                'font.size': 20,
                'xtick.labelsize': 20,
                'ytick.labelsize': 20,
                'mathtext.fontset': 'stix',
                #'font.family': 'sans-serif',
                #'font.sans-serif': 'Arial',
                'axes.linewidth': 2.0
            }
        rcParams.update(params)
        
        fig = plt.figure(figsize=(8.1, 6.5))
        
        
            
        gridx, gridy, gridz = np.loadtxt(name+'.data', unpack=True)
        gridx=gridx*0.52917721
        gridy=gridy*0.52917721
        
        N = int(len(gridz)**.5)
        #z=z*0.078415972
        Z = gridz.reshape(N, N)
        Z=Z.T
        fg=plt.imshow(Z, origin='lower',extent=(np.amin(gridx), np.amax(gridx), np.amin(gridy),  np.amax(gridy)),
                cmap=plt.cm.jet,aspect='equal', interpolation='bicubic')
        #rot = tr.Affine2D().rotate_deg(180)
        #fg.set_transform(fg.get_transform()+rot)
        #plt.xlim(-6,6)
        #plt.ylim(-2,9.5)
        #plt.colorbar()
        ####################################
        if os.path.exists('geometry.in'):
            print(" geometry file found")
            geometry=open('geometry.in','r')
            geometry.close
        
        n_line=0
        lines=geometry.readlines()
        ii=0
        for line in lines:
        
            if line.rfind('atom')!=-1:
               
       
               coord[ii,:]= float64(split_line(line)[1:4])
        
               ii=ii+1
        
       
        if coord is not None:
        
           x = coord[:, 0]
           y = coord[:, 1]
           z = coord[:, 2]
           a = np.argsort(z)
        
           plt.plot(x[a], y[a], 'wo', markersize=2, mew=2, color='white')
        ##################################
        plt.axis(aspect='image')
        plt.title('TERS image of benzene\nz=1\AA, $\sigma_{x}=\sigma_{y}=1.3$ \AA, $\sigma_{z}=5$ \AA' )
        plt.xlabel('Distance, \AA')
        plt.ylabel('Distance, \AA')
           
        cbar = plt.colorbar(fg, ticks=[Z.min(), Z.max()])
        cbar.set_label('Intensity', rotation=270,  labelpad=20, y=0.5)
        cbar.ax.set_yticklabels(['low', 'high'])  # vertically oriented colorbar
        plt.tight_layout()
        plt.savefig(name+'.png',transparent=True, dpi=400)
        plt.show()
        
if __name__ == "__main__":
    main()
