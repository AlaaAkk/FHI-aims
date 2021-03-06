"""
 !  Alaa Akkoush (Fritz Haber Institute)
 !  HISTORY
 !  February 2020 
 This script calculates born effective charges using finite difference of polarization with displacements.
 1.The user chooses which atoms need to be dispaced 
 2.The user can choose the direction of displacement x, y or z (-c 1/2/3)
 3.The user can choose the finite difference poles (defauls between 0 and 0.0025) -d d1 d2 
 4.The user can choose the grid size (-n nx nx nz)

 example: 
 Displacing Mg atoms along z direction with a grid 1x1x2, finite difference between 0.005 and 0.01 AA:
 python BEC-aims.py -r path-to-aims-excutable  Mg -n 1 1 2 -c 3 -d 0.005 0.01
 
 Notes: 1.The provided structure should have orthogonal unit cell.
        2.You can choose not to run aims but displace folders and run seperately.

"""
from scipy import constants, signal
from numpy import (
    array,
    zeros,
    ones,
    ravel,
    float64,
    append,
    sqrt,
    arange,
    identity,
    newaxis,
    delete,
    linalg,
    sum,
)
from pylab import (
    savetxt,
    transpose,
    eig,
    argsort,
    sort,
    sign,
    pi,
    dot,
    sum,
    linspace,
    argmin,
    r_,
    convolve,
)
import copy, os, shutil

#from ase.io import read, write
#from ase.build import bulk
import numpy as np
from sklearn import preprocessing

USAGE = """%prog [options] <name> 

<name> will be used as prefix for all output files.
<mode> selects preprocessing (and calculation if a path to FHI-aims binary is 
       is provides) and postprocessing - calcultion of vibrational modes from 
       forces obtained from DFT calculations.
       
"""


def split_line(lines):
    """Split input line"""
    line_array = array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


def replace_submission(template_job, name, counter, filename):
    """
     Prepare submission script
     Only jobname and output file will be replaced
     Feel free to add more
    """
    template_job = template_job.replace("<jobname>", name + "_" + str(counter))
    template_job = template_job.replace("<outfile>", filename)
    job_out = open("job.sh", "w")
    job_out.write(template_job)
    job_out.close()


def main():
    import optparse
    from numpy import sum

    # Parse coxand line
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option(
        "-i", "--info", action="store_true", help="Set up/ Calculate vibrations & quit"
    )

    parser.add_option(
        "-s",
        "--suffix",
        action="store",
        help="Call suffix for binary e.g. 'mpirun -np 12 '",
        default="",
    )
    parser.add_option(
        "-t",
        "--submit",
        action="store",
        help="""\
                    Path to submission script, string <jobname>
                    will be replaced by name + counter, string
                            <outfile> will be replaced by filename""",
    )
    parser.add_option("-n", action="store", type="int", nargs=3, dest="grid")
    parser.add_option(
        "-r", "--run", action="store", help="path to FHI-aims binary", default=""
    )
    parser.add_option(
        "-d",
        "--delta",
        action="store",
        type="float",
        nargs=2,
        dest="delta",
        help="Displacement",
        default=[0.00, 0.0025],
    )

    parser.add_option(
        "-c", "--direction", action="store", type="int", help="direction", default=3
    )
    options, args = parser.parse_args()
    if options.info:
        print(__doc__)
        sys.exit(0)

    AIMS_CALL = options.suffix + " " + options.run
    name = args[0]
    deltas = options.delta
    c = options.direction
    n = options.grid
    c_zero = 1.0 / abs(deltas[1] - deltas[0])
    coeff = array([-1, 1])
    if c != 1 and c != 2 and c != 3:
        parser.error("directions can be either 1,2 or 3")
    run_aims = False
    if options.run != "":
        run_aims = True

    submit_script = options.submit is not None

    # Asign all filenames
    inputgeomerty = "geometry.in." + name
    inputcontrol = "control.in." + name

    f = open("control.in", "r")  # read control.in template
    template_control = f.read()
    f.close

    if submit_script:
        f = open(options.submit, "r")  # read submission script template
        template_job = f.read()
        f.close

    folder = ""  # Duxy

    filename = "aims.out"

    # constants
    C = 1.6021766e-19  # in coulomb
    p = array([])
    p1 = array([])
    p2 = array([])
    p3 = array([])
    rec_lat=zeros(shape=(3,3))
    if os.path.exists("geometry.in"):
        geometry = open("geometry.in", "r")
        geometry.close
    n_line = 0
    lines = geometry.readlines()
    m=0
    for l in lines:
        if l.rfind('lattice_vector') != -1:
            rec_lat[m] = l.split()[1:4]
            m=m+1
    #find the reciprocal lattice vectors
    rec_lat= np.linalg.pinv(rec_lat).transpose()

    # Normalize the lattice vectors
    X_normalized = preprocessing.normalize(rec_lat, norm="l2")
    # SOlving 3 eq 3 unknowns to find linear combination
    x = np.linalg.solve(X_normalized, np.array([1, 0, 0]))
    y = np.linalg.solve(X_normalized, np.array([0, 1, 0]))
    z = np.linalg.solve(X_normalized, np.array([0, 0, 1]))
    R = np.array([x, y, z])

    i = 0
    for delta in deltas:
        geo = []

        folder = name + "_disp_" + str(delta) + "_{}_{}_{}".format(n[0], n[1], n[2])
        if c == 1:
            ii = 0
            for line in lines:
                if line.rfind(name) != -1:

                    d = float(line.split()[1]) + delta
                    new = (
                        line.split()[0],
                        str(d),
                        line.split()[2],
                        line.split()[3],
                        line.split()[4],
                    )
                    # new = new.split('(')[1].split(')')[0]
                    geo.append(new)
                    ii = ii + 1

                else:
                    new2 = line.split()
                    geo.append(new2)
        elif c == 2:
            ii = 0
            for line in lines:
                if line.rfind(name) != -1:

                    d = float(line.split()[2]) + delta
                    new = (
                        line.split()[0],
                        line.split()[1],
                        str(d),
                        line.split()[3],
                        line.split()[4],
                    )
                    geo.append(new)
                    ii = ii + 1

                else:
                    new2 = line.split()
                    geo.append(new2)
        elif c == 3:
            ii = 0
            for line in lines:
                if line.rfind(name) != -1:

                    d = float(line.split()[3]) + delta
                    new = (
                        line.split()[0],
                        line.split()[1],
                        line.split()[2],
                        str(d),
                        line.split()[4],
                    )
                    geo.append(new)
                    ii = ii + 1

                else:
                    new2 = line.split()
                    geo.append(new2)

        if not os.path.exists(folder):
            os.mkdir(folder)
        with open(folder + "/geometry.in", "w+") as f:
            for item in geo:

                print >> f, " ".join(item)
        f.close
        new_control = open(folder + "/control.in", "w")
        new_control.write(
            template_control
            + "KS_method serial \n"
            + "output polarization    "
            + str(1)
            + " {} {} {}\n".format(n[0], n[1], n[2])
            + "output polarization    "
            + str(2)
            + " {} {} {}\n".format(n[0], n[1], n[2])
            + "output polarization    "
            + str(3)
            + " {} {} {}\n".format(n[0], n[1], n[2])
        )

        new_control.close()
        os.chdir(folder)
        # Change directoy

        if run_aims:
            os.system(
                AIMS_CALL + " > " + filename
            )  # Run aims and pipe the output into a file named 'filename'
        if submit_script:
            replace_submission(template_job, name, 0, filename)
        # os.system('sbatch -N 4 job.sh') # Mind the environment variables <======= Modify according to you system

        os.chdir("..")
        # checking existence of aims.out
        if os.path.exists(folder + "/" + filename):
            data1 = open(folder + "/" + filename)
            liness = data1.readlines()
            for line in liness:
                if (
                    line.rfind(
                        "Detailed listing before/after branch mapping in terms of the polarization quantum P0="
                    )
                    != -1
                ):
                    q0 = float64(split_line(line)[-2:-1])  # polarization quantum
                if (
                    line.rfind(
                        "- Directive    1 in direction of rec. latt. vec.  "
                        + '1'
                        + " yields the full polarization      :"
                    )
                    != -1
                ):
                    p_1 = float64(split_line(line)[-2:-1])  #
                    if round(abs(p_1) / abs(q0)) != 0:
                       p_1 = p_1 - q0
                    p1 = append(p1, p_1)
                if (
                    line.rfind(
                        "- Directive    2 in direction of rec. latt. vec.  "
                        + '2'
                        + " yields the full polarization      :"
                    )
                    != -1
                ):
                    p_2 = float64(split_line(line)[-2:-1])  #
                    if round(abs(p_2) / abs(q0)) != 0:
                       p_2 = p_2 - q0
                    p2 = append(p2, p_2)
                if (
                    line.rfind(
                        "- Directive    3 in direction of rec. latt. vec.  "
                        + '3'
                        + " yields the full polarization      :"
                    )
                    != -1
                ):
                    p_3 = float64(split_line(line)[-2:-1])  #
                    if round(abs(p_3) / abs(q0)) != 0:
                       p_3 = p_3 - q0
                    p3 = append(p3, p_3)
                if line.rfind('| Unit cell volume ')!=-1:
                    volume=float64(split_line(line)[-2])


    P_1 = np.array([p1[0], p2[0], p3[0]])
    P_2 = np.array([p1[1], p2[1], p3[1]])
    delta_1 = np.dot(P_1, R)
    delta_2 = np.dot(P_2, R)
    p = [delta_1[c - 1], delta_2[c - 1]]
    p = array(p)
    print("polarization for disp " + str(deltas[0]) + " AA is : " + str(p[0]))
    print("polarization for disp " + str(deltas[1]) + " AA is : " + str(p[1]))

    # Change unit to e:
    born_factor = (volume * 1e-20) / (1 * C)
    I = (p[1] - p[0]) / abs(deltas[1] - deltas[0])  # Finite difference of polarization
    print(
        "Born effective charge for " + str(name) + " is:  " + str(I * born_factor / ii)
    )


if __name__ == "__main__":
    main()
