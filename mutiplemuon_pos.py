from math import *
import numpy as np
import os
import subprocess


"""
this method converts a theta, phi pair to its (x, y, z) on a unit sphere
"""
def changecord(theta, phi):
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)

    return np.array([x, y, z])


"""
This method converts a position (x, y, z) on earth to depth, latitude and longitude
depth: distance in km, negative if underground, positive if above
latitude: south pole as -90, north pole as 90
longitude: use x axis as cut of eastern and western, for western longitude is 0 to 180, for eastern is 0 to -180

the comment below can be used to test it.
"""
def convertcord(x, y, z):
    R = sqrt(pow(x, 2) + pow(y,2) + pow(z,2))
    theta = acos(z/R)
    try:
        phi = atan(y/x)
    except ZeroDivisionError:
        phi = pi/2
    if(x < 0 and y > 0):
        phi = pi + phi
    if(x <= 0 and y < 0):
        phi = phi - pi

    depth = R/1000 - 6371

    latitude = degrees(theta) - 90

    longitude = -1 * degrees(phi)

    return depth, latitude, longitude
"""
print(convertcord(6371000, 0, 0)) # should return 0, 0, 0
print(convertcord(6370000/2 * sqrt(2), 6370000/2 * sqrt(2), 0)) # should return -1, 0, -45
print(convertcord(0, 0, 6371000 - 2000)) # position of icecube, since theta = 0, phi does not matter (?)
exit()
"""




"""
This method calculate the magnetic field based on the position (x, y, z) on earth by calling wmm_point.exe
The return value of wmm_point.exe is in B_r, B_theta, B_phi form and will be converted to (Bx, By, Bz) in earth frame
and return

https://www.ngdc.noaa.gov/geomag/WMM/calculators.shtml

"""
def getB(x, y, z):

    """calling wmm_point.exe and extracting its output"""
    depth1, latitude1, longitude1 = convertcord(x, y, z)

    main = "./wmm_point.exe"

    latitude = latitude1
    longitude = longitude1
    height = depth1
    date = 2019.5

    p = subprocess.Popen(main, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)

    if(height < -10):
        rs, out = p.communicate(os.linesep.join(["c", str(latitude), str(longitude), str(height), "c", str(date), "c", "c"]))
    else:
        rs, out = p.communicate(os.linesep.join(["c", str(latitude), str(longitude), str(height), str(date), "c", "c"]))

    pre = str(rs).split('Secular Change')[1]

    array = pre.split('+/-')

    array = [x for x in array if x != ' ']

    # todo: check the conversion of coordinate
    # the return values of this method are north, east, vertical direction with NORTH pole as positive z
    # north: + N and -S; east: +E and -W; vertical: +D and -U
    # our coordinate use SOUTH pole as positive z axis, and use r, theta, phi
    # there are some conversions performed
    B_theta = float(array[2].split('X\t=\t')[1])        # North Component, so it is theta direction
    B_phi = -1 * float(array[3].split('Y\t=\t')[1])     # -1 * East Component, so its in phi direction
    B_r = -1 * float(array[4].split('Z\t=\t')[1])       # -1 * Vertical Component, so it is in r direction

    # now B_phi B_theta and B_r are B field in SOUTH pole as positive z axis cord
    # next step: convert them to (Bx, By, Bz) as well
    # get theta and phi for imput position (x, y, z)
    R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
    try:
        theta = acos(z / R)
    except ZeroDivisionError:
        theta = 0

    try:
        phi = atan(y / x)
    except ZeroDivisionError:
        phi = pi/2
    if (x < 0 and y > 0):
        phi = pi + phi
    if (x <= 0 and y < 0):
        phi = phi + pi
    if (x > 0 and y < 0):
        phi = phi + 2 * pi

    # projecting the (Br, Bthata, Bphi) on to (Bx, By, Bz)
    thetap = theta + pi / 2
    phip = phi + pi / 2

    B_x = B_theta * sin(thetap) * cos(phi) + B_phi * cos(phip) + B_r * sin(theta) * cos(phi)

    B_y = B_theta * sin(thetap) * sin(phi) + B_phi * sin(phip) + B_r * sin(theta) * sin(phi)

    B_z = B_r * cos(theta) + B_theta * cos(thetap)

    # multiplied by e-9 since it is in nT
    return np.array([B_x, B_y, B_z]) * pow(10, -9)

# tests
# print(getB(6371000, 0, 0))
# exit()


"""==================================== muon class ========================================"""


class muon(object):

    """
    Constructor
    Note that input v_dir is the direction we received from ice cube detector
    """
    def __init__(self, v_d, initial_E):
        self.initialvdir = v_d/np.linalg.norm(v_d) * -1         # initial direction of velocity, unit vec
        self.v_dir = v_d * -1
        self.v_dir = self.v_dir / np.linalg.norm(self.v_dir)    # current direction of velocity, unit vec
        self.initialE = initial_E                               # initial energy
        self.E = initial_E                                      # current energy, in eV

    # update reference, used to determine the time interval we use
    # it means how many pieces we want a circle to be cut into
    updateref = 5 * 10e-8

    # final direction of velocity
    finalvdir = np.array([0,0,0])

    # some constants
    const_c = 299792458  # m/s
    const_mass = 105.658 * np.power(10, 6)  # ev/c^2
    const_e = 1.602 * pow(10, -19)

    # current radius of muon path, with respect to the B field vector
    Radius = 0

    # magnitude of velocity
    v_total = 0
    # current position, initially set to be 0.0001 on x and y to avoid dividing by 0 issue when computing phi and theta
    pos = np.array([0.0001, 0.0001, 6371000 - 2000])

    # the position in previous step, used to do position correction in the final result (coming out of earth issue)
    prepos = np.array([])
    preenergy = 0       # energy in previous step
    totalDis = 0        # total distance traveled
    predis = 0          # total distance traveled in previous step
    totalAng = 0        # total deflection angle. change in deflection angle is tiny the final correction is ignored

    """
    convert energy in J to energy in eV
    """
    def convertev(self, E) :
        return E * self.const_e


    """
    get energy required for muon traveling distance dis
    using an initial energy E_ini in eV and dis in meter, according to the integration of 
    dE/dx = 0.259GeV/wme + 0.363/wme * E(GeV)
    dE/dx is eV per cm
    
    so E is (A e^bx  - a)/b, A is some constant
    plug x = 0, E = initial E, so that A = E_ini * b + a
    """
    def getnewE(self, E, dis):


        b = 0.000363

        a = 0.259

        density = 0.9167
        # todo: check this: the unit for a and b are /mwe, which is meter water equivalence, as the result
        # todo: density is used here to convert into ice.
        # todo: however on the original paper a and b are plotted for ice, so I don't know if density is needed
        # todo: paper is at  FIG.21 from arXiv:hep-ph/0407075v3

        # convert E in GeV to use this eq
        A = (E / pow(10, 9) * b * density) + a * density
        # convert E_final back in eV
        E = (A * np.power(e, b * dis * density) - a * density) / (b * density) * pow(10, 9)

        return E

    """
    calculate the velocity using energy
    v = sqrt[ c^2 - (mc^2 * c / E)^2 ]
    """
    def calV(self, E) :
        temp = self.convertev(self.const_mass) * self.const_c
        E = self.convertev(E)
        temp2 = temp/E
        temp3 = np.power(temp2, 2)
        temp4 = np.power(self.const_c, 2)
        return np.power(temp4 - temp3, 0.5)

    """
    calculate radius based on the energy E, v_dir, and B vectors
    cosT is used to project v to vertical component with respect to B vector
    so |v_total| * cosT = |v_vertical| and |v_total| * sinT = |v_parallel|
    
    according to R = sqrt[E^2 - (mc^2)^2] * cos^2(T)  /  [q * c * |v x B|]
    """
    def calR(self, v_dir, E, B, cosT):
        temp1 = np.power(self.convertev(self.const_mass), 2)
        temp2 = sqrt(np.power(self.convertev(E), 2) - temp1) * pow(cosT, 2)
        temp3 = self.const_e * self.const_c * np.linalg.norm(np.cross(v_dir, B))
        if temp3 == 0:
            return 0
        return temp2/temp3



    """-------------------------------- main part ------------------------------------"""
    def updatepos(self, t, p):

        """------------------------------- first order approximation ------------------------------------"""
        # get B field at current position, here -1 is multiplied since it is a backward simulation
        # where in reality muon is traveled forward
        const_B = getB(self.pos[0], self.pos[1], self.pos[2]) * -1

        # get current velocity
        self.v_total = self.calV(self.E)

        # doing 2 cross product, the result temp2 is the vertical component of v_dir to B field
        temp1 = np.cross(const_B, self.v_dir)
        temp2 = np.cross(temp1, const_B)

        # if v_dir is exactly parallel
        if np.linalg.norm(temp2) == 0:
            v_perpendiculardir = np.array([0,0,0])
            cosT = 0
        # if not, normalize vertical dir to unit vector (called temp3), get cos theta where theta is angle between v_dir
        # and v_perpendicular, then change the length of vertical dir to the actual ratio
        else:
            temp3 = temp2/np.linalg.norm(temp2)
            cosT = np.dot(temp3, self.v_dir)
            v_perpendiculardir = cosT * temp3

        v_perpendicular = self.v_total * cosT  # scalar representing the vertical speed component

        v_paralleldir = self.v_dir - v_perpendiculardir  # vector representing direction parallel with B

        # calculate current radius
        R = self.calR(self.v_dir, self.E, const_B, cosT)
        if R == 0:
            v_omg = 0
            omg = 0
        # get angular velocity v_omg and the rotation angle during time period t
        else:
            v_omg = v_perpendicular/R
            omg = v_omg * t

        # normalize B to be the unit vector as the rotation axis
        rotaxis = const_B/np.linalg.norm(const_B)

        # construct the rotational matrix
        M = np.array([[cos(omg) + (1 - cos(omg)) * np.power(rotaxis[0], 2),
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[1] - sin(omg) * rotaxis[2],
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[2] + sin(omg) * rotaxis[1]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[1] + sin(omg) * rotaxis[2],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[1], 2),
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] - sin(omg) * rotaxis[0]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[2] - sin(omg) * rotaxis[1],
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] + sin(omg) * rotaxis[0],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[2], 2)]])

        # perform rotation on perpendicular direction, and record the old v_perp
        old_perp = v_perpendiculardir
        v_perpendiculardir = np.dot(M, v_perpendiculardir)

        v_dir_new = v_perpendiculardir + v_paralleldir  # new velocity direction

        Dis = self.v_total * t  # total distance traveled in this time period

        E_new = self.getnewE(self.initialE, self.totalDis + Dis)    # new energy after this time period

        v_total_new = self.calV(E_new)  # new speed after this time period

        pos_new = self.pos + v_dir_new * v_total_new * t    # new position after this time period

        const_B_new = getB(pos_new[0], pos_new[1], pos_new[2]) * -1     # new B field after this time period

        """------------------------------------- higher order correction -------------------------------------"""
        # todo: please check if the proedure is resonable !!!!!!!!!

        """ 
        1. based on new energy, position, B, and v_dir, calculate new rotation angular velocity
        by the same procedure
        """
        temp_a = np.cross(v_dir_new, const_B_new)
        temp_b = np.cross(const_B_new, temp_a)

        if np.linalg.norm(temp_b) == 0 :
            cosT_new = 0
        else:
            temp_c = temp_b / np.linalg.norm(temp_b)
            cosT_new = np.dot(temp_c, v_dir_new)

        v_perp_new = v_total_new * cosT_new

        R2 = self.calR(v_dir_new, E_new, const_B_new, cosT_new)
        if R2 == 0:
            v_omg2 = 0
        else:
            v_omg2 = v_perp_new / R2

        """
        2. calculate averaged v_dir_next as new self.v_dir
           here v_omg is what we get in part 1, v_omg2 is what we get in part 2
           this works like improved euler method
           
           then new rotation matrix is constructed and correct the v_dir
        """
        v_omg_final = (v_omg + v_omg2) / 2

        omg = v_omg_final * t

        # use the new rotation angle, rotate the initial v_dir again
        M = np.array([[cos(omg) + (1 - cos(omg)) * np.power(rotaxis[0], 2),
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[1] - sin(omg) * rotaxis[2],
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[2] + sin(omg) * rotaxis[1]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[1] + sin(omg) * rotaxis[2],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[1], 2),
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] - sin(omg) * rotaxis[0]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[2] - sin(omg) * rotaxis[1],
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] + sin(omg) * rotaxis[0],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[2], 2)]])

        v_perp = np.dot(M, old_perp)

        # the new v_dir is set up, the old one is recorded to do correction for other info
        v_old_dir = self.v_dir
        self.v_dir = v_perp + v_paralleldir

        """3. correction for new position"""
        # todo: check the procedure
        # since we determine in previous part that v_dir in next step is self.v_dir and speed is v_total_new
        # in beginning they are v_old_dir and self.v_total
        # the averaged v_dir is get as new_v
        tempv = (v_old_dir * self.v_total + self.v_dir * v_total_new)

        new_v = tempv / np.linalg.norm(tempv)
        # the average speed is get from average of energy in the beginning and energy we predicted to have in the end
        avg_v = self.calV((self.E + E_new) / 2)

        # update position and pre_position
        self.prepos = self.pos
        self.pos = self.pos + new_v * avg_v * t

        # update total distance and pre_distance
        Dis = avg_v * t
        self.predis = self.totalDis
        self.totalDis = self.totalDis + Dis

        # update speed by new energy
        self.v_total = self.calV(E_new)

        self.totalAng = self.totalAng + omg

        self.finalvdir = self.v_dir  # update finalvdir (seems to be a redundant var since it is just v_dir ?)

        self.Radius = R

        self.preenergy = self.E

        # get new energy based on total distnace traveled and initial energy
        self.E = self.getnewE(self.initialE, self.totalDis)

        """4. get new update factor"""

        # get update factor according to how many pieces we want the circle to be cut into
        if v_perp_new != 0:
            timet = self.Radius * self.updateref / v_perp_new
        # if exactly parallel
        else:
            timet = 10e-5

        # we dont want time period to be so large since the muon will just go super far away out of earth in the
        # last step, then it is hard to do final position correction
        # so the largest distance to travel is limited to 1000m
        if Dis >= 1000:
            timet = 1000/self.v_total
        if Dis <= 1000 and Dis >= 900:
            timet = 1000/self.v_total
        
        return timet



"""
for a single muon, eject it and get the returned final info
"""
def main(v_dir, p, E):
    # cerate a new muon object, with initial v_dir and energy
    x = muon(v_dir, E)

    # set up initial update time, since in the first iteration we do not have the computed one
    # it should be as small as possible
    # todo: is this too small? I am not sure if this will cause any floating point issue
    t = 1./100000000000000


    while 1:
        # return value is the time period for next iteration
        # t is time and p is the file we write on
        t = x.updatepos(t, p)

        # check the end of iteration
        if np.linalg.norm(x.pos) >= 6371000 :
            break

    try:
        # get deflection angle based on the dot product of initial v_dir and final v_dir
        defang = acos(np.dot(x.initialvdir, x.finalvdir)/(np.linalg.norm(x.initialvdir) * np.linalg.norm(x.finalvdir)))
    except:
        defang = x.totalAng

    return x.totalAng, defang, x.totalDis, x.E, x.pos, x.prepos, x.predis, x.preenergy





"""
This is the start case for multiple muons with different initial v_dir and initial E
"""
def start():

    power = np.linspace(9, 12, 31)  # energy power, here is 0.1 energy step size

    E = []  # array for energy level

    for i in range(len(power)):
        E.append(pow(10, power[i]))

    E = np.array(E)

    # some print method to check E and file name
    print(np.log10(E))
    print(E/pow(10, 9))
    for j in range(0, len(E)):
      print("muon_2D_def_E_" + str(round(power[j] - 9, 1)))


    """
    get theta array and phi array
    """
    phi = np.linspace(0, 2*pi, 73)  # assign 73 phi uniformly on 0 and 2pi

    print(np.degrees(phi))  # check

    costheta = np.linspace(cos(radians(0)), cos(radians(90)), 28)   # assign costheta from 1 to 0 uniformly

    theta = np.arccos(costheta)  # convert back to theta in rad

    print(np.degrees(theta))     # check

    muon_v_dir = [[] for i in range(len(phi))]  # 2D array of initial direction of velocity

    for i in range(0, len(theta)):
        for j in range(0, len(phi)):
            # for each theta and phi combination perform a coordinate transformation to (x, y, z) of v_ini
            # notice that -1 is used to convert ejecting direction to receiving direction
            # since in the constructor we assume v_dir is the final velocity direction we received
            muon_v_dir[i].append(-1 * changecord(theta[i], phi[j]))

    # doing simulattion for all energy levels
    for j in range(0, len(E)):
        print("energy: " + str(E[j])) # print to check
        p = open("muon_2D_def_E_" + str(round(power[j] - 9, 1)) + "_pos.txt", "w")
        # doing simulation for all initial velocity directions in this energy level
        for i in range(0, len(muon_v_dir)):
            for k in range(0, len(muon_v_dir[i])):

                # get return value and write to output file
                totalang1, totalang2, totaldis, energy, pos, prepos, predis, preenergy = main(muon_v_dir[i][k], p, E[j])
                p.write(str(theta[i]) + " " + str(phi[k]) + " " + str(totalang2) + " " + str(totaldis) + " " + str(energy) + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(prepos[0]) + " " + str(prepos[1]) + " " + str(prepos[2]) + " " + str(predis) + " " + str(preenergy) + '\n')

                # print finished progress and some info to check
                print("finish theta: " + str(degrees(theta[i])) + " phi: " + str(degrees(phi[k])) + "\n")
                print(str(theta[i]) + " " + str(phi[k]) + " " + str(totalang2) + " " + str(totaldis) + " " + str(energy) + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(prepos[0]) + " " + str(prepos[1]) + " " + str(prepos[2]) + " " + str(predis) + " " + str(preenergy) + '\n')

        p.close()


# start the simulation
start()
