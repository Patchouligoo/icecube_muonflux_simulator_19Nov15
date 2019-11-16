from math import *
import numpy as np
from scipy.interpolate import UnivariateSpline
np.set_printoptions(threshold=np.nan)

"""
===================================file header======================================
The file generated by the simulation has only position, direction, energy info, and 
so on and the final position is highly inaccurate, because in the last step the muon
will go far away beyond earth surface which causes error in E and distance

So this file will read all files generated by simulator, doing correction on final
position and energy, then for each file (each file corresponds to an energy level
with 28 rows and 73 columns representing different initial angle), the final dir and
energy will be converted to flux

The result for each file is a 28 x 73 matrix representing flux corresponds to the 
total flux map of a certain initial energy level. for energy step 0.1, 31 files will
be generated.

Then they will be used to do numerical integration in integration.py
=======================================end==========================================
"""




"""
gget energy loss per cm : dE/dx using :
dE/dx = 0.259GeV/wme + 0.363/wme * E(GeV)
dE/dx is eV per cm
"""
def getEloss(E):

    E = E / pow(10, 9)
    totalloss = (0.000363 * E + 0.259) / 100 * 0.9167 * pow(10, 9)


    return totalloss



"""
Get positive and negative muon ratio according to arXiv:1206.6710v2

In this simulation the flux of negative muon is regarded as 1, and positive muon flux
is icreased by a factor according to their energy
This method will generate the ratio factor for positive muon
"""
def getratio(E):
    E = E/pow(10, 9)

    if E > 12 * pow(10, 4):
        return 1.4

    x = [97.82625165318657, 133.7675299099177, 200.77823789573986, 275.72398329251484, 352.1110743386039,
         476.5383449564929, 635.6351053049032,
         847.8271150190255, 1130.7166822094173, 1508.3206821038357, 2011.9989657768444, 2663.46265783033,
         4027.8127093079784, 5166.515006984569,
         7653.189919920697, 10210.089678867402, 13622.670782258743, 18175.964940003934, 24249.15815739906,
         32355.013124451427, 43166.44787435703,
         57592.638961903715, 76834.79570848323, 102514.57239865992, 119973.19012320979]
    y = [1.2611490711739954, 1.2677263867121407, 1.2750750472164714, 1.2804580767900025, 1.288172621958387,
         1.2983324829410647, 1.3073741724173016, 1.3172274171282663,
         1.3309291522102389, 1.3380521941185124, 1.3456334407705683, 1.3577839825843034, 1.3653796334056616,
         1.3677881117146828, 1.3778173204886548, 1.3816894036440746,
         1.3821837526278595, 1.3824920444548021, 1.3854757462098406, 1.385063818608925, 1.3876403768935077,
         1.3890337107410278, 1.3926461654031443, 1.3935212709862572, 1.3945520928716182]

    spl = UnivariateSpline(x, y, k= 4 ,s= 0.0001)

    return spl(E)



"""
This method return the new energy for initial energy E and new distance traveled forward as dis

the energy loss per cm is dE/dx, which we get from getEloss(E)
a simple forward euler method is performed to get new energy after dis
since the dis must be less than 1000m due to our set up in simulator, it should be accurate
"""
def getnewE_p(E, dis):

    deltax = 0.1
    num = int(dis/deltax) + 10
    deltax = dis/num


    for i in range(0, num):
        E_loss1 = getEloss(E)
        E_temp = E + 100 * deltax * E_loss1
        E_loss2 = getEloss(E_temp)
        E_loss = (E_loss1 + E_loss2)/2

        E = E + 100 * deltax * E_loss
    return E

"""
Same as previous one, but this time we get PREVIOUS !!! energy for energy E and dis traveled BACKWARD !!!

which means we over run during the correction and need to go back
"""
def getnewE_m(E, dis):
    deltax = 0.1
    num = int(dis / deltax) + 10
    deltax = dis / num


    for i in range(0, num):
        E_loss1 = getEloss(E)
        E_temp = E - 100 * deltax * E_loss1
        E_loss2 = getEloss(E_temp)
        E_loss = (E_loss1 + E_loss2)/2

        E = E - 100 * deltax * E_loss


    return E


"""============================================main body==============================================="""

numtheta = 28
numphi = 73

# used to control how much rows we process
numberofthetainterest = numtheta

# used to count how much files for pos and neg muons we have processed, used to set output file name
counter_pos = 0
counter_neg = 0

interval = 0.1  # energy step in log

# from 0 to 61, 31 for positive muons and 31 for negative muons
for num in range(0, 2 * int((12 - 9)/interval + 1)):

    check = 0   # used to determine if we are dealing with positive or negative muon file to change flux ratio

    # we deal with all positive muon files first ! 31 iterations
    if num < int((12 - 9)/interval + 1):
        check = 1
        # open file and output file
        f = open("./file/muon_2D_def_E_" + str(round(num  * interval, 1)) + "_pos.txt", "r")
        p = open("muon_pos_" + str(counter_pos) + ".txt", 'w')
        counter_pos = counter_pos + 1
    # then we deal with the negative muon file
    else:
        check = 2
        f = open("./file/muon_2D_def_E_" + str(round((num - int((12 - 9)/interval + 1)) * interval, 1)) + "_neg.txt", "r")
        p = open("muon_neg_" + str(counter_neg) + ".txt", 'w')
        counter_neg = counter_neg + 1

    # 2d list to hold the 28 x 73 matrix of muon flux
    array_phi_theta_flux = []

    # this loop iterates through all the (theta, phi) pairs in this file, for each of them correction
    # on final posiiton and energy is performed and muon flux is obtained
    for i in range(0, numtheta):
        arrayflux = []
        c = 0

        if i > numberofthetainterest:
            break
        else:
            for j in range(0, numphi):

                temp = f.readline().split()

                # get all the useful info from file, remember each line in that file includes the pos, prepos,
                # energy, distance ..... all the info for a muon with certain initial firection
                dis = float(temp[3])
                predis = float(temp[11])
                pos = np.array([float(temp[5]), float(temp[6]), float(temp[7])])
                prepos = np.array([float(temp[8]), float(temp[9]), float(temp[10])])
                preenergy = float(temp[12])

                """-----------------------correction of final step--------------------------"""
                norm = np.linalg.norm(pos - prepos)     # length of final step
                dir = (pos - prepos) / norm             # direction unit vector of final step

                try:
                    # this algorithm works like a binary search. Until the difference between radius r of final
                    # position and radius of earth is less than 0.001m, the algorithm will use current step size / 2
                    # as step for next time. If difference is greater than 0, then we go back, if it is less than 0
                    # then we go forward.
                    step = norm / 2

                    while abs(np.linalg.norm(prepos) - 6371000) > 0.001:

                        # go back! we get out of earth surface too far
                        if np.linalg.norm(prepos) - 6371000 > 0:
                            prepos = prepos - dir * step
                            predis = predis - step
                            preenergy = getnewE_m(preenergy, step)  # decrease energy using FIRST method
                            step = step / 2

                        # keep going! we are still inside the earth
                        else:
                            prepos = prepos + dir * step
                            predis = predis + step
                            preenergy = getnewE_p(preenergy, step)  # increase energy using SECOND method
                            step = step / 2

                    print("theta: " + str(degrees(float(temp[0]))) + " E: " + str(preenergy) + " dis: " + str(predis))

                    # in the end, prepos is the correct position and preenergy is the correct energy
                    pos = prepos
                    energy = preenergy
                except:
                    energy = 0
                    print("error! in " + str(j))
                    exit()
                """----------------------------------end------------------------------------"""

                # get the normal vector to the surface of earth
                pos_normal = pos / np.linalg.norm(pos)
                # get the cos theta, where theta is with respect to the normal vector
                coss = np.dot(dir, pos_normal) / (np.linalg.norm(dir) * np.linalg.norm(pos_normal))

                # if it is positive muon, there is an extra getratio()
                if check == 1:
                    Eflux = getratio(energy) * pow(energy / pow(10, 13), -3.7)
                # if it is negative muon
                else:
                    Eflux = pow(energy / pow(10, 13), -3.7)
                # add the angular distribution
                arrayflux.append(pow(coss, 2) * Eflux)

            """append np array to 2D list"""
            array_phi_theta_flux.append(np.array(arrayflux))

    # write the output files
    for i in range(0, len(array_phi_theta_flux)):
        for j in range(0, len(array_phi_theta_flux[i])):
            p.write(str(array_phi_theta_flux[i][j]) + " ")

        p.write('\n')

    f.close()
    p.close()




