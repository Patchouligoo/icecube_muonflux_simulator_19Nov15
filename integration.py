from math import *
import numpy as np
np.set_printoptions(threshold=np.nan)

"""
===================================file header======================================
This method reads data from corrected flux result from correction.py and doing a 
numerical integration based on energy to get the total flux of all energy levels

The output is the total flux for all energy range, and it is written to another
file and will get plotted in plot.py
=======================================end==========================================
"""

# get the energy list in log scale
delta = 0.1
maxpow = 12
numofsplit = int((maxpow - 9)/delta + 1)
E = np.logspace(9, 12, numofsplit)
E = np.log10(E)

print(E)

# total flux will be one 28 x 73 matrix
total_flux = []
counter = 0

# (0, 1) if we want pos only and (0, 2) if pos and neg
# first deal with pos and then neg
for i in range(0, 2):

    # the algorithm is that for each (theta, phi) pair we open all the files with different energy and doing
    # numerical integration one by one only for it. After finishing one pair we go for another
    for theta in range(0, 28):
        phi_totalflux = []
        for phi in range(0, 73):
            print("doing: " + str(counter))
            counter = counter + 1
            flux_array = []

            """read all the flux corresponding to this (phi, theta) pair and append it in an array"""
            for j in range(0, len(E)):
                # if i is one then we are dealing with pos, else neg
                if i is 0:
                    f = open("muon_pos_" + str(j) + ".txt", 'r')
                else:
                    f = open("muon_neg_" + str(j) + ".txt", 'r')

                # locate the current flux we are looking at
                temp_theta = 0
                while temp_theta < theta:
                    f.readline()
                    temp_theta = temp_theta + 1
                temp = float(f.readline().split()[phi])

                flux_array.append(temp)

                f.close()

            """now do the integration, result is the total flux for all energy range on this (theta, phi) pair"""
            result = 0
            for k in range(0, len(E) - 1):
                height = 0.5 * (flux_array[k] + flux_array[k + 1])
                dx = E[k + 1] - E[k]
                area = height * dx
                result = result + area/log10(e) * (pow(10, E[k]) + pow(10, E[k + 1]))/2

            # append it to the list representing total flux of a given theta, which is the row for final output
            phi_totalflux.append(result)

        # append all the pos results
        if i == 0:
            total_flux.append(phi_totalflux)
        # add the neg result to the pos result
        else:
            total_flux[theta] = np.array(total_flux[theta]) + np.array(phi_totalflux)

# write to output file
p = open("total_flux_0.1.txt", 'w')
for i in range(0, len(total_flux)):
    for j in range(0, len(total_flux[i])):
        p.write(str(total_flux[i][j]) + " ")

    p.write('\n')

p.close()