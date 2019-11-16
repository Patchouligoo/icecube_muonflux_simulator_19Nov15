import matplotlib.pyplot as plt
from math import *
import numpy as np
np.set_printoptions(threshold=np.nan)


"""------------------------------------set front size----------------------------------------"""
SMALL_SIZE = 8 * 2
MEDIUM_SIZE = 10 * 2
BIGGER_SIZE = 12 * 2

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


"""------------------------------------------------------------------------------------------"""



f = open("total_flux_0.1.txt", 'r')

numtheta = 28
numphi = 73

# get x and y range
phi = np.linspace(0, 2*pi, 73)
costheta = np.linspace(1, cos(radians(90)), 28)

# read and store the output
totalflux = []
for i in range(0, numtheta):
    temp_arr = f.readline().split()
    arrphi = []
    for j in range(0, numphi):
        arrphi.append(float(temp_arr[j]))

    totalflux.append(np.array(arrphi))


totalflux_ori = totalflux   # copy of original data

# doing row normalization
totalflux = np.array(totalflux)
for i in range(0, numtheta):
    totalflux[i] = totalflux[i]/np.average(totalflux[i])


"""   plotting the relative flux  """
fig = plt.figure()
ax = fig.add_subplot(111)
f = plt.imshow(totalflux, extent=(np.amin(phi), np.amax(phi), np.amin(costheta[0:numtheta]), np.amax(costheta[0:numtheta])), aspect = 'auto')
plt.xlabel('phi')
plt.ylabel('cos theta')
plt.colorbar(f)
ax.set_title("relative flux, all E")
plt.clim(0.999, 1.001)


"""plot the log scale absolute flux map  """
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
f2 = plt.imshow(np.log10(totalflux_ori), extent=(np.amin(phi), np.amax(phi), np.amin(costheta[0:numtheta]), np.amax(costheta[0:numtheta])), aspect = 'auto')
plt.xlabel('phi')
plt.ylabel('cos theta')
plt.colorbar(f2)
ax2.set_title("total flux, all E, log scale")


""" plot for 1D graphs of flux"""
oneDplotlist = []
counterlist = []

for i in range(len(totalflux)):
    if i%5 is 0:
    #if i >= 25:
        counterlist.append(i)
        oneDplotlist.append(totalflux[i])

fig = plt.figure()

for i in range(len(oneDplotlist)):
    plt.plot(phi, oneDplotlist[i], label="theta = " + str(np.round(np.degrees(np.arccos(costheta[counterlist[i]])), 2)))

plt.legend()
plt.xlabel("phi")
plt.ylabel("relative flux ratio")
plt.show()












"""-----------------------some methods to check the integration results, NOT USED ------------------------------"""


"""
def getnewE_fromintegral(E, dis):
    #print("integral " + str(E) + str(dis))
    b = 0.000363

    a = 0.259

    density = 0.9167

    A = (E / pow(10, 9) * b * density) + a * density
    # print(exp(dis))

    E = (A * np.power(e, b * dis * density) - a * density) / (b * density) * pow(10, 9)

    return E


def getflux():
    costheta = np.linspace(1, cos(radians(90)), 28)

    theta = np.arccos(costheta)

    #print(np.degrees(theta))
    r = 6371000 - 2000
    R = 6371000

    iniE = np.logspace(9, 12, 1000)

    iniE_log = np.log10(iniE)

    pre_flux = []

    flux_total = [0 for i in range(len(theta))]

    check = 0

    for j in range(len(iniE)):
        flux = []
        flux2 = []
        for i in range(0, len(theta)):
            thetaideal = theta[i]
            dis = sqrt(pow(R, 2) - pow(r, 2) * pow(sin(pi - thetaideal), 2)) + r * cos(pi - thetaideal)
            E = getnewE_fromintegral(iniE[j], dis)

            # todo: add integration of all E ??????

            cosa = (pow(dis, 2) + pow(R, 2) - pow(r, 2)) / (2 * dis * R)

            flux_temp = pow(E / pow(10, 13), -3.7) * pow(cosa, 2)

            flux.append(flux_temp)

        if check != 0:

            for n in range(len(flux)):
                height = 0.5 * (pre_flux[n] + flux[n])
                dx = iniE_log[j] - iniE_log[j - 1]

                area = height * dx

                result = area / log10(e) * (pow(10, iniE_log[j]) + pow(10, iniE_log[j - 1])) / 2

                flux_total[n] = flux_total[n] + result


        check = 1
        pre_flux = flux

    return costheta, flux_total
"""


"""-------------------------------------------------------------------------"""

"""
fixed_theta = []

for i in range(numtheta):
    flux = 0
    for j in range(numphi):
        flux = flux + totalflux_ori[i][j]
    fixed_theta.append(flux)


fig = plt.figure()

#plt.plot(costheta[0:numtheta], np.array(fixed_theta), label='simulated data')


plt.plot(getflux()[0], np.array(fixed_theta)/(np.array(getflux()[1]) * 2.4 * numphi), label='ratio between real and simulation')

print(fixed_theta)
#plt.yscale('log')
plt.xlabel('cos theta')
plt.ylabel('flux')

plt.legend()

plt.show()
"""
