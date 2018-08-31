#!/usr/bin/python3


# NOMBRE DEL ARCHIVO CON LOS DATOS DE LA HEXPIDER
filename = "hexpider.txt"

# CARGA LIBRERIAS PARA LEER DATOS Y DIBUJAR
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
from matplotlib.patches import FancyArrowPatch




# DEFINE LA FUNCION PARA DIBUJAR EN EL TORO
def torify(x, y):
    x -= x[0] - x[0]%1.0
    y -= y[0] - y[0]%1.0

    retX = [x[0]]
    retY = [y[0]]

    i = 1
    while i < len(x):
        if (0.0 <= x[i] < 1.0) and (0.0 <= y[i] < 1.0):
            retX.append(x[i])
            retY.append(y[i])

        elif (x[i] < 0.0) and (0.0 <= y[i] < 1.0):
            transx = x[i] - x[i]%1.0
            lam = 0
            if x[i] != retX[-1]: # They must be 0.0
                lam = x[i]/(x[i] - retX[-1])

            retX.append(0.0)
            retY.append(retY[-1] * lam + (1-lam)*y[i])
            retX.append(np.nan)
            retY.append(np.nan)
            retX.append(1.0)
            retY.append(retY[-1] * lam + (1-lam)*y[i])

            x[i:] -= transx
        elif (x[i] >= 1.0 ) and (0.0 <= y[i] < 1.0):
            transx = x[i] - x[i]%1.0
            lam = 0
            if x[i] != retX[-1]: # They must be 1.0
                lam = (1-x[i])/(retX[-1] - x[i])

            retX.append(1.0)
            retY.append(retY[-1] * lam + (1-lam)*y[i])
            retX.append(np.nan)
            retY.append(np.nan)
            retX.append(0.0)
            retY.append(retY[-1] * lam + (1-lam)*y[i])

            x[i:] -= transx

        elif (y[i] < 0.0) and (0.0 <= x[i] < 1.0):
            transx = y[i] - y[i]%1.0
            lam = 0
            if retY[-1] != y[i]:
                lam = y[i]/(y[i] - retY[-1])

            retY.append(0.0)
            retX.append(retX[-1] * lam + (1-lam)*x[i])
            retX.append(np.nan)
            retY.append(np.nan)
            retY.append(1.0)
            retX.append(retX[-1] * lam + (1-lam)*x[i])

            y[i:] -= transx

        elif (y[i] >= 1.0 ) and (0.0 <= x[i] < 1.0):
            transx = y[i] - y[i]%1.0
            lam = 0
            if retY[-1] != y[i]:
                lam = (1-y[i])/(retY[-1] - y[i])

            retY.append(1.0)
            retX.append(retX[-1] * lam + (1-lam)*x[i])
            retX.append(np.nan)
            retY.append(np.nan)
            retY.append(0.0)
            retX.append(retX[-1] * lam + (1-lam)*x[i])
            y[i:] -= transx
        else:
            print("shit")
        i += 1

    return np.array(retX), np.array(retY)


# CODIGO QUE LEE Y DIBUJA
filename = "output/hexpider.txt"

with open(filename, 'r') as infile:
    desfases = set( float(l.strip().split('\t')[0]) for l in infile.readlines() )
print(desfases)
for r in desfases:
    with open(filename, 'r') as infile:
        for line in [line.strip().split('\t') for line in infile.readlines()]:
            rr,x,y = tuple(map(float, line[:3]))
            if rr!=r:
                continue;
            line = line[2:]
            neurons = np.array(list(int(n) for n in line[::2]), dtype=int)
            events = np.array(list(float(n) for n in line[1::2]))

            # Event times
            phi1 = events[neurons == 1]
            phi2 = events[neurons == 0]
            phi3 = events[neurons == 2]

            # keep minimim size
            m = min( len(phi1), len(phi2), len(phi3) )
            phi2 = phi2[:m] - phi1[:m]
            phi3 = phi3[:m] - phi1[:m]

            # Inter event times
            activations = phi1[1:] - phi1[:-1]

            m = min(m, len(activations))

            d12 = phi2[:m]/activations[:m]
            d13 = phi3[:m]/activations[:m]

            d12, d13 = torify(d12, d13)
            base_line = plt.plot(d12, d13, c='k', alpha=0.1)
            #plt.plot(d12[-1], d13[-1], 'o', c='k')#c=base_line[0].get_color())

plt.axis('scaled')

plt.xlim([0,1])
plt.ylim([0,1])
# plt.savefig('/Users/alvaro/{}.png'.format(r))
plt.show()
