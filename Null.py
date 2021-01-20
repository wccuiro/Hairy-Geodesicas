import black_hole as bh
import matplotlib.pyplot as plt
import numpy as np
import math as ma

import os

# creando todos los objetos
miBH = bh.black_hole(1, 12.527, 1.52)  # branch negativo
miParticula = bh.particula_null(1, 12.527, 1.52, -0.025, 2.6e-6)  # (bh,energia,J)

# condiciones iniciales de los objetos para geodesicas tipo time like
x_n, y_n = miParticula.cond_init(0.4)

# geodesicas null
theta_end = 8.5*ma.pi
s = 10000
h = (theta_end-0)/s
theta = np.linspace(0, theta_end, s)
r = bh.RK(x_n, -y_n, h, s, miParticula)

# para que sea en polares sin necesidad de transformar
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.spines['polar'].set_visible(False)  # para quitarle el borde

ax.plot(theta, r, dashes=[4, 1])

#ax.legend((l, l1, l2), (f'E = {miParticula.energia}', f'E = {miParticula1.energia}',f'E = {miParticula2.energia}'), loc='upper right', shadow=True)
ax.legend([f'E = {miParticula.energia}'], loc='upper right', shadow=True)


#ax.set_title("Desviación de Einstein", va='bottom')

BH_schwarzschild = plt.Circle(
    (0, 0), miBH.horizonte(), transform=ax.transData._b, color='dimgrey')
BH_hairy = plt.Circle((0, 0), miBH.horizonte_hairy(),
                      transform=ax.transData._b, color='black')


ax.add_artist(BH_schwarzschild)
ax.add_artist(BH_hairy)

ax.set_rmax(3)  # lim para time_like
ax.set_rticks([0.08, 1, 2, 3])  # Less radial ticks
ax.set_xticks([])



ax.set_rlabel_position(-50)
ax.grid(True)

try:
    os.mkdir('./Figures')
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}')
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Null_orbits')
except FileExistsError:
    print("Listo")

plt.savefig(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Null_orbits/b={miParticula.b}.svg')
