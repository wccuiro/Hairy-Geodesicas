import black_hole as bh
import matplotlib.pyplot as plt
import math as ma
import numpy as np

import os

# creando todos los objetos
miBH = bh.black_hole(1, 12.52655373, 1.4)  # branch negativo
miParticula = bh.particula_time_like(1, 12.52655373, 1.4, 1, 2.6e-6)  # (bh,energia,J)
# miParticula1 = bh.particula_time_like(-40, 3.252719443, 1.76, 5, 7e-7)  # (bh,energia,J)
# miParticula2 = bh.particula_time_like(-40, 3.252719443, 1.76, 10, 7e-7)  # (bh,energia,J)

#print(miBH.horizonte(),miBH.horizonte_hairy(),miBH.horizonte_hairy_x())
print(miBH.masa())

# condiciones iniciales de los objetos para geodesicas tipo time like
x_n, y_n = miParticula.cond_init(2)
print(x_n, y_n)

# x_n1, y_n1 = miParticula1.cond_init(2)
# x_n2, y_n2 = miParticula2.cond_init(2)

# geodesicas null
theta_end = 1.5*ma.pi
s = 10000
h = (theta_end-0)/s
theta = np.linspace(0, theta_end, s)
r = bh.RK(x_n, -y_n, h, s, miParticula)
print(r)

# r1 = bh.RK(x_n1, y_n1, h, s, miParticula1)
# r2 = bh.RK(x_n2, y_n2, h, s, miParticula2)


# para que sea en polares sin necesidad de transformar
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.spines['polar'].set_visible(False)  # para quitarle el borde

ax.plot(theta, r, dashes=[4, 1])

# l1, = ax.plot(theta, r, dashes=[4, 1])
# l2, = ax.plot(theta, r1, dashes=[2, 1])
# l3, = ax.plot(theta, r2, dashes=[1, 1])

ax.legend([f'E = {miParticula.energia}'], loc='upper right', shadow=True)
# ax.legend((l1, l2, l3), (f'E = {miParticula.energia}', f'E = {miParticula1.energia}', f'E = {miParticula2.energia}'), loc='upper right', shadow=True)

BH_schwarzschild = plt.Circle(
    (0, 0), miBH.horizonte(), transform=ax.transData._b, color='dimgrey')
BH_hairy = plt.Circle((0, 0), miBH.horizonte_hairy(),
                      transform=ax.transData._b, color='black')

ax.add_artist(BH_schwarzschild)
ax.add_artist(BH_hairy)

ax.set_rmax(0.3)  # lim para time_like
ax.set_rticks([0.08, 0.15, 0.3])  # Less radial ticks
ax.set_xticks([])

#ax.set_theta_zero_location("E", offset=-(theta[9289]/2+np.pi)/np.pi*180)

ax.set_rlabel_position(100)
ax.grid(True)

print(r.max())



try:
    os.mkdir('./Figures')
except FileExistsError:
    pass
try:
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}')
except FileExistsError:
    pass
try:
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/timeLike_orbits')
except FileExistsError:
    print("Listo")

plt.savefig(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/timeLike_orbits/E={miParticula.energia}_J={miParticula.J}.svg')
