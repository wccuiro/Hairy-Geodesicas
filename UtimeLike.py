import black_hole as bh
import matplotlib.pyplot as plt
import numpy as np
import math as ma

import os

# creando todos los objetos
miBH = bh.black_hole(-40, 3.252719443, 1.76)  # branch negativo
miParticula = bh.particula_time_like(-40,3.252719443,1.76,0, 7e-7) #(bh,b)

# caracteristicas de los objetos para grafica de potencial
r_apoyo=np.linspace(-5e11,5e11,5)
r_schwarzschild=miBH.horizonte()*np.ones(len(r_apoyo))
r_hairy=miBH.horizonte_hairy()*np.ones(len(r_apoyo))
print(miBH.horizonte(),miBH.horizonte_hairy(),miBH.horizonte_hairy_x())

# tipo time like
r_min,U_min=miParticula.minimo()
r_nul,U_nul=miParticula.U_nulo()
r_max,U_max=miParticula.maximo()
r,U=miParticula.potencial()

# # fill para time like
# plt.fill_between(r_apoyo, U_max, 50,facecolor='aquamarine')
# plt.fill_between(r_apoyo, U_nul, U_max,facecolor='violet')
# plt.fill_between(r_apoyo, U_min, U_nul,facecolor='greenyellow')
# plt.fill_between(r_apoyo, -50, U_min,facecolor='coral')

plt.plot(r,U,r_schwarzschild,r_apoyo,r_hairy,r_apoyo)

plt.xlim(0,0.5)
plt.ylim(-0.5,1.5)

plt.grid(True)


try:
    os.mkdir('./Figures')
except FileExistsError:
    pass
try:
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}')
except FileExistsError:
    pass
try:
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Time_Like_Potencial')
except FileExistsError:
    pass

plt.savefig(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Time_Like_Potencial/J={miParticula.J}.svg')