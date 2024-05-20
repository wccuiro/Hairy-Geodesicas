import black_hole as bh
import matplotlib.pyplot as plt
import numpy as np
import math as ma
import matplotlib.font_manager

import os

#fig, ax = plt.subplots()

# creando todos los objetos
miBH = bh.black_hole(1, 12.52655373, 1.52)  # branch negativo
miParticula = bh.particula_null(1, 12.52655373, 1.52, 0.1) #(bh,b)

# caracteristicas de los objetos para grafica de potencial
r_apoyo=np.linspace(-5e11,5e11,5)
r_schwarzschild=miBH.horizonte()*np.ones(len(r_apoyo))
r_hairy=miBH.horizonte_hairy()*np.ones(len(r_apoyo))
print(miBH.horizonte(),miBH.horizonte_hairy(),miBH.horizonte_hairy_x())

# tipo null
r_null,U_null=miParticula.potencial()
rn_max,Un_max=miParticula.maximo()
print(rn_max, Un_max)

## fill para null
#plt.fill_between(r_apoyo, -5e11, 0,facecolor='aquamarine')
#plt.fill_between(r_apoyo, 0, Un_max,facecolor='violet')
#plt.fill_between(r_apoyo, Un_max, 50e11,facecolor='greenyellow')

#U maximo
plt.plot(r_apoyo, Un_max*np.ones(len(r_apoyo)), color='blue', dashes=[7,2])
plt.plot(rn_max, Un_max, color='black', marker='D')



plt.xlim(0,1)      #lim para null
plt.ylim(-0.1e11,1.5e11)   #lim para null

plt.xlabel('r',labelpad=-5, fontsize=25)
plt.ylabel(r'$\mathcal{V}$',labelpad=15, fontsize=25, rotation='horizontal', position=(0,0.46))

plt.plot(r_null, U_null, color='red', linewidth=1)
#plt.plot(r_schwarzschild,r_apoyo)
plt.plot(r_hairy, r_apoyo,  dashes=[4, 1], linewidth=3, color='black')

plt.annotate('Event Horizon',
            xy=(65, 290), xycoords='figure points', fontsize=10, fontstyle='italic', fontfamily='serif')

plt.annotate('Region - I',
            xy=(190, 270), xycoords='figure points', fontsize=20, fontstyle='italic', fontfamily='serif')

plt.annotate('Region - II',
            xy=(190, 170), xycoords='figure points', fontsize=20, fontstyle='italic', fontfamily='serif')

plt.title(r'$\dfrac{G_{N}M}{c^2}=0.04 [AU],\alpha = 1 [AU]^{-2}, \eta = 12.527 [AU]^{-1}$', pad=20)

#plt.grid(True)

try:
    os.mkdir('./Figures')
except:
    pass
try:
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}')
except FileExistsError:
    print("Listo")

plt.savefig(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Null_Potencial.svg')