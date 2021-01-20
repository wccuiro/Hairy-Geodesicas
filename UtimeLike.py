import black_hole as bh
import matplotlib.pyplot as plt
import numpy as np

# creando todos los objetos
miBH = bh.black_hole(1, 12.527, 1.52)  # branch negativo
miParticula = bh.particula_null(1,3.252719443,1.52,0.225) #(bh,b)

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

# fill para time like
plt.fill_between(r_apoyo, U_max, 50,facecolor='aquamarine')
plt.fill_between(r_apoyo, U_nul, U_max,facecolor='violet')
plt.fill_between(r_apoyo, U_min, U_nul,facecolor='greenyellow')
plt.fill_between(r_apoyo, -50, U_min,facecolor='coral')

plt.plot(r,U,r_schwarzschild,r_apoyo,r_hairy,r_apoyo)

plt.grid(True)


try:
    os.mkdir('./Figures')
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}')
    os.mkdir(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Time_Like_Potencial')
except FileExistsError:
    print("Listo")

plt.savefig(f'./Figures/alpha={miParticula.alpha}_eta={miParticula.eta}_nu={miParticula.nu}/Time_Like_Potencial/J={miParticula.J}.svg')