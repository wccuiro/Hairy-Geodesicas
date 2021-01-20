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

# tipo null
r_null,U_null=miParticula.potencial()
rn_max,Un_max=miParticula.maximo()

# fill para null
plt.fill_between(r_apoyo, -5e11, 0,facecolor='aquamarine')
plt.fill_between(r_apoyo, 0, Un_max,facecolor='violet')
plt.fill_between(r_apoyo, Un_max, 50e11,facecolor='greenyellow')

plt.xlim(0,1)      #lim para null
plt.ylim(-0.1e11,1.14e11)   #lim para null

plt.ylim(-0.4,1.4)   #lim para null_like

plt.plot(r_null,U_null,r_schwarzschild,r_apoyo,r_hairy,r_apoyo)

plt.grid(True)

plt.show()