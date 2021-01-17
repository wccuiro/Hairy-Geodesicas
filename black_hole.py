import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math as ma


class black_hole():
    def __init__(self, alpha0, eta0, nu0):
        self.__c = 6.407e4
        self.__G = 39.309
        self.alpha = alpha0
        self.eta = eta0
        self.nu = nu0

    def masa(self):
        if self.alpha < 0:
            '''calcular masa branch positivo'''
            M = -self.__c ** 2 / self.__G * \
                (3 * self.eta ** 2 + self.alpha) / self.eta ** 3 / 6
            return M
        elif self.alpha > 0:
            # calcular masa branch negativo
            M = self.__c ** 2 / self.__G * \
                (3 * self.eta ** 2 + self.alpha) / self.eta ** 3 / 6
            return M
        else:
            print("alpha no puede ser 0")

    def horizonte_hairy_x(self):
        '''calcular x horizonte'''
        if self.alpha > 0:
            x_inicial = 0.00004
            return fsolve(self.f, x_inicial)
        elif self.alpha < 0:
            x_inicial = 26
            return fsolve(self.f, x_inicial)

    def horizonte_hairy(self):
        if self.alpha > 0:
            x_inicial = 0.00004
            x_h = fsolve(self.f, x_inicial)
            return self.omega(x_h)**0.5
        elif self.alpha < 0:
            x_inicial = 26
            x_h = fsolve(self.f, x_inicial)
            return self.omega(x_h)**0.5

    def horizonte(self):
        r_h = 2*self.__G*self.masa()/self.__c**2
        return r_h

    def omega(self, x):
        O = self.nu ** 2 * x ** (self.nu - 1) / \
            self.eta ** 2 / (x ** self.nu - 1) ** 2
        return O

    def f(self, x):
        f = self.alpha * (1 / (self.nu ** 2 - 4) - x ** 2 * (1 + x ** (-self.nu) / (self.nu - 2) - x ** self.nu / (
            self.nu + 2)) / self.nu ** 2) + x * self.eta ** 2 * (x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1)
        return f

    def diff_omega(self, x):
        dO = self.nu ** 2 * x ** (self.nu - 1) * (self.nu - 1) / x / self.eta ** 2 / (x ** self.nu - 1) ** 2 - \
            2 * self.nu ** 3 * \
            x ** (self.nu - 1) / self.eta ** 2 / \
            (x ** self.nu - 1) ** 3 * x ** self.nu / x
        return dO

    def diff_f(self, x):
        df = self.alpha * (-2 * x * (1 + x ** (-self.nu) / (self.nu - 2) - x ** self.nu / (self.nu + 2)) / self.nu ** 2 - x ** 2 * (-x ** (-self.nu) * self.nu / x / (self.nu - 2) - x ** self.nu * self.nu / x / (self.nu + 2)) / self.nu ** 2) + self.eta ** 2 * (
            x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1) + 2 * self.eta ** 2 * (x ** self.nu - 1) / self.nu / x ** (self.nu - 1) * x ** self.nu - self.eta ** 2 * (x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1) * (self.nu - 1)
        return df

    def estado(self):
        return "La masa del BH es: ", self.masa(), "\nEl horizonte hairy es: ", self.horizonte_hairy(), "\nEl horizonte tipo Schwarzschild es: ", self.horizonte()


class particula_time_like(black_hole):
    def __init__(self, alpha0, eta0, nu0, energia0, J0):
        super().__init__(alpha0, eta0, nu0)
        self.__c = 6.407e4
        self.__G = 39.309
        self.energia = energia0
        self.J = J0

    def U_potencial(self, x):
        U = self.omega(x)*self.f(x)*(1+self.J**2*self.__c**2/(self.omega(x)))-1
        return U

    def potencial(self):
        if self.alpha < 0:
            x = np.linspace(1, 100, 100000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            print(r, U)
            return r, U
        elif self.alpha > 0:
            x = np.linspace(self.horizonte_hairy(), 1, 10000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        else:
            print("alpha no puede ser 0")

    def diff_U(self, x):
        dU = self.diff_f(x)*self.omega(x)+self.diff_omega(x) * \
            self.f(x)+(self.J*self.__c)**2*self.diff_f(x)
        return dU

    def maximo(self):
        x_inicial = float(self.horizonte_hairy_x())
        x_max = fsolve(self.diff_U, x_inicial)
        U_max = self.omega(x_max)*self.f(x_max)*(1+self.J **
                                                 2*self.__c**2/(self.omega(x_max)))-1
        return self.omega(x_max)**(0.5), U_max

    def minimo(self):
        x_inicial = float(self.horizonte_hairy_x())+0.8
        x_min = fsolve(self.diff_U, x_inicial)
        U_min = self.omega(x_min)*self.f(x_min)*(1+self.J **
                                                 2*self.__c**2/(self.omega(x_min)))-1
        return self.omega(x_min)**(0.5), U_min

    def U_nulo(self):
        x_inicial = float(self.horizonte_hairy_x())+0.6
        x_nulo = fsolve(self.U_potencial, x_inicial)
        return self.omega(x_nulo)**(0.5), self.U_potencial(x_nulo)

    def cond_init(self, r):
        x_initial_guess = 0.83
        def O_r(x): return self.nu ** 2 * x ** (self.nu - 1) / \
            self.eta ** 2 / (x ** self.nu - 1) ** 2-r**2
        x_radio = fsolve(O_r, x_initial_guess)
        x_prima = ma.sqrt((self.energia-self.U_potencial(x_radio)
                           )/(self.eta**2*self.J**2*self.__c**2))
        return x_radio, x_prima

    def funH(self, x):
        H = self.diff_U(x)/(2*self.J**2*self.__c**2*self.eta**2)
        return H

    def x_turn_r(self, x):
        return self.omega(x)**0.5


class particula_null(black_hole):
    def __init__(self, alpha0, eta0, nu0, b0):
        super().__init__(alpha0, eta0, nu0)
        self.__c = 6.407e4
        self.__G = 39.309
        self.b = b0

    def U_potencial(self, x):
        U = self.__c**2*self.f(x)
        return U

    def potencial(self):
        if self.alpha < 0:
            x = np.linspace(1, 1000000, 10000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        elif self.alpha > 0:
            x = np.linspace(self.horizonte_hairy(), 1, 10000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        else:
            print("alpha no puede ser 0")

    def maximo(self):
        x_inicial = float(self.horizonte_hairy_x())
        x_max = fsolve(self.diff_f, x_inicial)
        return self.omega(x_max)**(0.5), self.U_potencial(x_max)

    def cond_init(self, r_x):
        r = ma.sqrt(r_x**2+self.b**2)
        theta = ma.atan(self.b/r_x)
        x_initial_guess = 0.8
        def O_r(x): return self.nu ** 2 * x ** (self.nu - 1) / \
            self.eta ** 2 / (x ** self.nu - 1) ** 2-r**2
        x_radio = fsolve(O_r, x_initial_guess)
        x_prima = ma.sqrt((1/self.b**2-self.f(x_radio))/(self.eta**2))
        return x_radio, x_prima, theta

    def x_turn_r(self, x):
        return self.omega(x)**0.5

    def funH(self, x):
        H = self.diff_f(x)/(2*self.eta**2)
        return H


def RK(x_n, y_n, h, s, particula):
    r = np.array([])

    # punto de referencia
    for i in range(s):
        r_n = particula.omega(x_n)**0.5
        r = np.append(r, r_n)

        k_0 = h*y_n
        l_0 = h*(-particula.funH(x_n))

        k_1 = h*(y_n+l_0/2)
        l_1 = h*(-particula.funH(x_n+0.5*k_0))

        k_2 = h*(y_n+l_1/2)
        l_2 = h*(-particula.funH(x_n+0.5*k_1))

        k_3 = h*(y_n+l_2)
        l_3 = h*(-particula.funH(x_n+k_2))

        y_n1 = y_n+(l_0+2*l_1+2*l_2+l_3)/6
        x_n1 = x_n+(k_0+2*k_1+2*k_2+k_3)/6

        y_n = y_n1
        x_n = x_n1
    return(r)


# creando todos los objetos
miBH = black_hole(1, 12.527, 1.52)  # branch negativo
miParticula = particula_time_like(1, 12.527, 1.52, -0.025, 2.6e-6)  # (bh,energia,J)


# miParticula=particula_null(1,3.252719443,1.52,0.225) #(bh,b)

# caracteristicas de los objetos para grafica de potencial
# generales
# r_apoyo=np.linspace(-5e11,5e11,5)
# r_schwarzschild=miBH.horizonte()*np.ones(len(r_apoyo))
# r_hairy=miBH.horizonte_hairy()*np.ones(len(r_apoyo))
# print(miBH.horizonte(),miBH.horizonte_hairy(),miBH.horizonte_hairy_x())

# tipo time like
# r_min,U_min=miParticula.minimo()
# r_nul,U_nul=miParticula.U_nulo()
# r_max,U_max=miParticula.maximo()
# r,U=miParticula.potencial()


# tipo null
# r_null,U_null=miFoton.potencial()
# rn_max,Un_max=miFoton.maximo()

# condiciones iniciales de los objetos para geodesicas tipo time like
x_n, y_n = miParticula.cond_init(0.4)

# geodesicas null
theta_end = 15*ma.pi
s = 10000
h = (theta_end-0)/s
theta = np.linspace(0, theta_end, s)
r = RK(x_n, -y_n, h, s, miParticula)

# print(r[9289],theta[9289])

# para que sea en polares sin necesidad de transformar
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.spines['polar'].set_visible(False)  # para quitarle el borde feo

# ax.set_yticks([])  #le quita todas las etiquetas en r
# ax.set_xticks([])  #le quita todas las etiquetas de los angulos
ax.plot(theta, r, dashes=[4, 1])

#ax.legend((l, l1, l2), (f'E = {miParticula.energia}', f'E = {miParticula1.energia}',f'E = {miParticula2.energia}'), loc='upper right', shadow=True)
ax.legend([f'E = {miParticula.energia}'], loc='upper right', shadow=True)
# ax.set_rmax(r.max())
# ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
# ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
# ax.grid(True)

#ax.set_title("Desviación de Einstein", va='bottom')


#ax.plot(theta, r)
# ax.grid(True)
#fig, ax = plt.subplots()
BH_schwarzschild = plt.Circle(
    (0, 0), miBH.horizonte(), transform=ax.transData._b, color='dimgrey')
BH_hairy = plt.Circle((0, 0), miBH.horizonte_hairy(),
                      transform=ax.transData._b, color='black')

# x=r*np.cos(theta)
# y=r*np.sin(theta)

# ax.set_xlim(-0.6,0.6)
# ax.set_ylim(-0.6,0.6)
ax.add_artist(BH_schwarzschild)
ax.add_artist(BH_hairy)
# ax.plot(x,y)

# fill para time like
#plt.fill_between(r_apoyo, U_max, 50,facecolor='aquamarine')
#plt.fill_between(r_apoyo, U_nul, U_max,facecolor='violet')
#plt.fill_between(r_apoyo, U_min, U_nul,facecolor='greenyellow')
#plt.fill_between(r_apoyo, -50, U_min,facecolor='coral')

# fill para null
#plt.fill_between(r_apoyo, -5e11, 0,facecolor='aquamarine')
#plt.fill_between(r_apoyo, 0, Un_max,facecolor='violet')
#plt.fill_between(r_apoyo, Un_max, 50e11,facecolor='greenyellow')

# plt.xlim(0,1)      #lim para null
# plt.ylim(-0.1e11,1.14e11)   #lim para null

ax.set_rmax(3)  # lim para time_like
ax.set_rticks([0.08, 1, 2, 3])  # Less radial ticks
ax.set_xticks([])
#ax.set_xticks([0, theta[9289]-2*ma.pi])
# plt.ylim(-0.4,1.4)   #lim para null_like
#ax.set_theta_zero_location("E", offset=-(theta[9289]/2+np.pi)/np.pi*180)
# plt.plot(r,U,r_schwarzschild,r_apoyo,r_hairy,r_apoyo)
# plt.plot(r_null,U_null,r_schwarzschild,r_apoyo,r_hairy,r_apoyo)
ax.set_rlabel_position(-50)
ax.grid(True)

if isinstance(miParticula, particula_time_like):
    plt.savefig(
        f'timeLike_alpha-{miParticula.alpha}_eta-{miParticula.eta}_nu-{miParticula.nu}_E-{miParticula.energia}_J-{miParticula.J}.svg')

elif isinstance(miParticula, particula_null):
    plt.savefig(
        f'null_alpha-{miParticula.alpha}_eta-{miParticula.eta}_nu-{miParticula.nu}_b-{miParticula.b}.svg')
