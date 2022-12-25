"""
Codigo: Luis Lucas García
Grado en Física - Prácticas de Mecánica Newtoniana
Práctica 2 - Órbitas
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Variables de plt y otros

size=(6, 6)
ql = 1200
sizeU = (6, 4)
sizeM = (6, 6)

#Constantes físicas

G = 6.67e-11
M = 5.972e24
m = 200.0
k = 1
g = 9.81
m2 = 1

"""

Resolución de las ecuaciones de movimiento

"""

def grav(z, t):
    
    z1, z2, z3, z4 = z
    
    dzdt = [z3, z4, -G*M*z1/np.sqrt(z1**2 + z2**2)**3,
            -G*M*z2/np.sqrt(z1**2 + z2**2)**3]
    
    return dzdt

#Uso de la rutina odeint

x_0 = 6670000.0
y_0 = 0.0
vx_0 = 0.0
vy_0 = 1.*np.sqrt(G*M/x_0)

z0 = [x_0, y_0, vx_0, vy_0]
nt = np.linspace(0, 100000, 25000)

z = odeint(grav, z0, nt)

x = z[:,0]
y = z[:,1]
rMedio = sum((x**2 + y**2)**0.5)/len(x)
print(rMedio/1000)

#Resolución de forma analítica

#Obtenemos el ángulo entre v y r

vVect = [vx_0, vy_0]
rVect = [x_0, y_0]
angulo = np.arccos(np.dot(rVect, vVect)/(((x_0**2 + y_0**2)**0.5)+((vx_0**2 + vy_0**2)**0.5)))

#Resolvemos analíticamente

n = m*M/(m + M)
L = m*((x_0**2 + y_0**2)**0.5)*((vx_0**2 + vy_0**2)**0.5)*np.sin(angulo)
p = L**2/(G*M*m*n)
E = -G*M*m/np.sqrt(x_0**2 + y_0**2) + 0.5*m*(vx_0**2 + vy_0**2)
e = (abs(1 + (2*E*L**2)/(n*(G*M*m)**2)))**0.5

def gravA(theta):
    
    return p/(1 + e*np.cos(theta))

thetas = np.linspace(0, 2*np.pi, 2000)
r = np.array([gravA(i) for i in thetas])

xA = r*np.cos(thetas)
yA = r*np.sin(thetas)
rMedioA = sum((xA**2 + yA**2)**0.5)/len(xA)
print(rMedioA/1000)

#Gráficas

plt.figure(figsize=size)
plt.plot(x, y, label="Solución numérica", linewidth=3)
plt.plot(0, 0, "*r")
plt.plot(xA, yA, "--g", label="Solución analítica")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.legend(loc="lower left")
#plt.savefig("Circ.png", dpi=ql)

#Animación

fig, ax = plt.subplots()
ax.set_xlim(min(x)-1e6, max(x)+1e6)
ax.set_ylim(min(y)-1e6, max(y)+1e6)
ax.plot(0, 0, "*r")
line, = ax.plot(x[:1], y[:1])
a = 100


def Orb(i):
    
    line.set_data(x[:a*i+1], y[:a*i+1])
    
    return line

ani = animation.FuncAnimation(fig, Orb, frames=len(x)//a, interval=5)
#ani.save("orbHiper.gif")

"""

Dibujando el potencial efectivo

"""

def UEff(r):
    
    return ((L**2)/(2*n*(r**2)) - (G*M*m)/r)

rEff = np.linspace(0.5e7, 1e8, 25000)
UEffs = UEff(rEff)
EEff = [E for i in range(len(rEff))]

plt.figure(figsize=sizeU)
plt.plot(rEff, UEffs)
plt.plot(rEff, EEff)
plt.xlabel("r(m)")
plt.ylabel("UEff(J)")
#plt.savefig("UEffC.png", dpi=ql)

#Dibujamos el potencial efectivo y la órbita

plt.figure(figsize=(9, 4))
plt.subplot(1, 2, 1)
plt.plot(x, y, label="Solución numérica", linewidth=3)
plt.plot(0, 0, "*r")
plt.plot(xA, yA, "--g", label="Solución analítica")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.legend(loc="lower right")
plt.subplot(1, 2, 2)
plt.plot(rEff, UEffs)
plt.plot(rEff, EEff)
plt.xlabel("r(m)")
plt.ylabel("E(J)")
#plt.savefig("dobleCirc.png", dpi=ql)

"""

Potencial del muelle

"""

#Resolución por odeint

def muelle(z, t):
    
    x, y, vx, vy = z
    
    return [vx, vy, -k*x/m2, -k*y/m2]

xM_0 = 12
yM_0 = 0
vxM_0 = 20
vyM_0 = 2

m_0 = [xM_0, yM_0, vxM_0, vyM_0]
nM = np.linspace(0, 10, 1000)

solM = odeint(muelle, m_0, nM)

xM = solM[:,0]
yM = solM[:,1]
vxM = solM[:,2]
vyM = solM[:,3]

#Gráficas

plt.figure(figsize=sizeM)
plt.plot(xM, yM)
plt.xlabel("x(m)")
plt.ylabel("y(m)")
#plt.savefig("M.png", dpi=ql)

#Animación

figM, axM = plt.subplots()
axM.set_xlim(-max(xM)-2, max(xM)+2)
axM.set_ylim(-max(yM)-0.2, max(yM)+0.2)
lineM, = axM.plot(xM[:1], yM[:1])
aM = 2

def mueAnim(i):
    
    lineM.set_data(xM[:aM*i+1], yM[:aM*i+1])
    
    return lineM

aniM = animation.FuncAnimation(figM, mueAnim, 
                               frames=len(xM)//aM, interval=4)
#aniM.save("aniM.gif")

"""

Añadiendo al muelle un término con r**3

"""

def muelle2(z, t):
    
    x, y, vx, vy = z
    
    return [vx, vy, (-k*x*np.sqrt(x**2 + y**2)**3)/m2, 
            (-k*y*np.sqrt(x**2 + y**2)**3)/m2]

xM2_0 = 10
yM2_0 = 10
vxM2_0 = -100
vyM2_0 = 100

m2_0 = [xM2_0, yM2_0, vxM2_0, vyM2_0]
nM2 = np.linspace(0, 1, 1000)

solM2 = odeint(muelle2, m2_0, nM2)

xM2 = solM2[:,0]
yM2 = solM2[:,1]
vxM2 = solM2[:,2]
vyM2 = solM2[:,3]

#Gráficas

def uM(r):
    
    return (r**5 + 1/r**2)

plt.figure(figsize=sizeM)
plt.plot(xM2, yM2)
plt.xlabel("x(m)")
plt.ylabel("y(m)")
#plt.savefig("M2.png", dpi=ql)

r = np.linspace(0.1, 2.5, 100000)
UEM = uM(r)
plt.figure(figsize=sizeM)
plt.plot(r, UEM)
plt.xlabel("r(m)")
plt.ylabel("$U_eff (J)$")
plt.savefig("UEff.png", dpi=ql)

#Animación

figM2, axM2 = plt.subplots()
axM2.set_xlim(-max(xM2)-1, max(xM2)+1)
axM2.set_ylim(-max(yM2)-1, max(yM2)+1)
lineM2, = axM2.plot(xM2[:1], yM2[:1])
aM2 = 2

def mue2Anim(i):
    
    lineM2.set_data(xM2[:aM2*i+1], yM2[:aM2*i+1])
    
    return lineM2

aniM2 = animation.FuncAnimation(figM2, mue2Anim, 
                                frames=len(xM2)//aM2, interval=5)
#aniM2.save("aniM2.gif")

"""

Resolución del potencial del boletín

"""

def orbR(z, t, L):
    
    x, y, vx, vy = z
    
    return [vx, vy, -8*xR_0**2*x*L**2/(m**2*(x**2 + y**2)**3), 
            -8*xR_0**2*y*L**2/(m**2*(x**2 + y**2)**3)]

xR_0 = 6670000*2
yR_0 = 0
vxR_0 = 0
vyR_0 = np.sqrt(G*M/x_0)
L = m*xR_0*vyR_0

zR_0 = (xR_0, yR_0, vxR_0, vyR_0)
ntR = np.linspace(0, 1000, 25000)

zR = odeint(orbR, zR_0, ntR, args=(L,))

xR = zR[:,0]
yR = zR[:,1]
vxR = zR[:,2]
vyR = zR[:,3]

#Graficamos

plt.figure(figsize=size)
plt.plot(xR, yR)
plt.plot(0, 0, "*r")
plt.xlabel("x(m)")
plt.ylabel("y(m)")
#plt.savefig("orbR.png", dpi=ql)

#Animamos

figR, axR = plt.subplots()
axR.set_xlim(min(xR)-1e6, max(xR)+1e6)
axR.set_ylim(min(yR)-1e5, max(yR)+1e5)
axR.plot(0, 0, "*r")
lineR, = axR.plot(xR[:0], yR[:0])
aR = 50

def orbR(i):
    
    lineR.set_data(xR[:aR*i+1], yR[:aR*i+1])
    
    return lineR

aniR = animation.FuncAnimation(figR, orbR, frames=len(xR)//aR, interval=10)
#aniR.save("Raro.gif")