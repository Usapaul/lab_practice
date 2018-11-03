from scipy.integrate import quad
from numpy import inf, loadtxt, linspace
from math import pi, fabs
from pashas_mod import Planck_fun
import matplotlib.pyplot as plt

m_i = 1.31 #for ice particles
m_c = 2.5 #for carbon
m_si =3.3 #for silicon

R = 7*10e10
T = 5.0*10e3
M = 2*10e33
a = 10e-4
p = 2.5 

def beta(R, T, M, Q, a, p):
	"""
	Gravity force to light pressure force fractoin.

	R - radius  (star)
	T - temperature  (star)
	M - mass  (star)
	Q - Plank factor of light pressure
	a - radius (dust)
	p - density (dust)
	"""
	return 2.12*10e-8*R**2*T**4*Q/(M*a*p)

def Q_mean(T, m, a):
	'''
	Plank factor of light pressure

	x - wavelength
	a - radius (dust)
	m - refractive index (dust)
	'''
	f1 = lambda x: pi*Planck_fun(10e8*x, T)
	f2 = lambda x: (8/3)*(2*pi*a/x)**4*abs(((m**2 - 1)/(m**2 + 2))**2)*f1(x)
	return quad(f2, 0, inf)[0]/quad(f1, 0, inf)[0]

x, y = [], []

#for l in linspace(10e-8,3*10e-7, 200):
#	x.append(l)
#	y.append(Planck_fun(10e8*l, T))
#plt.plot(x, y)
#plt.show()

for T in linspace(1000, 50000, 10):
	print(T)
	x.append(T)
	y.append(Q_mean(T, m_si, a))

plt.plot(x, y)
plt.show()


#print(beta(R, T, M, Q_mean(T, m_i, a), a, p))