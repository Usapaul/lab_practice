from scipy.integrate import quad
from scipy.interpolate import interp1d
from pickle import dump
from numpy import inf, loadtxt, linspace, genfromtxt, delete
from math import pi, fabs, log10
from pashas_mod import Planck_fun
from os import remove, system
from time import sleep
import matplotlib.pyplot as plt
from matplotlib import style

style.use('bmh')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

T = 6580                       
M = 2.5616456E+33                                 
R = 8.3641537E+10
p = 2.5 


def Q_mean(T):
	'''
	Plank factor of light pressure

	x - wavelength
	a - radius (dust)
	m - refractive index (dust)
	'''
	A = loadtxt('a.dat',  unpack = True)
	B = []
	for j in range(0, len(A)):
		a = A[j]
		print('Work with a = %1.3f'%a)
		lambd, Q = loadtxt('bigInput.dat',  usecols = (0, j+1),  unpack = True)
		indexes = [i for i, x in enumerate(list(Q)) if x == 0.]
		Q = delete(Q, indexes)
		lambd = delete(lambd, indexes)
		Qpextrapolate = interp1d(lambd, Q, fill_value="extrapolate")
		def Qextr(l):
			if Qpextrapolate(l)>0:
				return Qpextrapolate(l)
			else:
				return 0
		#X = linspace(lambd[0], lambd[-1], 10000)
		#q = [Qextr(x) for x in X]
		#plt.plot(X, q, label = '$a = {}  \mu m$'.format(a))
		#plt.show()
		f1 = lambda x: pi*Planck_fun(10e3*x, T)
		f2 = lambda x: Qextr(x)*f1(x)
		def Qpr_mean(T):
			return quad(f2, lambd[0], lambd[-1])[0]/quad(f1, lambd[0], lambd[-1])[0]
		Qint = Qpr_mean(T)
		print(Qint)
		B.append(beta(R, T, M, Qint, a, p))	
	dump([A,B], open('AB.obj', 'wb'))
	plt.semilogx(A, B)
	plt.xlim(0, 1)
	plt.ylabel('beta')
	plt.xlabel('$a, \mu m $')
	plt.savefig('F5V_c.png')
	plt.show()

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
	return 2.12*1e-8*R**2*T**4*Q/(M*a*p)

def main():
	x, y = [], []
	i = 0
	Q_mean(T)

if __name__ == '__main__':
	main()
