from scipy.integrate import quad
from scipy.interpolate import interp1d
from numpy import inf, loadtxt, linspace
from math import pi, fabs, log10
from pashas_mod import Planck_fun
from os import remove, system
from time import sleep
import matplotlib.pyplot as plt

m_re = 1.33
m_im = -0.02

R = 2*10e11
T = 1.5*10e4
M = 2*10e33
#10e-7 10e-4cm
a = 10e-7
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
	ls = linspace(a/40.5, 1e-1, 50)
	print(ls)
	try:
		remove('mie.out')
	except:
		print('File mie.out not found... Continue')
	f = open('mie.dat', 'w')

	#WRITE INPUT FILE TO FORTRAN MODULE
	m_re_str ="%1.2f"%m_re
	m_im_str ="%1.2f"%m_im 
	f.write(m_re_str+' '*(16-len(m_re_str))+m_im_str+' '*(17-len(m_im_str)) +
										 'ri: refractive index\n')
	f.write("50                               nx: number of sizes\n")
	i = 0
	for lamb in ls:
		x = 2*pi*a/lamb
		if i == 0:
			l_str = '%1.13f'%x
			f.write(l_str+' '*(34-len(l_str))+'x: size 1\n')
			i = 1
		else:
			l_str = '%1.13f'%x
			f.write(l_str+' '*(34-len(l_str))+'x: size ...\n')
	f.close()
	
	#RUN FORTRAN MODULE	
	system('./a.out')
	sleep(0.01)
	#READ mie.out
	lambd, Qpr = loadtxt('mie.out', skiprows = 5, usecols = [0,5], unpack = True,  comments = '-')
	#print(lambd, Qpr)
	Qpextrapolate = interp1d(lambd, Qpr, fill_value="extrapolate")
	def Qpr(l):
		if Qpextrapolate(l)>0:
			return Qpextrapolate(l)
		else:
			return 0
	#x = []
	#y = []
	#for l in linspace(0,100, 1000):
	#	x.append(l)
	#	y.append(Qpr(l))
	#plt.plot(x, y)
	#plt.show()
	f1 = lambda x: pi*Planck_fun(10e8*x, T)
	f2 = lambda x: Qpr(x)*f1(x)
	#print(Qpr(23072))
	#print('asdf')Q
	#print(quad(f2, 0, 23072)[0])
	#print(quad(f2, 0, 23072)[0]/quad(f1, 0, 23072)[0])
	return quad(f2, 0, inf)[0]/quad(f1, 0, inf)[0]

def main():
	x, y = [], []
	i = 0
	for a in linspace(10e-4, 1e-7, 25):
		print(i)
		x.append(log10(a))
		y.append(Q_mean(T, m_re, a))
		i+=1
		#y.append(beta(R, T, M, Q_mean(T, m_c, a), a, p))
	print(x)
	print(y)
	plt.plot(x, y, 'o')
	plt.plot(x, y)
	plt.show()

if __name__ == '__main__':
	main()
#x, y = [], []

#for l in linspace(10e-8,3*10e-7, 200):
#	x.append(l)
#	y.append(Planck_fun(10e8*l, T))
#plt.plot(x, y)
#plt.show()

#for T in linspace(0.1, 1000, 10):
#	print(T)
#	x.append(T)
#	y.append(beta(R, T, M, Q_mean(T, m_i, a), a, p))
#plt.plot(x, y)
#plt.show()

#print(Q_mean(T, m_c, a))
	

#print(Qpr(123411))
#plt.plot(x, y)
#plt.show()