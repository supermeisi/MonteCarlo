import os
from os import path
from math import pi
import math
from math import sin, cos
from random import random, uniform, seed
import matplotlib.pyplot as plt

if not path.exists("figs"):
  os.mkdir("figs")

seed(1)

m_mu = 0.1056583755 #muon mass
p0 = 1

N = 100 #number of particles

par = [0.473115591, 0.631038719, 0.906404498, 0.012995717, 0.0041280992, 98.7685322]

#calculating beta value from momentum
def beta(m, p):
  return p/math.sqrt(pow(m, 2) + pow(p, 2))

#calculating refractive index
def rindex(l):
  ll = (l/1000.)**2
  rindex = math.sqrt(par[0]/(1.0-(par[3]/ll)) + par[1]/(1.0-(par[4]/ll)) + par[2]/(1.0-(par[5]/ll)) + 1.0)
  return rindex

#calculating cherenkov angle
def thetac(wlen, m, p):
  thetac = math.acos(1/(rindex(wlen)*beta(m,p)))
  return thetac

#number of photons per track length for specific wavelength
def franck_tamm(wlen, m, p):
  alpha = 1./137. 
  z = 1.
  n = rindex(wlen)
  dndx = 2*math.pi*alpha*pow(z,2)*(1/pow(wlen,2) - 1/(pow(n,2)*pow(beta(m,p),2)*pow(wlen,2)))*1e9
  return dndx

#wavelength distribution
def wavelengh(wlen_min, wlen_max):
  while True:
    wlen = uniform(wlen_min, wlen_max)
    value = uniform(0,1000)
    if value < franck_tamm(wlen, m_mu, p0):
      break
  return wlen

#muon generator
def muongen(x0, y0, dx, dy, phi_min, phi_max):
  x0 = x0 - dx/2 + dx*random()
  y0 = y0 - dy/2 + dy*random()
  theta = math.acos((1-random())**(1/3))
  phi = uniform(phi_min,phi_max)
  u = 0
  x = 0
  while True:
    x = 500 + 9500*random()
    u = random()*1e-4
    if u < 1/x**2.5:
      break
  mom = math.sqrt(x**2 - m_mu**2)  
  return x0, y0, theta, phi, mom

arr_x = []
arr_y = []
arr_theta = []
arr_phi = []
arr_mom = []

arr_thetac = []
arr_wlen = []

#simulating particle tracks
for i in range(N):
  x0, y0, theta0, phi0, mom0 = muongen(0, 0, 10, 10, -3, 3)

  print(x0, y0, theta0, phi0, mom0)

  arr_x.append(x0)
  arr_y.append(y0)
  arr_theta.append(theta0)
  arr_phi.append(phi0)
  arr_mom.append(mom0)

  for j in range(1000):
    wlen = wavelengh(300,800)
    arr_wlen.append(wlen)
    thetac0 = thetac(wlen, m_mu, mom0)
    arr_thetac.append(thetac0)
    phic0 = 2*pi*random()
    #direction of Cherenkov photons in lab system
    pxc = cos(phi0)*(cos(phic0)*cos(theta0)*sin(thetac0)+sin(theta0)*cos(thetac0))-sin(phi0)*sin(phic0)*sin(thetac0)
    pyc = sin(phi0)*(cos(phic0)*cos(theta0)*sin(thetac0)+sin(theta0)*cos(thetac0))+cos(phi0)*sin(phic0)*sin(thetac0)
    pzc = cos(theta0)*cos(thetac0)-cos(phic0)*sin(theta0)*sin(thetac0)
    print(pxc, pyc, pzc)

n_bins = 10
plt.hist(arr_x, bins=n_bins)
plt.savefig("figs/x.pdf")
plt.close()
plt.hist(arr_y, bins=n_bins)
plt.savefig("figs/y.pdf")
plt.close()
plt.hist(arr_theta, bins=n_bins)
plt.savefig("figs/theta.pdf")
plt.close()
plt.hist(arr_phi, bins=n_bins)
plt.savefig("figs/phi.pdf")
plt.close()
plt.hist(arr_mom, bins=n_bins)
plt.savefig("figs/mom.pdf")
plt.close()
plt.hist(arr_thetac, bins=n_bins)
plt.savefig("figs/thetac.pdf")
plt.close()
plt.hist(arr_wlen, bins=n_bins)
plt.savefig("figs/wlen.pdf")
plt.close()