#!/usr/bin/python

import os
from os import path
from math import pi
import math
from math import sin, cos, sqrt, asin, acos, tan, atan, radians
from random import random, uniform, seed, gauss
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import threading
import scipy.integrate as integrate

if not path.exists("figs"):
  os.mkdir("figs")

seed(1)

multithread = False

m_mu = 0.1056583755 #muon mass
p0 = 1

N = 1 #number of particles

n_threads = 1 #number of threads

par = [0.473115591, 0.631038719, 0.906404498, 0.012995717, 0.0041280992, 98.7685322]

sx = 0
sy = 0
stheta = 0
sphi = 0
R = 0.8

n_roms = 288

wlen_min = 300
wlen_max = 300

#mathematical functions
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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
def wavelengh(wlen_min, wlen_max, m, p):
  while True:
    wlen = uniform(wlen_min, wlen_max)
    value = uniform(0,1000)
    if value < franck_tamm(wlen, m, p):
      break
  return wlen

#muon generator
def muongen(x0, y0, dx, dy, phi_min, phi_max, theta_min, theta_max):
  x0 = x0 - dx/2 + dx*random()
  y0 = y0 - dy/2 + dy*random()
  theta = math.acos((1-(theta_min/radians(90)+random()*(theta_max/radians(90)-theta_min/radians(90))))**(1/3))
  phi = uniform(phi_min,phi_max)
  u = 0
  x = 0
  while True:
    x = 0.7 + 9.3*random()
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
arr_phi_prime = []
arr_thetac_reco = []

def simulation():
  #simulating particle tracks
  for i in range(int(N/n_threads)):
    theta0 = radians(5)
    phi0 = radians(0)
    x0 = sin(theta0)*cos(phi0)
    y0 = sin(theta0)*sin(phi0)
    mom0 = 1
    
    arr_x.append(x0)
    arr_y.append(y0)
    arr_theta.append(theta0)
    arr_phi.append(phi0)

    p2 = [cos(phi0),sin(phi0)]

    nphotons, error = integrate.quad(lambda x: franck_tamm(x,m_mu,mom0), wlen_min, wlen_max)
    
    n_hits = 0

    for j in range(n_roms):
      wlen = wavelengh(wlen_min, wlen_max, m_mu, mom0)
      thetac0 = thetac(wlen, m_mu, mom0)
      arr_thetac.append(thetac0)
      angle_fel = radians(360)/n_roms*j
      fel_x = R*cos(angle_fel)
      fel_y = R*sin(angle_fel)
      dx = fel_x - x0
      dy = fel_y - y0
      dr = [dx, dy]
      phirel = angle(p2,dr)
      norm = [cos(angle_fel),sin(angle_fel)]
      alpha = angle(norm,dr)
      A = sin(theta0)*cos(phirel)
      B = A**2 + (cos(theta0))**2
      phi = acos(A*cos(thetac0)/B+sqrt(((cos(theta0))**2-(cos(thetac0))**2)/B+(A*cos(thetac0)/B)**2))
      phi_prime = atan(tan(phi)/cos(alpha))
      if(phi_prime < 21*pi/180 or phi_prime > 41*pi/180):
        continue      
      angle_total = math.asin(1/rindex(wlen))
      theta_photon = pi/2-phi_prime
      if theta_photon < angle_total:
        continue
      arr_phi_prime.append(phi_prime)
      dx = fel_x - gauss(x0,sx)
      dy = fel_y - gauss(y0,sy)
      dr = [dx,dy]
      phirel = angle(p2,dr)
      alpha = angle(norm,dr)
      phi = atan(tan(phi)*cos(alpha))
      theta0 = gauss(theta0,stheta)
      phirel = angle(p2,dr)
      thetac0 = acos(sin(theta0)*cos(phirel)*cos(phi)+cos(theta0)*sin(phi))
      print(thetac0)
      arr_thetac_reco.append(thetac0)
      n_hits += 1

  print(n_hits)

if multithread:
  threads = []

  for index in range(n_threads):
    x = threading.Thread(target=simulation)
    threads.append(x)
    x.start()


  for index, thread in enumerate(threads):
    thread.join()
else:
  simulation()

n_bins = int(sqrt(N))
plt.hist(arr_x, bins=n_bins)
plt.xlabel('x [mm]')
plt.ylabel('Entries')
plt.savefig("figs/x.pdf")
plt.savefig("figs/x.png")
plt.close()
plt.hist(arr_y, bins=n_bins)
plt.xlabel('y [mm]')
plt.ylabel('Entries')
plt.savefig("figs/y.pdf")
plt.savefig("figs/y.png")
plt.close()
plt.hist(arr_theta, bins=n_bins)
plt.xlabel('$\\theta_p$ [rad]')
plt.ylabel('Entries')
plt.savefig("figs/theta.pdf")
plt.savefig("figs/theta.png")
plt.close()
plt.hist(arr_phi, bins=n_bins)
plt.xlabel('$\phi_p$ [rad]')
plt.ylabel('Entries')
plt.savefig("figs/phi.pdf")
plt.savefig("figs/phi.png")
plt.close()
plt.hist(arr_mom, bins=n_bins)
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('Entries')
plt.savefig("figs/mom.pdf")
plt.savefig("figs/mom.png")
plt.close()
plt.hist(arr_thetac, bins=n_bins)
plt.xlabel('$\\theta_c0$ [rad]')
plt.ylabel('Entries')
plt.savefig("figs/thetac.pdf")
plt.savefig("figs/thetac.png")
plt.close()
plt.hist(arr_wlen, bins=n_bins)
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('Entries')
plt.savefig("figs/wlen.pdf")
plt.savefig("figs/wlen.png")
plt.close()
plt.hist(arr_phi_prime, bins=n_bins)
plt.xlabel('$\\varphi_c$ [rad]')
plt.ylabel('Entries')
plt.savefig("figs/phi_prime.pdf")
plt.savefig("figs/phi_prime.png")
plt.close()
plt.hist(arr_thetac_reco, bins=n_bins)
mu, std = norm.fit(arr_thetac_reco)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
ax = plt.gca()
std = round(std*1e3,4)
plt.text(0.7,0.8,str('$\sigma$ = '+str(std)+' mrad'),transform=ax.transAxes)
plt.xlabel('$\\theta_{c,reco}$ [rad]')
plt.ylabel('Entries')
plt.savefig("figs/thetac_reco.pdf")
plt.savefig("figs/thetac_reco.png")
plt.close()