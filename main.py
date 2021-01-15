import os
from os import path
import math
from random import random, uniform, seed
import matplotlib.pyplot as plt

if not path.exists("figs"):
  os.mkdir("figs")

seed(1)

m_mu = 105.658 #muon mass

#muon generator
def muongen(x0, y0, dx, dy, phi_min, phi_max):
  x0 = x0 - dx/2 + dx*random()
  y0 = y0 - dy/2 + dy*random()
  theta = math.acos((1-random())**(1/3))
  phi = uniform(phi_min,phi_max)
  u = 0
  x = 0
  while(1):
    x = 500 + 9500*random()
    u = random()*1e-4
    if u < 1/x**2.5:
      break
  mom = math.sqrt(x**2 - m_mu**2)  
  return x0, y0, theta, phi, mom

x = []
y = []
theta = []
phi = []
mom = []

#simulating particle tracks
for i in range(100):
  x0, y0, theta0, phi0, mom0 = muongen(0, 0, 10, 10, -3, 3)

  print(x0, y0, theta0, phi0, mom0)

  x.append(x0)
  y.append(y0)
  theta.append(theta0)
  phi.append(phi0)
  mom.append(mom0)

plt.hist(x)
plt.savefig("figs/x.pdf")
plt.close()
plt.hist(y)
plt.savefig("figs/y.pdf")
plt.close()
plt.hist(theta)
plt.savefig("figs/theta.pdf")
plt.close()
plt.hist(phi)
plt.savefig("figs/phi.pdf")
plt.close()
plt.hist(mom)
plt.savefig("figs/mom.pdf")
plt.close()