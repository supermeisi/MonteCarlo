import math
from random import random, uniform

def muongen(x0, y0, dx, dy, phi_min, phi_max):
  x0 = x0 - dx/2 + dx*random()
  y0 = y0 - dy/2 + dy*random()
  theta = math.acos(math.sqrt(1-3*u))/1.746549412046033
  phi = uniform(phi_min,phi_max)
  return x0, y0, theta, phi

