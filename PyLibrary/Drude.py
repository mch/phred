#!/usr/bin/python

# Plots the propagation and attenuation constants for a wave in a
# Drude media over a range of frequencies

from Numeric import *
import biggles

mu = 1.256637061435917e-06
eps_0 = 8.854185336732027e-12

# Plasma Frequency
w_p = 2e15 * 2 * pi

# Collision Frequency
v = 5.7e13 # * 2 * pi

# Frequency of interest
#w = arrayrange(300, 750, 20) * 1e12 * 2 * pi;
w = (3e8 / arrayrange(0.4e-6, 0.8e-6, 0.05e-6)) * 2 * pi;
#print w

k_0 = w * sqrt(mu * eps_0)

# Relative dielectric constant
eps_r = 1 + w_p**2 / (- w**2 + 1j*w*v)

eps = eps_r
#print eps_r

# Attenuation constant
alpha = k_0 * 1j * sqrt(abs(eps.real) / 2) \
        * (sqrt(1 + (eps.imag / eps.real)**2) - 1)**0.5

# Propagation constant
beta = k_0 * 1j * sqrt(abs(eps.real) / 2) \
       * (sqrt(1 + (eps.imag / eps.real)**2) + 1)**0.5

p = biggles.FramedPlot()
p.xlabel = "Frequency (THz)"
p.ylabel = "Relative Permativity"
#p.ylabel = "Attenuation"
#p.add( biggles.Curve(w / (1e12 * 2 * pi), alpha, color="red") )
p.add( biggles.Curve(w / (2 * pi), eps_r.real, color="red") )
p.add( biggles.Curve(w / (2 * pi), eps_r.imag, color="green") )
p.show()  
