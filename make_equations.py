from sympy import *
from sympy.functions import *
from coordsys import CoordinateSystem

r = symbols('r', real = True, positive = True)
theta, z = symbols('theta z', real = True)

a,tau = symbols('a tau', real = True, positive = True)
sigma = symbols('sigma', real = True)

x,y,z = symbols('x y z', real = True)

#Cartesian
#crds = [x,y,z]
#crdDefs = Matrix([x,y,z])

# Cylindrical
#crds = [r,theta,z]
#crdDefs = Matrix([r*cos(theta),r*sin(theta),z])

#Spherical
phi = symbols('phi', real = True, positive = True)
crds = [r,theta,phi]
crdDefs = Matrix([r*cos(theta)*sin(phi),r*sin(theta)*sin(phi),r*cos(phi)])

#Toroidal
#phi = symbols('phi', real = True)
#crds = [sigma,tau,phi]
#crdDefs = Matrix([cos(phi)*sinh(tau), sin(phi)*sinh(tau), sin(sigma)])*a/(cosh(tau)-cos(sigma))

# #############################

cds = CoordinateSystem(crds,crdDefs)

import matplotlib.pyplot as plt

def latex_print(text,vpos,hpos = 0.5):
	plt.text(hpos,vpos,'$%s$' %text,size = 'x-large')

print('####### Volume Element Coeff.')
latex_print(latex(cds.dVc),0.8)

print('####### Unit Vectors')
for i in range(3):
  for j in range(3):
    latex_print(latex(cds.uvecs[i,j]),0.6-i*0.1,0.2+0.2*j)

print('####### Metric Tensor')
#latex_print(latex(cds.metricT),0.4)

print('####### Christoffel Symbols')
#latex_print(latex(cds.christoffel2),0.2)

rho, c, k, T, t = symbols('rho c k T t', real=True, positive=True)

consE = rho*c*diff(T(t),t) #- grad(-k*grad(T))
latex_print(latex(consE),0.2)

plt.show()
