import sympy as spy
from sympy import init_printing
from sympy import *
from sympy.functions import *

init_printing()

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
#phi = symbols('phi', real = True, positive = True)
#crds = [r,theta,phi]
#crdDefs = Matrix([r*cos(theta)*sin(phi),r*sin(theta)*sin(phi),r*cos(phi)])

#Toroidal
phi = symbols('phi', real = True)
crds = [sigma,tau,phi]
crdDefs = Matrix([cos(phi)*sinh(tau), sin(phi)*sinh(tau), sin(sigma)])*a/(cosh(tau)-cos(sigma))

# #############################

Rn = len(crds)
print('####### Co-ordinates')
print(latex(crdDefs))

#Compute the derivatives of each of the co-ordinate definitions with respect to each co-ordinate
dcrds = [[simplify(diff(crdDef,crd)) for crdDef in crdDefs] for crd in crds]

#Compute the scale factors
h = zeros(Rn,1)
for i in range(Rn):
	for j in range(Rn):
		h[i] += dcrds[i][j]**2
	h[i] = simplify(root(h[i],2))

print('####### Scale Factors')
print(latex(h))

# Volume element
dVc = 1
for i in range(Rn):
	dVc *= h[i]

print('####### Volume Element Coeff.')
print(latex(dVc))

#Define the unit vectors for the co-ordinate system: http://mathworld.wolfram.com/UnitVector.html
uvecs = zeros(Rn)
for i in range(Rn):
	for j in range(Rn):
		uvecs[i,j] = simplify(dcrds[i][j]/h[i])

print('####### Unit Vectors')
print(latex(uvecs))

#Metric Tensors
g = zeros(Rn)
for i in range(Rn):
	for j in range(Rn):
		for k in range(Rn):
			g[i,j] += dcrds[i][k]*dcrds[j][k]
		g[i,j] = simplify(g[i,j])

print('####### Metric Tensor')
print(latex(g))
gi = g**-1

#Christoffel Symbol of second kind
Tau2 = [zeros(Rn) for i in range(Rn)]
for i in range(Rn):
	for j in range(Rn):
		for k in range(Rn):
			expr = simplify((diff(g[i,k],crds[j])
				+	diff(g[j,k],crds[i])
				-	diff(g[i,j],crds[k]))/2)
			for m in range(Rn):
				Tau2[m][i,j] += gi[k,m]*expr
				#Tau2[m][i,j] = simplify(Tau2[m][i,j])
Tau2 = simplify(Tau2)

print('####### Christoffel Symbols')
print(latex(Tau2))

# ########################
rho, c, k, T, t = symbols('rho c k T t', real=True, positive=True)

#consE = rho*c*diff(T,t) - grad(-k*grad(T))
