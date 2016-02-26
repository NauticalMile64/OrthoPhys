from sympy import *
from sympy.functions import *

class CoordinateSystem:

	def __init__(self,crds,crdDefs):

		self.crdDefs = crdDefs
		self.n = n = len(crdDefs)

		#Compute the derivatives of each of the co-ordinate definitions with respect to each co-ordinate
		dcrds = [[simplify(diff(crdDef,crd)) for crdDef in crdDefs] for crd in crds]

		#Compute the scale factors
		h = zeros(n,1)
		for i in range(n):
			for j in range(n):
				h[i] += dcrds[i][j]**2
			h[i] = simplify(root(h[i],2))
		self.scaleFactors = h

		# Volume element
		dVc = 1
		for i in range(n):
			dVc *= h[i]
		self.dVc = dVc

		#Define the unit vectors for the co-ordinate system: http://mathworld.wolfram.com/UnitVector.html
		uvecs = zeros(n)
		for i in range(n):
			for j in range(n):
				uvecs[i,j] = simplify(dcrds[i][j]/h[i])
		self.uvecs = uvecs

		#Metric Tensors
		g = zeros(n)
		for i in range(n):
			for j in range(n):
				for k in range(n):
					g[i,j] += dcrds[i][k]*dcrds[j][k]
				g[i,j] = simplify(g[i,j])

		gi = g**-1
		self.metricT = g

		#Christoffel Symbol of second kind
		Tau2 = [zeros(n) for i in range(n)]
		for i in range(n):
			for j in range(n):
				for k in range(n):
					expr = simplify((diff(g[i,k],crds[j])
					+	diff(g[j,k],crds[i])
					-	diff(g[i,j],crds[k]))/2)
					for m in range(n):
						Tau2[m][i,j] += gi[k,m]*expr
		self.christoffel2 = simplify(Tau2)
