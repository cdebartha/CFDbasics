import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify as ld

sp.init_printing()

x, nu, t = sp.symbols('x nu t')
phi = (sp.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sp.exp(-(x - 4 * t - 2 * sp.pi)**2 / (4 * nu * (t + 1))))
phiprime = phi.diff(x)

u = -2 * nu * (phiprime / phi) + 4
u_lamb = ld((t, x, nu), u)

nx = 101
L = 2.0 * np.pi
dx = L / (nx - 1)
nu = 0.07
nt = 100
sigma = 0.1
dt = sigma * dx**2 / nu

x = np.linspace(0.0, L, nx)

u0 = np.array([u_lamb(0.0, xi, nu) for xi in x])



