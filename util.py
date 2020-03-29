import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify as ld

sp.init_printing()

x = sp.symbols('x')
phi = sp.cos(x)**2*sp.sin(x)**3/(4*x**5*sp.exp(x))
phiprime = phi.diff(x)

ans = ld((x), phiprime)

print(ans(2.2))
