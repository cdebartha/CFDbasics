import numpy as np
import matplotlib.pyplot as pp

def euler_step(u,f,dt,*args):
    u_new = u + dt * f(u,*args)
    return u_new

def rhs_rocket(u, ms, mp, mpdot, ve, rho, A, cd, g):
    v, h = u
    rhs = np.array([(-(ms+mp)*g+ \
     mpdot*ve-0.5*rho*v*abs(v)*A*cd)/(ms+mp), v])
    return rhs

ms = 50
g = 9.81
rho = 1.091
r = 0.5
A = np.pi*r**2
ve = 325
cd = 0.15
mp0 = 100

dt = 0.1
T = 50
N = int(T/dt) + 1
M = int(5.0/dt)

t = np.linspace(0,T,N)
mpdot = np.concatenate([20.0*np.ones(M),np.zeros(N-M)])
mp = np.zeros(N)
mp[:M] = mp0 - mpdot[:M]*t[:M] 

#v = np.zeros(N)
#for n in range(N-1):
#    v[n+1] = v[n] + dt*(-(ms+mp[n])*g+ \
#     mpdot[n]*ve-0.5*rho*v[n]*abs(v[n])*A*cd)/(ms+mp[n])
    
#h = np.zeros(N)
#for n in range(N-1):
#    h[n+1] = h[n] + dt*v[n]

v0 = 0.0
h0 = 0.0

u = np.empty((N,2))
u[0] = np.array([v0,h0])

for n in range(N-1):
    u[n+1] = euler_step(u[n], rhs_rocket, dt, ms, mp[n], mpdot[n], ve, rho, A, cd, g)

v = u[:,0]
h = u[:,1]

pp.figure(figsize=(9.0, 4.0))
pp.title('Height of the rocket vs time')
pp.xlabel('t')
pp.ylabel('h')
pp.grid()
pp.plot(t, h, color='C0', linestyle='-', linewidth=2);

pp.figure(figsize=(9.0, 4.0))
pp.title('Velocity of the rocket vs time')
pp.xlabel('t')
pp.ylabel('v')
pp.grid()
pp.plot(t, v, color='C0', linestyle='-', linewidth=2);