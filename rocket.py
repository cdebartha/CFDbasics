import numpy as np
import matplotlib.pyplot as pp

def euler_step(u,f,dt,*args):
    u_new = u + dt * f(u,*args)
    return u_new

def rk2_step(u,f,dt,*args):
    u_star = u + 0.5 * dt * f(u,*args)
    u_new = u + dt * f(u_star,*args)
    return u_new

def rhs_rocket(u, ms, mp, mpdot, ve, rho, A, cd, g):
    v, h = u
    rhs = np.array([(-(ms+mp)*g + mpdot*ve - 0.5*rho*v*abs(v)*A*cd)/(ms+mp) , v])
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
T = 10
N = int(T/dt) + 1

t = np.linspace(0,T,N)

mpdot = np.zeros(N)
mask = np.where(t < 5.0)
mpdot[mask] = 20.0

mp = np.zeros(N)
mp[mask] = mp0 - mpdot[mask]*t[mask] 

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
u2 = np.empty((N,2))
u[0] = np.array([v0,h0])
u2[0] = np.array([v0,h0])

for n in range(N-1):
    u[n+1] = euler_step(u[n], rhs_rocket, dt, ms, mp[n], mpdot[n], ve, rho, A, cd, g)
    u2[n+1] = rk2_step(u2[n], rhs_rocket, dt, ms, mp[n], mpdot[n], ve, rho, A, cd, g)

v = u[:,0]
h = u[:,1]

v2 = u2[:,0]
h2 = u2[:,1]

pp.figure(figsize=(9.0, 4.0))
pp.title('Height of the rocket vs time')
pp.xlabel('t')
pp.ylabel('h')
pp.grid()
pp.plot(t, h, color='b', linestyle='-', linewidth=2);
pp.plot(t, h2, color='r', linestyle='-', linewidth=2);

pp.figure(figsize=(9.0, 4.0))
pp.title('Velocity of the rocket vs time')
pp.xlabel('t')
pp.ylabel('v')
pp.grid()
pp.plot(t, v, color='b', linestyle='-');
pp.plot(t, v2, color='r', linestyle='-');