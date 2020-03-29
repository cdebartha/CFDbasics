import numpy as np
import matplotlib.pyplot as pp
import matplotlib.animation as am

def F(Vmax, rho_max, rho):
    F = Vmax*rho*(1-rho/rho_max)
    return F

def Vel(Vmax, rho_max, rho):
    V = Vmax*(1-rho/rho_max)
    return V

def convection(u0, dx, dt, nt, Vmax, rho_max):
    u_hist = [u0.copy()]
    u = u0.copy()
    for n in range(nt):
        un = u.copy()
        #Update boundary condition
        u[0] = 20
        #Update other points
        u[1:] = un[1:] - dt / dx * (F(Vmax, rho_max, u[1:]) - F(Vmax, rho_max, u[:-1]))
        u_hist.append(u.copy())
    return u_hist

def update_plot(n, u_hist, fig, line):
    fig.suptitle('Time step {:0>2}'.format(n))
    line.set_ydata(u_hist[n])

#parameters
Vmax = 136.0
L = 11.0
rho_max = 250.0
nx = 51
dt = 0.001
dx = L/(nx-1)
nt = 101

#discretize grid
x = np.linspace(0.0,L,nx)

#set initial conditions
rho0 = np.ones(nx)*20
rho0[10:20] = 50

#compute density history of numerical solution
rho_hist = convection(rho0, dx, dt, nt, Vmax, rho_max)

#compute velocity history
v0 = Vel(Vmax, rho_max, rho0)
v_hist = [v0.copy()]
for n in range(1,nt):
    v_hist.append(Vel(Vmax, rho_max, rho_hist[n]))

fig_rho = pp.figure(figsize=(6.0, 4.0))
pp.xlabel('x')
pp.ylabel('rho')
pp.grid()
line_rho = pp.plot(x, rho0, color='b', linestyle='-', linewidth=2)[0]
pp.xlim(0.0, L)
pp.ylim(0.0, 60)
fig_rho.tight_layout()

anim_rho = am.FuncAnimation(fig_rho, update_plot, frames=nt, fargs=(rho_hist, fig_rho, line_rho), interval=100)

fig_v = pp.figure(figsize=(6.0, 4.0))
pp.xlabel('x')
pp.ylabel('vel')
pp.grid()
line_v = pp.plot(x, v0, color='b', linestyle='-', linewidth=2)[0]
pp.xlim(0.0, L)
pp.ylim(100, Vmax)
fig_v.tight_layout()

anim_v = am.FuncAnimation(fig_v, update_plot, frames=nt, fargs=(v_hist, fig_v, line_v), interval=100)

