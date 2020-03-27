import numpy as np

ms = 50
g = 9.81
rho = 1.091
r = 0.5
A = np.pi*r**2
ve = 325
cd = 0.15
mp0 = 100
delta_t = 0.1
total_time = 50
nt = int(total_time/delta_t) + 1
mpt = int(5.0/delta_t)

t = np.linspace(0,total_time,nt)
mpdot = np.concatenate([20.0*np.ones(mpt),np.zeros(nt-mpt)])
mp = np.zeros(nt)
mp[:mpt] = mp0 - mpdot[:mpt]*t[:mpt] 

v = np.zeros(nt)
for n in range(nt-1):
    v[n+1] = v[n] + delta_t*(-(ms+mp[n])*g+ \
     mpdot[n]*ve-0.5*rho*v[n]*abs(v[n])*A*cd)/(ms+mp[n])
    
h = np.zeros(nt)
for n in range(nt-1):
    h[n+1] = h[n] + delta_t*v[n]

