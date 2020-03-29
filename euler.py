import numpy as np

def euler_step(u,f,dt,*args):
    u_new = u + dt * f(u,*args)
    return u_new

def rk2_step(u,f,dt,*args):
    u_star = u + 0.5 * dt * f(u,*args)
    u_new = u + dt * f(u_star,*args)
    return u_new