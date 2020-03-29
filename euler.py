import numpy as np

def euler_step(u,f,dt,*args):
    u_new = u + dt * f(u,*args)
    return u_new