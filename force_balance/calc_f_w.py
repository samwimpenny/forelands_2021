''' Calculates the force needed
    to break an asperity on a fault
    for a given effective coefficient
    of friction, dip, area and
    centroid depth.

    Dependencies: numpy '''

# Modules
import numpy as np
import matplotlib.pyplot as plt

# Variables:
dip   = 30
theta = dip * (np.pi/180)
rho   = 2850
g     = 9.81
W     = 5e3
mu    = 0.4

# Defining calculation function
def calc_F(T1,mu,theta,rho,g,W):
    T2 = T1 + ( W*np.sin(np.pi*0.5-theta) )
    F  = (rho*g*mu)*(T2**2 - T1**2) \
        / (np.sin(2*theta) - mu*(1+np.cos(2*theta)))
    return F

# Creating grids to contour:
T1 = np.linspace(0,40e4,num=20)
F = calc_F(T1,mu,theta,rho,g,W)
C = T1 + 0.5*W*np.sin(np.pi*0.5-theta)

out = np.column_stack((F,C))
np.savetxt('F_W_'+str(mu)+'_'+str(W/1e3)+'_'+str(dip)+'.out',out)
