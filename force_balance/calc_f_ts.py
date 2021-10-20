''' Calculates the strength of a seismogenic
    layer cut by faults with a given dip and
    effective friction using the analytical
    solution given in Turcotte and Schubert
    2002.

    Dependencies: numpy, math'''


# Modules
from numpy import column_stack,savetxt,linspace
from math import sin,cos

n   = 50
Ts  = linspace(0.1*1e3,50*1e3,n)
mu  = 0.05
rho = 2850
g   = 9.81
dip = 30
ddip= dip * (np.pi/180)

# Calculating effective friction
F = rho*g*mu*(Ts**2) / ( sin(2*ddip) - mu*(1+cos(2*ddip)) )

# Saving
# print(out)
out = column_stack((F,Ts))
savetxt('F_Ts_'+str(mu)+'_'+str(dip)+'.out',out)
