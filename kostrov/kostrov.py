''' Script for performing a summation of a set of
    earthquake moment tensors to give the principle
	axes of the strain rate tensor.

    Dependencies: seismic.py, numpy, sys '''

# Uses convention x = east, y = north and z = down.

# Modules
import numpy as np
import seismic as seis
import sys

# Loading information
fdat = np.loadtxt(sys.argv[1]) #[strike, dip, rake, M0 (Nm)]
V = float(sys.argv[2]) #[m^3, kostrov volume]
t = float(sys.argv[3]) #[yr, time of catalogue]
mu = 30e9 #[Pa], Rigidity modulus

# Creating store array
mt_sum=np.zeros((3,3))

# Summing moment tensors
for i in range(0,len(fdat[:,0])):

    # Convert strike/dip/rake to unit MT & multiply by moment
    mt = seis.sdr2smt(fdat[i,0],fdat[i,1],fdat[i,2]) * fdat[i,3]
    mt_sum = mt_sum + mt

# Converting fro up-south-east to down-north-east
# [muu mus mue] -> [mdd mdn mde]
# [mus mss mse] -> [mnd mnn mne]
# [mue mse mee] -> [mde mne mee]
mt_sum = np.array([ [mt_sum[0,0],mt_sum[0,1],-1*mt_sum[0,2]],
                    [mt_sum[1,0],mt_sum[1,1],-1*mt_sum[1,2]],
                    [-1*mt_sum[2,0],-1*mt_sum[2,1],mt_sum[2,2]] ])

# Calculating cumulative moment
# M0 = (sum(Mij)^2)^0.5
Msum = np.linalg.norm(mt_sum,'fro') / np.sqrt(2)
print('Unit Moment tensor: \n', mt_sum/Msum)
print('Cumulative Moment:', Msum)

# Converting to strain rate tensor averaged over volume and obs time (1/yr)
e_bar = mt_sum / (2*mu*V*t)

# Creating horizontal strain rate tensor
# [Exx,Exy]
# [Eyx,Eyy]
# where x = east and y = north
e_hor = np.array([ [e_bar[2,2],e_bar[2,1]],
                   [e_bar[1,2],e_bar[1,1]] ])

print('Horizontal Strain Rate Tensor: \n', e_hor)

# Calculating principle strains (+ve compressional)
# e1 = max strain
# e2 = min strain
e1 = 0.5*(e_hor[0,0] + e_hor[1,1]) + (( (0.25*(e_hor[0,0] - e_hor[1,1]))**2 + (e_hor[0,1]**2) ) ** 0.5)
e2 = 0.5*(e_hor[0,0] + e_hor[1,1]) - (( (0.25*(e_hor[0,0] - e_hor[1,1]))**2 + (e_hor[0,1]**2) ) ** 0.5)

# Calculating rotation of e2 relative to north
theta = 0.5*np.arctan2(2*e_hor[0,1],e_hor[0,0]-e_hor[1,1])

# Printing to window
print('e1,e2,theta:', e1,e2,theta*180./np.pi)

# Saving and exporting
np.savetxt('kostrov.out',np.column_stack((e1,e2,theta)),fmt='%5.3e')
