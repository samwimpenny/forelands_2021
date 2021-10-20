### A bunch of functions
### for doing various things
### with moment tensors.

# Modules
import os
from numpy import sin,cos,arcsin,arccos,tan,arctan,arctan2,degrees,radians,pi
from numpy import sqrt
import numpy as np


def mw2m0(mw):
	from numpy import power
	inter=1.5*(float(mw)+6.)
	m0=power(10,inter)
	return m0
def m02mw(m0):
	from numpy import log10
	inter=log10(m0)
	mw=(2./3)*inter-6.
	return mw
def CMTmw2m0(mw):
	from numpy import power
	inter=1.5*(float(mw)+10.7)
	m0=power(10,inter)
	return m0
def CMTm02mw(m0):
	from numpy import log10
	inter=log10(m0)
	mw=(2./3)*inter-10.7
	return mw

def slipvec(strike,dip,rake):
	""" 
	Function to find slip vector azimuth from strike, dip,rake (all in degrees)
	Returns azimuth of slip vector in degrees.
	
	"""
	if rake==180.:
		azimuth=strike-180.
	else:
		#Separates horizontal component of slip vector 
		# into strike-parallel and strike-perpendicular components
		strike_par=cosd(rake)
		strike_perp=sind(rake)*cosd(dip)
		#Find angle of slip vector from strike (0 is strike-parallel)
		angle=arctan2(strike_perp,strike_par)
		azimuth=strike-degrees(angle)
	if azimuth<0.:
		azimuth=azimuth+360.
	elif azimuth>360.:
		azimuth=azimuth-360.
	return azimuth

def svplunge(strike, dip, rake):
	"""
	Function to find plunge of slip-vector from strike, rake and dip.
	Arguments: strike, dip, rake (all in degrees).
	"""
	### Change angles to radians
	sr = radians(strike)
	dr = radians(dip)
	rr = radians(rake)
	### 
	per_sv = sin(rr)*cos(dr)
	par_sv = cos(rr)
	horr_sv = sqrt(per_sv**2 + par_sv**2)
	plunge = arccos(horr_sv)
	return degrees(plunge)

def svproj(strike, dip, rake, proj_azimuth):
	### Change angles to radians
	sr = radians(strike)
	dr = radians(dip)
	rr = radians(rake)
	### 
	per_sv = sin(rr) * cos(dr)
	par_sv = cos(rr)
	horr_sv = sqrt(per_sv**2 + par_sv**2)
	ver_sv = sin(rr) * sin(dr)
	###
	az = slipvec(strike, dip, rake)
	azdiff = abs(proj_azimuth - az)
	if azdiff > 180.:
		azdiff = 360. - azdiff
	if azdiff > 90.:
		azdiff = 180. - azdiff
	new_hor_sv = horr_sv * cos(radians(azdiff))
	proj_plunge = arctan(ver_sv/new_hor_sv)
	
	return degrees(proj_plunge)
	

def taxis(strike, dip, rake):
	""" 
	Calculates the trend and plunge of the T-axis of an earthquake from 
	strike, dip and rake (all in degrees).
	Returns trend (azimuth, degrees) and plunge (downwards positive, degrees)
	Andy Howell, Jul 2015
	"""
	### Change angles to radians
	sr = radians(strike)
	dr = radians(dip)
	rr = radians(rake)
	### Length of circle of radius sqrt(2)/2 (T- axis always 45deg from plane)
	circ_length = sin(pi/4.)
	### Separate into strike-parallel and strike-perpendicular components
	strike_par = circ_length * cos(rr)
	### Horizontal (y-) component of strike-perpendicular component
	strike_perp_h = circ_length * (sin(rr)*cos(dr) - sin(dr))
	### Vertical (z-) component of strike-perpendicular component
	strike_perp_v = circ_length * (sin(rr)*sin(dr) + cos(dr))
	
	### Find magnitude of horizontal component of T-axis vector
	mag_hor = sqrt(strike_par**2 + strike_perp_h**2)
	### Find plunge of T-axis
	plunge1 = arctan(strike_perp_v/mag_hor)
	### Determine whether plunge is in an up- or down-direction
	if plunge1 <= 0.:
		### If plunge negative, make positive (cartesian is upwards-positive)
		plungedeg = degrees(-1. * plunge1)
		### Find horizontal azimuth of T-axis (same direction as slip-vector)
		az_rel_strike = degrees(arctan2(strike_perp_h, strike_par))
		### Subtract from strike to get slip-vector azimuth
		azimuth=strike-az_rel_strike
		### Write as bearing
		if azimuth<0.:
			azimuth=azimuth+360.
		elif azimuth>360.:
			azimuth=azimuth-360.
			
	else:
		### If plunge positive, turn to degrees
		plungedeg = degrees(plunge1)
		### Find (back-) azimuth in degrees
		az_rel_strike_rev = degrees(arctan2(strike_perp_h, strike_par))
		### Flip azimuth by 180deg to get downwards-positive azimuth
		if az_rel_strike_rev <= 0.:
			az_rel_strike = az_rel_strike_rev + 180.
		else:
			az_rel_strike = az_rel_strike_rev - 180.
		
		azimuth=strike-az_rel_strike
		if azimuth<0.:
			azimuth=azimuth+360.
		elif azimuth>360.:
			azimuth=azimuth-360.
	return azimuth, plungedeg

def paxis(strike, dip, rake):
	""" 
	Calculates the trend and plunge of the P-axis of an earthquake from 
	strike, dip and rake (all in degrees).
	Returns trend (azimuth, degrees) and plunge (downwards positive, degrees)
	Pasted from T-axis code and altered
	Andy Howell, Jul 2015
	"""
	### Change angles to radians
	sr = radians(strike)
	dr = radians(dip)
	rr = radians(rake)
	### Length of circle of radius sqrt(2)/2 (P- axis always 45deg from plane)
	circ_length = sin(pi/4.)
	### Separate into strike-parallel and strike-perpendicular components
	strike_par = circ_length * cos(rr)
	### Horizontal (y-) component of strike-perpendicular component
	strike_perp_h = circ_length * (sin(rr)*cos(dr) + sin(dr))
	### Vertical (z-) component of strike-perpendicular component
	strike_perp_v = circ_length * (sin(rr)*sin(dr) - cos(dr))
	
	### Find magnitude of horizontal component of P-axis vector
	mag_hor = sqrt(strike_par**2 + strike_perp_h**2)
	### Find plunge of P-axis
	plunge1 = arctan(strike_perp_v/mag_hor)
	### Determine whether plunge is in an up- or down-direction
	if plunge1 <= 0.:
		### If plunge negative, make positive (cartesian is upwards-positive)
		plungedeg = degrees(-1. * plunge1)
		### Find horizontal azimuth of P-axis (same direction as slip-vector)
		az_rel_strike = degrees(arctan2(strike_perp_h, strike_par))
		### Subtract from strike to get slip-vector azimuth
		azimuth=strike-az_rel_strike
		### Write as bearing
		if azimuth<0.:
			azimuth=azimuth+360.
		elif azimuth>360.:
			azimuth=azimuth-360.
			
	else:
		### If plunge positive, turn to degrees
		plungedeg = degrees(plunge1)
		### Find (back-) azimuth in degrees
		az_rel_strike_rev = degrees(arctan2(strike_perp_h, strike_par))
		### Flip azimuth by 180deg to get downwards-positive azimuth
		if az_rel_strike_rev <= 0.:
			az_rel_strike = az_rel_strike_rev + 180.
		else:
			az_rel_strike = az_rel_strike_rev - 180.
		
		azimuth=strike-az_rel_strike
		if azimuth<0.:
			azimuth=azimuth+360.
		elif azimuth>360.:
			azimuth=azimuth-360.
	return azimuth, plungedeg
	

def strikedipsv2rake(strike,dip,slipvec):
	""" 
	Program to perform the 'inverse' of slipvec.
	Arguments: strike, dip azimuth of slip vector (all degrees)
	Returns rake (degrees)
	"""
	
	
	angle = strike - slipvec
	if angle < -180.:
		angle = 360. + angle
	elif angle > 180.:
		angle = angle - 360.
	
	if angle == 90.:
		rake = 90.
	elif angle == -90.:
		rake = -90.
	else:
		strike_par = cosd(angle)
		strike_perp = sind(angle)/cosd(dip)
		rake = degrees(arctan2(strike_perp,strike_par))
	return rake

def slipdip2rake(strike,dip,slipvec):
	""" 
	Does exactly the same as the previous function (strikedipsv2rake)
	but has a stupid name!
	Kept here to avoid going through changing scripts
	
	Program to perform the 'inverse' of slipvec.
	Arguments: strike, dip azimuth of slip vector (all degrees)
	Returns rake (degrees)
	"""
	
	
	angle = strike - slipvec
	if angle < -180.:
		angle = 360. + angle
	elif angle > 180.:
		angle = angle - 360.
	
	if angle == 90.:
		rake = 90.
	elif angle == -90.:
		rake = -90.
	else:
		strike_par = cosd(angle)
		strike_perp = sind(angle)/cosd(dip)
		rake = degrees(arctan2(strike_perp,strike_par))
	return rake

def otherplane(strike1,dip1,rake1):
	""" 
	Function to calculate strike, dip, rake of auxiliary plane and
	Inputs: array containing strike,dip,rake of fault plane (all degrees)
	Outputs are strike,dip,rake of other plane
	"""
	### Convert to radians
	s1=radians(strike1); d1=radians(dip1); r1=radians(rake1)
	### Calculate new strike,dip,rake (routine from Shearer, seismology)
	d2=arccos(sin(r1)*sin(d1))
	
	sr2=cos(d1)/sin(d2)
	cr2=-sin(d1)*cos(r1)/sin(d2)
	r2=arctan2(sr2,cr2)
	
	s12=cos(r1)/sin(d2)
	c12=-1./(tan(d1)*tan(d2))
	s2=s1-arctan2(s12,c12)
	
	### Convert back to degrees
	strike2=degrees(s2); dip2=degrees(d2); rake2=degrees(r2)
	
	### Consistency with strike convention:
	if dip2>90.:
		strike2=strike2+180.
		dip2=180.-dip2
		rake2=360.-rake2
		
	### Check strike in range of interest
	if strike2>360.:
		strike2=strike2-360.
	
	return strike2,dip2,rake2

def sdr2mtvect(strike_deg, dip_deg, rake_deg):
	"""
	Convert strike, dip, rake (degrees) into mt...
	mt = [ mrr, mtt, mpp, mtp, mrp, mrt ] where r=UP, t=S, p=E
	Pikeyed off PythMT module on github.
	"""
	### Convert degrees to radians
	strike = np.radians(strike_deg)
	dip = np.radians(dip_deg)
	rake = np.radians(rake_deg)
	
	is2 = 1/sqrt(2.0)

	mrr =  is2*( sin(2*dip) * sin(rake) );

	mtt = -is2*( sin(dip)  * cos(rake) * sin(2*strike) +     
					sin(2*dip) * sin(rake) * sin(strike)**2 );

	mpp =  is2*( sin(dip)  * cos(rake) * sin(2*strike) -
					sin(2*dip) * sin(rake) * cos(strike)**2 );

	mtp = -is2*( sin(dip)  * cos(rake) * cos(2*strike) + 
					0.5*sin(2*dip) * sin(rake) * sin(2*strike) );

	mrp =  is2*( cos(dip)  * cos(rake) * sin(strike)  -
					cos(2*dip) * sin(rake) * cos(strike));

	mrt = -is2*( cos(dip)  * cos(rake) * cos(strike)  +    
					cos(2*dip) * sin(rake) * sin(strike) );

	return np.array([ mrr, mtt, mpp, mtp, mrp, mrt ] )

def sdr2smt( strike, dip, rake ):
	"""
	Convert strike, dip, rake (radians) into unit norm mt matrix...
	mt = [ mrr, mrt, mrp ;
			mrt, mtt, mtp ;
			mrp, mtp, mpp ],  where r=UP, t=S, p=E  
	"""
	# get in vect form
	mtvect = sdr2mtvect( strike, dip, rake )

	# convert to matrix form
	return voigt2mat( mtvect )

def voigt2mat( v ):
	"""
	convert from vector with voigt notation to a matrix
	"""
	return np.array([[ v[0], v[5], v[4] ],
				     [ v[5], v[1], v[3] ],
				     [ v[4], v[3], v[2] ]])

def rotate_moment(orig_mt, azimuth):
	"""
	Function to rotate moment tensor so that its x-axis is aligned with the 
	direction to the station from the event. orig_mt is a 3 x 3 numpy array.
	azimuth is a bearing.
	"""
	x_orig = 90.
	angle = azimuth - x_orig
	
	cosa = np.cos(np.radian(angle))
	sina = np.sin(np.radian(angle))
	clockwise_rotation = np.array([ [cosa, sina, 0], 
									[-sina, cosa, 0],
									[0, 0, 1] ])
	
def cosd(x):
  y = np.cos(np.radians(x))
  return y

def sind(x):
  y = np.sin(np.radians(x))
  return y
  
def arccosd(x):
  y = np.degrees(np.arccos(x))
  return y

def arcsind(x):
  y = np.degrees(np.arcsin(x))
  return y
  
def tand(x):
  y = np.tan(np.radians(x))
  return y

def arctand(x):
  y = np.degrees(np.arctan(x))
  return y
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


	
	
