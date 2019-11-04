from math import *
import numpy as np

#constants
ProtonMass=1.6725e-24
Boltzmann=1.38066e-16
Hydrogen_MassFrac=0.76
SolarMetallicity=0.02
#set unit system
UnitMass_in_g            = 1.989e43    # 1.e10 solar masses
UnitVelocity_in_cm_per_s = 1e5         # 1 km/s
UnitLength_in_cm         = 3.085678e21 # 1 kpc
#derived
UnitTime_in_s            = UnitLength_in_cm / UnitVelocity_in_cm_per_s
UnitTime_in_Gyr          = UnitTime_in_s / (3600.*24.*365.*10.**9.)
UnitDensity_in_cgs       = UnitMass_in_g / UnitLength_in_cm**3.
UnitEnergy_in_cgs        = UnitMass_in_g * UnitLength_in_cm**2. / UnitTime_in_s**2.
G                        = 6.672e-8 / UnitLength_in_cm**3. * UnitMass_in_g * UnitTime_in_s**2.
Hubble                   = 3.2407789e-18  * UnitTime_in_s


def GetTempAdiabatic(u, gamma):
        """
        RETURNS: temperature in Kelvin for adiabatic runs (i.e. neutral gas)
        INPUT: u     : thermal energy
               gamma : adiabatic index
        """
        MeanWeight = 4./(1.+3.*Hydrogen_MassFrac) * ProtonMass 
        return (gamma-1.)*u/Boltzmann*UnitEnergy_in_cgs/UnitMass_in_g * MeanWeight

def GetTemp(u, Nelec, gamma):
	"""
	RETURNS: temperature in Kelvin
	INPUT: u     : thermal energy
	       Nelec : electron abundance
	       gamma : adiabatic index	
    	"""
	MeanWeight = 4./(1.+3.*Hydrogen_MassFrac+4.*Hydrogen_MassFrac*Nelec) * ProtonMass
	return (gamma-1.)*u/Boltzmann*UnitEnergy_in_cgs/UnitMass_in_g * MeanWeight

def GetPressure(u, rho, gamma):
	"""
	RETURNS: pressure in simulation units
	INPUT: u     : thermal energy
	       rho   : density
	       gamma : adiabatic index	
    	"""
	return (gamma-1.)*rho*u

def GetEntropy(u, rho, gamma):
	"""
	RETURNS: entropic function in simulation units
	INPUT: u     : thermal energy
	       rho   : density
	       gamma : adiabatic index	
    	"""
	return GetPressure(u, rho, gamma) / rho**gamma
	
def GetRhoCrit():
	"""
	RETURNS: critical density
    	"""
	return 3.*Hubble**2./(8.*pi*G)

def GetnH(rho, ascale, h=0.7):
	"""
	RETURNS: physical hydrogen number density in h^2 cm^-3
	INPUT: rho    : density
	       ascale : scale factor	
	       h      : Hubble constant
    	"""
	return (rho*UnitDensity_in_cgs*h**2.) * (Hydrogen_MassFrac/ProtonMass) / ascale**3.

def GetTime(ascale, OmegaM=0.27,OmegaL=0.73,h=0.7):
	"""
	RETURNS: time for given cosmology and scale factor in simulation units
	INPUT: ascale : scale factor
	       OmegaM : Omega Matter
	       OmegaL : Omega Lambda
	       h      : Hubble constant
    	"""
	aml=(OmegaM/OmegaL)**(1./3.)
	return 1./(h*Hubble) * 2./(3.* (1.-OmegaM)**0.5) * np.log((ascale/aml)**1.5 + (1.+(ascale/aml)**3.)**0.5) * UnitTime_in_Gyr

def GetLookBackTime(ascale, OmegaM=0.27,OmegaL=0.73,h=0.7):
	"""
	RETURNS: lookback time in simulation units
	INPUT: ascale : scale factor
	       OmegaM : Omega Matter
	       OmegaL : Omega Lambda
	       h      : Hubble constant
    	"""

	return GetTime(1.,OmegaM, OmegaL, h) - GetTime(ascale,OmegaM, OmegaL, h)
	
def GetCoolingRate(coor, rho):
        """
        RETURNS: cooling rate in erg s^-1 cm^-3 h^3  
        INPUT: coor : cooling rate in code units
	       rho : density in code units 
        """

	return coor*rho * UnitEnergy_in_cgs *  UnitTime_in_s**(-1.0) * UnitLength_in_cm**(-3.0) 	
