from astropy.io import fits
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys
import os
from scipy.optimize import fsolve
from math import pi,log
from numpy import exp,sqrt
import ast 

def load_data(fits_path):
    hdulist = fits.open(fits_path)
    data = hdulist[0].data
    # droping out the stokes dimension
    data = np.ascontiguousarray(data[0])

    # in case NaN values exist on cube
    mask = np.isnan(data)
    if np.any(mask): data = ma.masked_array(data, mask=mask)

    # map to 0-1 intensity range
    #data -= data.min()
    #data /= data.max()

    if data.shape[0]==1:
        data = np.ascontiguousarray(data[0])
        if np.any(mask):
            mask = np.ascontiguousarray(mask[0])
            data = ma.masked_array(data, mask=mask)
    return data


def load_header(fits_path):
    hdulist = fits.open(fits_path)
    hdr = hdulist[0].header
    return hdr

#os.chdir('calculations')

file = open("params.dat", "r")
contents = file.read()
par = ast.literal_eval(contents)
file.close()


data12 = load_data(par['data12'])
data13 = load_data(par['data13'])
head13 = load_header(par['data12'])
data12_int = load_data(par['data12_int'])
data13_int = load_data(par['data13_int'])

#fits.writeto("new.fits",data=data12,header=head13,overwrite=True)

R = 80.0




h = 6.6260755E-27 # erg s
k = 1.380658E-16 # erg K-1
Tbg = 2.7

nu = float(par['nu']) # Hz
nu13 = float(par['nu13']) # Hz
mu_esu = 0.11E-18 # esu
B = float(par['Brot']) # Hz
Eu = float(par['Eu']) # K
J=float(par['J'])

gj = 2.0*J+1
gk = 1.0
gi = 1.0
S = J**2/(J*(2*J+1))

def Qrot(T):
	return k*T/(h*B)
def Jv(T,nu):
	return (h*nu/k)/(exp(h*nu/(k*T))-1.0)


data_tau0 = np.copy(data12)*0
data_Tex = np.copy(data12)*0
data_tauint = np.copy(data12)*0
data_N_CO = np.copy(data12)*0
data_N_H2 = np.copy(data12)*0
for y in range(0,data12.shape[0]): #data12.shape[0]
#for y in range(150,220): #data12.shape[0]
	#print(y)
	#print(y,data12.shape[0])
	for x in range(0,data12.shape[1]): #data12.shape[1]
	#for x in range(270,370): #data12.shape[1]
		#print("x",x,"y",y,"value",data12[y,x])
		Tpeak12 = data12[y,x]
		Tpeak13 = data13[y,x]

		Tint12 = data12_int[y,x]/1000.0
		Tint13 = data13_int[y,x]/1000.0
		if Tpeak12 > 0.2:
		#if Tpeak13 > 0.1 and Tint13 > 0 and Tint12 > 0:
			def tau0rat(x):
				return Tpeak13/Tpeak12 - (1.0-exp(-x/R))/(1.0-exp(-x))
			tau0 = fsolve(tau0rat, 10.0)[0]
			
			def tauintrat(x):
				return Tint13/Tint12 - (1.0-exp(-x/R))/(1.0-exp(-x))
			tauint = fsolve(tauintrat, 10.0)[0]

			def Tex_calc(x):
				return Tpeak12 - (Jv(x,nu)-Jv(Tbg,nu))*(1.0-exp(-tau0))
			
			T0 = h*nu13/k
			#print("T0 = ",T0)

			Tex = fsolve(Tex_calc, 5.0)[0]
			Tex2 = T0/log(1.0+T0/(Tpeak12+T0/(exp(T0/Tbg)-1.0)))
			#print(Tex,Tex2)
			#print("tau0",tau0,"tauint",tauint,"Tex",Tex)

			Nthin = 3.0*h/(8.0*pi**3.0*S*mu_esu**2.0)*Qrot(Tex)/(gj*gk*gi)*(exp(Eu/(Tex)))/(exp(h*nu13/(k*Tex))-1.0)*tauint*1.0E5
			#print("Nthin=",Nthin)
			# /(Jv(Tex,nu)-Jv(Tbg,nu))*intTr

			tau_corr = tauint/(1.0-exp(-tauint))
			N = tau_corr*Nthin

			data_tau0[y,x] = tau0
			data_Tex[y,x] = Tex
			data_tauint[y,x] = tauint
			data_N_CO[y,x] = N
			data_N_H2[y,x] = N/8.0E-5

J_str = "_J="+str(int(J))+"-"+str(int(J-1))
fits.writeto("tau0"+J_str+".fits",data=data_tau0,header=head13,overwrite=True)
fits.writeto("Tex"+J_str+".fits",data=data_Tex,header=head13,overwrite=True)
fits.writeto("tauint"+J_str+".fits",data=data_tauint,header=head13,overwrite=True)
fits.writeto("N_CO"+J_str+".fits",data=data_N_CO,header=head13,overwrite=True)
fits.writeto("N_H2"+J_str+".fits",data=data_N_H2,header=head13,overwrite=True)