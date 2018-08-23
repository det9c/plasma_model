import matplotlib.pyplot as plt
import math
import sys
import scipy
import numpy as np
import random
import os
from get_ion_fraction import *


class gas:
      def __init__(self,mass,ip,zn,an,bn):
       self.mass = mass
       self.ip = ip
       self.an=an
       self.bn=bn
       self.zn=zn





hues=["black","blue","red","orange","yellow","green","purple","brown","pink"]
charge=1.0
avogadro=6.022e23
partial_pressures = np.zeros((7,1))
molecules=["n2","co2","h2o","co","ch4","nh3","h2"]
eq=1.602e-19


#Define product gas instances. Each has a mass,IP in eV,and fitting parameters from equations 26/27 in Redmer
n2=gas(28.0134,15.58,0.470340100314975,5.90991891377251,2.61299721101750)
#make list of objects to loop over later
products=[n2]



density_range=np.logspace(-6,1,8)
temp_range=np.linspace(300,12000,50)
for density in density_range:
  tvalues=[]
  econvals=[]
  ion_count=[]
  pvalues=[]
  elec_count=[]
  for temp in temp_range:
#      dens_value=float(density)
      number_of_molecules,partial_pressures=get_molecule_counts(density,temp)
      pressure_total=np.sum(partial_pressures)
      elec_total=0.0
      iproduct=-1
      for product_gas in products:
         iproduct+=1
         pval=pressure_total
         dmass=product_gas.mass/avogadro
         ionization_pot=product_gas.ip
         dnum_neutral,dnum_ion=get_ion_fraction(pval,temp,dmass,ionization_pot,number_of_molecules[iproduct,0])
#         print(dnum_neutral+dnum_ion,number_of_molecules[iproduct,0])
         elec_total+=dnum_ion
#      print(dnum_neutral,dnum_ion,number_of_molecules,dnum_neutral+dnum_ion)
      iproduct=-1
      econ=0.0
      pvalues.append(pressure_total)
      for product_gas in products:
         iproduct+=1
         pval=pressure_total
         dmass=product_gas.mass/avogadro
         ionization_pot=product_gas.ip
         dnum_neutral,dnum_ion=get_ion_fraction(pval,temp,dmass,ionization_pot,number_of_molecules[iproduct,0])
         if(elec_total != 0.0 and dnum_ion !=0.0):
           dk_elec_ion=elec_ion(temp,elec_total,dnum_ion,charge)
           zn1=product_gas.zn
           an1=product_gas.an
           bn1=product_gas.bn
           econ+=(1.0/dk_elec_ion)
         if(elec_total != 0.0 and dnum_neutral !=0.0):
           zn1=product_gas.zn
           an1=product_gas.an
           bn1=product_gas.bn
           dk_elec_neu=elec_neu(temp,zn1,an1,bn1,dnum_neutral,elec_total)
           econ+=1.0/(dk_elec_neu)
#         econ+=(1.0/dk_elec_ion)+1.0/(dk_elec_neu)
#           print(str(pressure_total)+" "+str(temp)+" "+str(econ)+"\n")
      if(econ != 0.0):
         econvals.append(eq**2/econ)
      else:
         econvals.append(0.0)
      tvalues.append(temp)
      eplot=econvals
      elec_count.append(elec_total/np.sum(number_of_molecules))
#      ion_count.append(ions_total)
  dlabel="D="+str(density)
#  if(np.sum(econvals)>0.0):
#     eplot=np.log10(econvals)
#  else:     
#      eplot=econvals
  plt.title("Electrical Conductivity of Nitrogen Plasma")
#  plt.plot(tvalues,np.log10(eplot),linewidth=1.0,ls='-',label=dlabel,markersize=4,marker="o")
#  plt.plot(tvalues,elec_count,linewidth=1.0,ls='-',label=dlabel,markersize=4,marker="o")
#  plt.ylim(-5,3)
  plt.plot(tvalues,eplot,linewidth=1.0,ls='-',label=dlabel,markersize=4,marker="o")
  plt.xlabel("Temperature (K)")
  plt.ylabel("Econ (S/m)")
  plt.legend(loc='best')
plt.show()







density_range=np.logspace(-6,-2,20)
temp_range=np.linspace(5000,10000,6)
ipass=-1
for temp in temp_range:
  ipass+=1
  dvalues=[]
  econvals=[]
  ion_count=[]
  pvalues=[]
  elec_count=[]
  for density in density_range:
#      dens_value=float(density)
      number_of_molecules,partial_pressures=get_molecule_counts(density,temp)
      pressure_total=np.sum(partial_pressures)
      elec_total=0.0
      iproduct=-1
      for product_gas in products:
         iproduct+=1
         pval=pressure_total
         dmass=product_gas.mass/avogadro
         ionization_pot=product_gas.ip
         dnum_neutral,dnum_ion=get_ion_fraction(pval,temp,dmass,ionization_pot,number_of_molecules[iproduct,0])
#         print(dnum_neutral+dnum_ion,number_of_molecules[iproduct,0])
         elec_total+=dnum_ion
#      print(dnum_neutral,dnum_ion,number_of_molecules,dnum_neutral+dnum_ion)
      iproduct=-1
      econ=0.0
      pvalues.append(pressure_total)
      for product_gas in products:
         iproduct+=1
         pval=pressure_total
         dmass=product_gas.mass/avogadro
         ionization_pot=product_gas.ip
         dnum_neutral,dnum_ion=get_ion_fraction(pval,temp,dmass,ionization_pot,number_of_molecules[iproduct,0])
         if(elec_total != 0.0 and dnum_ion !=0.0):
           dk_elec_ion=elec_ion(temp,elec_total,dnum_ion,charge)
           zn1=product_gas.zn
           an1=product_gas.an
           bn1=product_gas.bn
           econ+=(1.0/dk_elec_ion)
         if(elec_total != 0.0 and dnum_neutral !=0.0):
           zn1=product_gas.zn
           an1=product_gas.an
           bn1=product_gas.bn
           dk_elec_neu=elec_neu(temp,zn1,an1,bn1,dnum_neutral,elec_total)
           econ+=1.0/(dk_elec_neu)
#         econ+=(1.0/dk_elec_ion)+1.0/(dk_elec_neu)
#           print(str(pressure_total)+" "+str(temp)+" "+str(econ)+"\n")
      if(econ != 0.0):
         econvals.append(eq**2/econ)
      else:
         econvals.append(0.0)
      dvalues.append(density)
      eplot=econvals
      elec_count.append(elec_total/np.sum(number_of_molecules))
#      ion_count.append(ions_total)
  dlabel="T="+str(temp)
#  if(np.sum(econvals)>0.0):
#     eplot=np.log10(econvals)
#  else:
#      eplot=econvals
  plt.title("Electrical Conductivity of Nitrogen Plasma")
#  plt.plot(tvalues,np.log10(eplot),linewidth=1.0,ls='-',label=dlabel,markersize=4,marker="o")
#  plt.plot(tvalues,elec_count,linewidth=1.0,ls='-',label=dlabel,markersize=4,marker="o")
#  plt.ylim(-10,5)
  plt.plot(np.log10(dvalues),eplot,linewidth=1.0,ls='-',label=dlabel,markersize=6,marker="o",mfc=hues[ipass],color='black')
  plt.xlabel("log(Density)")# (g/cc)")
  plt.ylabel("Econ (S/m)")
#  plt.legend(loc='best')
#plt.show()



'''
exp_data=np.zeros((24,2))
with open("exp_data") as file:
      reference=file.readlines()


dens_exp=np.zeros((24,1))
pexp=101325 #exp pressure from paper is 1 atm or 101325 Pa
dkb=1.381e-23
i=0
while i<24:
      tmp=[]
      tmp=reference[i].split()
      exp_data[i,0]=float(tmp[0])
      exp_data[i,1]=float(tmp[1])*100.0
      dens_exp[i,1]= (0.01)**3*(28.01588d0/6.022e23)*pexp/(dkb*exp_data[i,0])
      i+=1
plt.plot(exp_data[0:7,0],exp_data[0:7,1],linewidth=0.0,label="Exp [D=0.00165 g/cc]",markersize=4,marker="s",color='black',mfc='red')
'''



exp_data2=np.zeros((11,5))
with open("exp_data2") as file:
      reference=file.readlines()

dkb=1.381e-23
dens_exp2=np.zeros((11,4))
i=1
while i<12:
      tmp=[]
      tmp=reference[i].split()
      exp_data2[i-1,0]=float(tmp[0])
      exp_data2[i-1,1]=float(tmp[1])*100.0
      exp_data2[i-1,2]=float(tmp[2])*100.0
      exp_data2[i-1,3]=float(tmp[3])*100.0
      exp_data2[i-1,4]=float(tmp[4])*100.0
      dens_exp2[i-1,0]= (0.01)**3*(28.01588/6.022e23)*1.0*101325/(dkb*exp_data2[i-1,0])
      dens_exp2[i-1,1]= (0.01)**3*(28.01588/6.022e23)*3.0*101325/(dkb*exp_data2[i-1,0])
      dens_exp2[i-1,2]= (0.01)**3*(28.01588/6.022e23)*10.0*101325/(dkb*exp_data2[i-1,0])
      dens_exp2[i-1,3]= (0.01)**3*(28.01588/6.022e23)*30.0*101325/(dkb*exp_data2[i-1,0])
      i+=1





i=4
ipass=-1
while i<10:
 ipass+=1
 dlabel="T(Exp)="+str(exp_data2[i,0])
 plt.plot(np.log10(dens_exp2[i,:]),exp_data2[i,1:5],linewidth=1.0,label=dlabel,markersize=6,marker="s",mfc=hues[ipass],color='black')
 i+=1
plt.legend(loc='best',fontsize=8)
plt.savefig("nitrogen.png")
plt.show()


