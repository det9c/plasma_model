import numpy as np

def get_ion_fraction(pressure,temp,dmass,d_ip_pot,dnum_total):
#    e_joule=energy*1.602e-19
    emass=9.109e-31
    hbar=1.054e-34
    h=6.634e-34
    eq=1.602e-19
    pi=3.14159
    dkb=1.381e-23
    dkb_ev=0.000085

#    pressure=20
#    temp=12000


# T in K  P in Pa mass in g/mol ip in ev
#      temp=300d0
#      pressure=101325d0
#      dmass=28.01588d0/6.022e23
#      d_ip_pot=13.6d0
# p=nkt 
#    dnum_total=pressure/(dkb*temp)
    dion_by_dneu=2.0*np.exp(-d_ip_pot/(dkb_ev*temp))*(2*pi*emass*dkb*temp/h**2)**1.5*(dkb*temp) /pressure
#    print(dion_by_dneu,dion_by_dneu/(1.0+dion_by_dneu))
#    dnum_neutral=dnum_total/(1.0e0+dion_to_neutral)
#    dnum_neutral=dnum_total-dnum_total*dion_by_dnum_total
    dnum_ion=dnum_total*dion_by_dneu/(1.0+dion_by_dneu)
    dnum_neutral=dnum_total-dnum_ion

    dnum_ion=dnum_ion*4000
    dnum_neutral=dnum_total-dnum_ion
#    if(dnum_ion > dnum_total):
#        print(dnum_ion,dnum_total,dnum_ion/3000)
#        print("error in counts")
#        dnum_ion=dnum_total*dion_by_dneu/(1.0+dion_by_dneu)
#        dnum_neutral=dnum_total-dnum_ion

#    print(dion_by_dnum_total,np.exp(-d_ip_pot/(dkb_ev*temp)),(2*pi*emass*dkb*temp/h**2)**1.5*(dkb*temp)*1.0e0/pressure)
#    print(-d_ip_pot,-d_ip_pot/(dkb_ev*temp),np.exp(-d_ip_pot/(dkb_ev*temp)),temp,pressure)
#    density=(0.01e0)**3*dnum_total*dmass
    return dnum_neutral,dnum_ion



def elec_ion(temp,dnum_elec,dnum_ion,z):
      ee_scatter=1.0 - 0.440*np.exp(-0.0458*z)/(z**0.369) #scale factor. eqn 20
      pi=3.14159e0
      eo=8.8541878e-12
      eq=1.602e-19
      dkb=1.381e-23
      emass=9.109e-31
      hbar=1.054e-34
      prefac=(2**5.5e0)*np.sqrt(pi)*6.0e0*eo**2
      prefac=prefac*dnum_elec/(3.0*z**2*eq**4*np.sqrt(emass)*dnum_ion)
      prefac=prefac*(dkb*temp)**1.5
      dkappa_sq=eq**2*dnum_elec/(eo*dkb*temp)
      bn=8.0e0*np.sqrt(6.0)*emass*dkb*temp/(hbar**2*dkappa_sq)
      dlambda=0.5e0*( np.log(1+bn)- (bn/(1.0+bn))  )
      dk_elec_ion=prefac*ee_scatter/dlambda
      return dk_elec_ion





def elec_neu(temp,z,ain,b_in,dnum_neutral,dnum_elec):
      pi=3.14159e0
      eo=8.8541878e-12
      eq=1.602e-19
      dkb=1.381e-23
      emass=9.109e-31
      hbar=1.054e-34
      h=6.634e-34

      e_joule=temp*dkb
      a=ain*.52917e-10
      b=b_in*.52917e-10

      xn=8.0e0*b*b*emass*e_joule/(hbar*hbar)
      term=np.log(xn+1)/2.0e0
      term=term-(3.0*a*a + 6.0*a*b - 2.0*b*b)/(6.0* a* a)
      term=term+(4.0*a*b + a*a)/(2.0*a*a*(xn+1))
      term=term-(b*b+a*b)/(a*a*(xn+1)**2)
      term=term+2.0e0*b*b/(3.0e0*a*a*(xn+1)**3)
 #     cross_sec=term*z**2*eq**4/(16.0e0*pi*eo**2*e_joule*e_joule)

      prefac=(2**5.5e0)*np.sqrt(pi)*6.0e0*eo**2
      prefac=prefac*dnum_elec*(dkb*temp)**1.5
      prefac=prefac/(3.0*z**2*eq**4*np.sqrt(emass)*dnum_neutral)
      

      argu=8.0e0*np.sqrt(6.0)*emass*dkb*temp*b*b/(hbar**2)
      term=np.log(argu+1)/2.0e0
      term=term-(3.0*a*a + 6.0*a*b - 2.0*b*b)/(6.0* a* a)
      term=term+(4.0*a*b + a*a)/(2.0*a*a*(argu+1))
      term=term-(b*b+a*b)/(a*a*(argu+1)**2)
      term=term+2.0e0*b*b/(3.0e0*a*a*(argu+1)**3)

      dk_elec_neu=prefac/term

      return dk_elec_neu




def get_molecule_counts(dens_value,temp):
    counts = np.zeros((1,1))
    partial_pressures = np.zeros((1,1))
    dkb=1.381e-23
    num_h2=6.022E+23*1.0*dens_value*100**3/28.0134
    counts[0,0]=num_h2
    partial_pressures[0,0]=num_h2*dkb*temp
    return counts,partial_pressures
    



















