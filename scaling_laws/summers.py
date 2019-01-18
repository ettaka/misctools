import scipy.misc
import math

def Bc2(T,Bc20,Tc0):
    return Bc20*(1.-(T/Tc0)**2.)*(1-0.31*(T/Tc0)**2.*(1.-1.77*math.log(T/Tc0)))
def Jc(B,T,C,Bc20,Tc0):
    return C/math.sqrt(B)*(1.-B/Bc2(T,Bc20,Tc0))**2.*(1.0-(T/Tc0)**2.)**2.

def getJcB(T,C,Bc20,Tc0):
    def JcB(B):
        return Jc(B,T,C,Bc20,Tc0)
    return JcB

JcB=getJcB(4.2,3.7845e10,27.6,17.)

print JcB(12.)*1e-6
print scipy.misc.derivative(JcB, 12.,1e-6)*1e-6

