# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np 
x = np.linspace(0,np.pi/8.,10000)
def E1p(x,a):
    E1x = a*x - np.arccos((1. - a*(1. - np.cos(x)))/np.sqrt(1. - 2.*a* \
(1. - a)*(1. - np.cos(x)))) 
    return E1x
def E2p(x,a):
    E2x = a*x - np.arccos((1. - a**2*(1. - np.cos(x)))/np.sqrt(1. - a**2*\
(1. - a**2)*(1. - np.cos(x))**2)) 
    return E2x
def E3p(x,a):
    y = x/2.
    E3p_d = np.sqrt((1. - 2.*a*np.sin(y)**2*(1. - (1. - a)*np.cos(y)**2))**2 \
+ 4.*a**2*np.sin(y)**2*np.cos(y)**2*(1. + (1. - a)*np.sin(y)**2)**2)
    E3p = a*x - np.arccos((1. - 2*a*np.sin(y)**2*(1. - (1. - a)*np.cos(y)**2))/E3p_d) 
    return E3p
for a in [0.2,0.5,0.8]:
    plt.figure(figsize=(10,6)) 
    plt.plot(x,E1p(x,a), label="Schéma 1")
    plt.plot(x,E2p(x,a), label="Schéma 2")
    plt.plot(x,E3p(x,a), label="Schéma 3") 
    plt.legend()
    plt.title(r"Comparaison des erreurs de phase, avec $\alpha = $"+str(a), fontsize=14) 
    plt.xlabel(r"$\xi = k \Delta x$",fontsize=12) 
    plt.ylabel(r"$E_{\varphi} (\xi)$",fontsize=12) 
    plt.show() 
for a in [0.2,0.5,0.8]:
    plt.figure(figsize=(10,6)) 
    plt.loglog(x,np.abs(E1p(x,a))+1e-12, label="Schéma 1")
    plt.loglog(x,np.abs(E2p(x,a))+1e-12, label="Schéma 2")
    plt.loglog(x,np.abs(E3p(x,a))+1e-12, label="Schéma 3") 
    plt.legend()
    plt.title(r"Comparaison des erreurs de phase, avec $\alpha = $"+str(a), fontsize=14) 
    plt.xlabel(r"$\xi = k \Delta x$",fontsize=12) 
    plt.ylabel(r"$E_{\varphi} (\xi)$",fontsize=12) 
    plt.xlim([2e-3,5.5*1e-1])
    plt.show() 
