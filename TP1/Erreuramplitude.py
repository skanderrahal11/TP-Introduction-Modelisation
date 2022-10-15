# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
x=np.linspace(0,np.pi/8,10000)
def E1(x,a):
    Ea1=1-np.sqrt(1-2*a*(1-a)*(1-np.cos(x)))
    return Ea1
def E2(x,a):
    Ea2=1-np.sqrt(1-a**2*(1-a**2)*(1-np.cos(x))**2)
    return Ea2
def E3(x,a):
    y=x/2
    u=1-2*a*np.sin(y)**2*(1-(1-a)*np.cos(y)**2)
    v=4*a**2*np.sin(y)**2*np.cos(y)**2*(1+(1-a)*np.sin(y)**2)**2
    Ea3=1-np.sqrt((u**2+v))
    return Ea3
for a in [0.2,0.5,0.8]:
    figure=plt.figure(figsize=(10,6))
    plt.plot(x,E1(x,a),label="Schéma 1")
    plt.plot(x,E2(x,a),label="Schéma 2")
    plt.plot(x,E3(x,a),label="Schéma 3")
    plt.legend()
    plt.title("Les erreurs d'amplitude pour alpha ="+str(a), fontsize=13)
    plt.xlabel(r"$\xi = k \Delta x$", fontsize=12)
    plt.ylabel(r"$E_{A} (\xi)$",fontsize=12)
    plt.show()
    
for a in [0.2,0.5,0.8]:
     plt.figure(figsize=(8.5,3.5))
     plt.semilogy(x,E1(x,a), label="Schéma 1") 
     plt.semilogy(x,E2(x,a), label="Schéma 2") 
     plt.semilogy(x,E3(x,a), label="Schéma 3") 
     plt.legend()
     plt.title(r"Comparaison des erreurs d'amplitude, avec $\alpha = $"+str(a), fontsize=14) 
     plt.xlabel(r"$\xi = k \Delta x$",fontsize=12) 
     plt.ylabel(r"$E_{A} (\xi)$",fontsize=12) 
     plt.show()
   
