# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import tp1_lib as tp

## Initialisation des données 
Lx = 1.0                     # Demi longueur de l'intervalle       
a = 2.0                      # Vitesse du signal 
tmax = 10.0                  # Temps max de la simulation (=T)

J= 201                       # Nombre de noeuds
X = np.linspace(-Lx,Lx,J)    # Abscisses des noeuds
dx = X[1]-X[0]               # Pas de maillage
cfl = 0.2                    # Nombre CFL 
dt = cfl*dx/a                # Valeur du pas de temps 

## Initialisation de la solution 
Uini = tp.Uinit2(X,J)         
U1 = Uini   # Schema 1
U2 = Uini   # Schema 2
U3 = Uini   # Schema 3


##Boucle en temps (tant que time < tmax) 
time = 0.0                     
while time < tmax :
    dt_reel = min(dt, tmax-time)   # dt_reel = dt sauf si tmax - time < dt 
    alpha = a*dt_reel/dx           # Nombre CFL effectif 
    print (' TIME ITERATION : ' , time , ' \ ' , tmax)
    U1 = tp.schema1(U1,J,alpha)    # Mise à jour de la solution schema 1
    U2 = tp.schema2(U2,J,alpha)    # Mise à jour de la solution schema 2
    U3 = tp.schema3(U3,J,alpha)    # Mise à jour de la solution schema 3
    time += dt_reel                # Incrémente le temps 
 
## Affichage des résultats
fig=plt.figure(1)
plt.title("Solution u en fonction de x pour CFL = "+str(cfl)+" et T = "+str(tmax))
plt.xlabel('x')
plt.ylabel('u')
plt.axis([-1 , 1 , -0.25, 1.25])
plt.plot(X,U1,'-g',X,U2,'-b',X,U3,'-r',X,Uini,'-k')
plt.legend(['Schema1','Schema2', 'Schema3', 'Exacte'],loc='best')
plt.grid(True)     
plt.show(block=False)
plt.pause(10.0)  ## Pause pour voir la figure avant sauvegarde
plt.savefig("u=f(x)"+"T="+str(tmax)+"_CFL ="+str(cfl)+".png")  ##sauve figure
plt.close('all')



