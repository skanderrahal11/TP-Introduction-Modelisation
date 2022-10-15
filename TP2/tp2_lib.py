import numpy as np
import scipy.sparse as spsp    ## Sparse matrix 
 

def indk(i,j, K):
    k= i +(K+1)*j
    return k

def vit(x,y,V0,L):
    V=[0,0]
    V[0]=-V0*np.sin(np.pi*x/L)*np.cos(np.pi*y/L)
    V[1]=V0*np.cos(np.pi*x/L)*np.sin(np.pi*y/L)
    return V

def source(x,y,L):
    f=256*(x/L)**2*(y/L)**2*(1-x/L)**2*(1-y/L)**2 
    return f

def MatA(L, V0, D, h, dt, K, indk, vit):

    K2 = (K+1)**2
    A = spsp.lil_matrix((K2, K2))   ## Declaration of A as a sparse matrix 
## Loop on the internal nodes
    for i in range(1,K):
        for j in  range(1,K):
            Vij=vit(h*i,h*j,L,V0)
            k = indk(i,j,K)
            ke = indk(i+1,j,K)
            ko = indk(i-1,j,K)
            kn = indk(i,j+1,K)
            ks = indk(i,j-1,K)
            A[k,k]=1-4*dt*D/h**2
            A[k,ke]=dt*D/h**2 -dt*Vij[0]/(2*h)
            A[k,ko]=dt*Vij[0]/(2*h)+dt*D/h**2
            A[k,kn]=dt*D/h**2-dt*Vij[1]/(2*h)
            A[k,ks]=dt*Vij[1]/(2*h)+dt*D/h**2
            
 ## Loop on the nodes located in x = 0 (Neumann Boundary Condition) 
    for j in range(1,K):
        V0j=vit(0,h*j,L,V0)
        k = indk(0,j,K)
        ke = indk(1,j,K)
        kn = indk(0,j+1,K)
        ks = indk(0,j-1,K)
        A[k,k]=1-3*dt*D/h**2
        A[k,kn]=dt*D/h**2 - dt*V0j[1]/(2*h)
        A[k,ks]=dt*D/h**2 + dt*V0j[1]/(2*h)
        A[k,ke]=dt*D/h**2
        
## Loop on the nodes located in x = L (Neumann Boundary Condition) 
    for j in range(1,K):
        VKj=vit(L,h*j,L,V0)
        k=indk(K,j,K)
        ko = indk(K-1,j,K)
        kn = indk(K,j+1,K)
        ks = indk(K,j-1,K)
        A[k,k]=1-3*dt*D/h**2
        A[k,kn]=dt*D/h**2 - dt*VKj[1]/(2*h)
        A[k,ks]=dt*D/h**2 + dt*VKj[1]/(2*h)
        A[k,ko]=dt*D/h**2
        
           
    return A


def secmb(K, dt, h, L, Tbas, Thaut, indk, source):

    S=np.zeros((K+1)**2)
    ##points intérieurs
    for i in range(1,K):
        for j in range(1,K):
            k=indk(i,j,K)
            Fk=source(i*h,j*h,L)
            S[k]=dt*Fk
    ##points frontière y=0 et y=L
    for i in range(K+1):
        S[i+K*(K+1)]=Thaut
    ##point frontière x=0
    for j in range(1,K):
        k=indk(0,j,K)
        Fk=source(0,j*h,L)
        S[k]=dt*Fk
    ##point frontière x=L
    for j in range(1,K):
        k=indk(K,j,K)
        Fk=source(L,j*h,L)
        S[k]=dt*Fk
        
    return S 
           
    
def matT(T, K, indk):
    T2D = np.zeros((K+1,K+1))
    for i in range(K+1):
        for j in range(K+1):
            T2D[i,j]=T[indk(i,j,K)]
        
    return T2D

def vecT(T2D,K,indk):
    T=np.zeros((K+1)**2)
    for i in range(K+1):
        for j in range(K+1):
            k=indk(i,j,K)
            T[k]=T2D[i,j]
    return T

def evolution(X, Y, Theta, K, dt, L, D, V0, Tbas, Thaut, source, vit):
       V=vit(X,Y,L,V0)
       X1=X+dt*V[0]+np.sqrt(2*D*dt)*np.random.normal(0,1,K)
       Y1=Y+dt*V[1]+np.sqrt(2*D*dt)*np.random.normal(0,1,K)
       Theta=(Y1>=1)+(0<Y1)*(1>Y1)*(Theta+dt*source(X,Y,L))
       X=np.minimum(L,np.maximum(X1,0))
       Y=np.minimum(L,np.maximum(Y1,0))
       return X,Y,Theta


def tmoy(X, Y, Theta, K, M, L):

    Tm = np.zeros([M,M])    ## mean temperature in a cell 
    nbp = np.zeros([M,M])   ## number of particles in a cell 
    e=L/M
    N=len(X)
    for k in range(0, N):
        i_cel=int(np.floor(X[k]/e))
        j_cel=int(np.floor(Y[k]/e))
        if X[k]==L :
            i_cel=M-1
        if Y[k]==L:
            j_cel=M-1
        nbp[i_cel,j_cel]+=1
        Tm[i_cel,j_cel]+=Theta[k]
    Tm=Tm/nbp
    return Tm
    

def posinit(K, L):
    X=np.random.uniform(0,L,K)
    Y=np.random.uniform(0,L,K)
    return X,Y
    
    


    
    
