# -*- coding: cp1252 -*-
import numpy as np

def Uinit(X,J):
    Uinit=np.zeros(J)
    for i in range(J):
        if X[i]>=(-1/4) and X[i]<=(1/4):
            Uinit[i]=256*(X[i]-1/4)**2 * (X[i]+1/4)**2
        else:
            Uinit[i]=0
    return Uinit

def schema1(U0,J,alpha):
    U1=np.zeros(J)
    U1[0]=U0[0]-alpha*(U0[0]-U0[J-2])
    for i in range(1,J):
        U1[i]=U0[i]-alpha*(U0[i]-U0[i-1])
    return U1

def schema2(U0,J,alpha):
    U1=np.zeros(J)
    U1[0]=U0[0]-(alpha/2)*(U0[1]-U0[J-2])+((alpha**2)/2)*(U0[1]-2*U0[0]+U0[J-2])
    U1[J-1]=U0[J-1]-(alpha/2)*(U0[1]-U0[J-2])+((alpha**2)/2)*(U0[1]-2*U0[J-1]+U0[J-2])
    for i in range(1,J-1):
        U1[i]=U0[i]-(alpha/2)*(U0[i+1]-U0[i-1])+((alpha**2)/2)*(U0[i+1]-2*U0[i]+U0[i-1])
    return U1

def schema3(U0,J,alpha):
    U1=np.zeros(J)
    U1[0]=U0[0]-alpha*(U0[0]-U0[J-2])-alpha*(1-alpha)/4*(U0[1]-U0[0]-U0[J-2]+U0[J-3])
    U1[1]=U0[1]-alpha*(U0[1]-U0[0])-alpha*(1-alpha)/4*(U0[2]-U0[1]-U0[0]+U0[J-2])
    U1[J-1]=U0[J-1]-alpha*(U0[J-1]-U0[J-2])-alpha*(1-alpha)/4*(U0[1]-U0[J-1]-U0[J-2]+U0[J-3])
    for i in range(2,J-1):
        U1[i]=U0[i]-alpha*(U0[i]-U0[i-1])-alpha*(1-alpha)/4*(U0[i+1]-U0[i]-U0[i-1]+U0[i-2])
    return U1

def Uinit2(X,J):
    Uinit2=np.zeros(J)
    for i in range(J):
        if X[i]>=(-1/4) and X[i]<=(1/4):
            Uinit2[i]=1
        else:
            Uinit2[i]=0
    return Uinit2