#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 09:38:05 2020

@author: priscila
"""
from math import exp
import matplotlib.pyplot as plt
from numpy import arange


h = 0.5
TOL = 0.1
beta = 0.5

r = 2
K = 10**5
T = 20.0
y0 = 1
t0 = 0
t = t0

beta = 0.5
alpha = 1.5

def ye(t):
    # exact solution 
    return K/(1+(K-1)*exp(-r*t))


def m_euler(tk, y,dt,f):
   #define modified euler discretization function
   k1 = f(tk,y)
   k2 = f(tk+dt,y+dt*k1)
   return 0.5*(k1+k2)

phi = m_euler

def f(t, y):  
   return r*(1 - y/K)*y

def erro(t,y,h):
    k1 = f(t,y)
    k2 = f(t+h/2,y+h*k1)
    k3 = f(t+h,y-h*k1+2*h*k2)
    
    return 1/6*(k1-2*k2+k3)


erro_list = []

print("t      y      h     erro")

hlist = [0.5,1, 1.5]

for h in hlist:
    
    yn = [y0]
    tn = [t0]

    while (tn[-1] < T):
#    erro_tol = erro(tn[-1],yn[-1],h)
    
#    while (erro_tol > TOL):
#        h = beta * h
#        erro_tol = erro(tn[-1],yn[-1],h)
    
        erro_tol = abs(yn[-1] - ye(tn[-1]))    
        erro_list.append(erro_tol)
        yn.append( yn[-1] + h*phi(tn[-1],yn[-1],h,f) )
        tn.append(tn[-1] + h)
 
        print(str(tn[-1]) + "  " + str(yn[-1]) + " " + str(erro_tol ))
        plt.xlabel('time t   (in units)')
        plt.ylabel('y  state variable')
        plt.plot(tn, yn, linestyle=(':'))
        print("erro minimo" + " "+ str(min(erro_list))) 
        print("erro máximo" +" "+ str(max(erro_list)))
#    if (tn[-1] < T):
#        h = h*alpha


#    else:
#         break
     

plt.plot(list(arange(0,T,0.25)),list(map(ye,list(arange(0,T,0.25)))),color='black',linestyle='-',linewidth=1.0)
plt.show()