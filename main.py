#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# main.py
# 
# Copyright 2020 Priscila Gutierres <priscila.gutierres@usp.br>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
# 
# 
#
 
###
# This program is intended to plot or make a table of four different 
# ode methods: euler, rk4, taylor 2 and   .
###
     
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import methods

#############################################
### defines function and initial condition
#############################################
def f(t, y):
   f0 =  y[1]
   f1 = -y[0]
    
   return np.array([f0,f1])

#############################################
### defines exactly solution
#############################################
def ye(t):
    # exact solution 
    return math.cos(t)


#############################################
#defines one step of the method
############################################
def oneStepMethod(t0,y0,T,n,phi):
     t_n = [t0];           # time interval: t in [t0,T]
     y_n = [np.array(y0)]  # initial condition
     dt   = (T-t0)/n        # time step

     while t_n[-1] < T:
         y_n.append( y_n[-1] + dt*phi(t_n[-1],y_n[-1],dt,f) )
         t_n.append(t_n[-1] + dt)
         dt = min(dt, T-t_n[-1])
     y_n = np.array(y_n)
    
     return t_n, (T-t0)/n,y_n[-1],y_n
 
############################################
### methods
############################################
def euler(t,y,dt,f):
   #define euler discretization function
   return f(t,y)
        
def rk4(t,y,dt,f):
   # define rk4 discretization function 
   k1 = f(t, y)
   k2 = f(t+dt/2, y + dt/2*k1)
   k3 = f(t+dt/2, y + dt/2*k2)
   k4 = f(t+dt, y + dt*k3)
   return 1/6*(k1 + 2*k2 + 2*k3 + k4) 
        
def m_euler(t,y,dt,f):
   #define modified euler discretization function
   k1 = f(t,y)
   k2 = f(t+dt,y+dt*k1)
   return 0.5*(k1+k2)
    
#############################################
#defines plot function
############################################
def plot_graphic(t_n, y_n, method,m):
    plt.plot(t_n, y_n[:,0], color='black', linestyle=(0,(1,1,3,1)),
             label = 'y_1(t)  (in y_1 units)')
    plt.xlabel('time t   (in units)')
    plt.ylabel('y  state variables')
    plt.title('Numerical Approximation of State Variables')
    plt.legend()
    plt.savefig(method +' '+ str(m) +' '+ 'y0'+'.png', format='png')

    plt.plot(t_n, y_n[:,1], c = 'k', label = 'y_2(t)  (in y_2 units)')
    plt.xlabel('time t   (in units)')
    plt.ylabel('y  state variables')
    plt.savefig(method +' '+ str(m) +' '+ 'y1'+'.png', format='png')
    plt.legend()
    

############################################
### main code
############################################
def main():
    
    if len(sys.argv) != 5:
        print("Not enough arguments given")
        sys.exit()
    
    T = int(sys.argv[2])
    n = int(sys.argv[3])
    m = int(sys.argv[4])
    

    if sys.argv[1] == 'euler':
        method = 'euler'
        phi = euler
        
    elif sys.argv[1] == 'rk4':
        method = 'rk4'
        phi = rk4
        
    elif sys.argv[1] == 'm_euler':
        method = 'm_euler'
        phi = m_euler
    
    
 

    t0=0; y0=[1,0]  # initial condition
 
    
    h=[0]*m  
    yn=[y0]*m
    
    for i in range(1,m+1): 
        n=16*2**(i-1)      
        
        tn, h[i-1],yn[i-1],Yn=oneStepMethod(t0,y0,T,n,phi);
        plot_graphic(tn,Yn,method,i)
        r = 0
        if i>1:
            r = abs(ye(T)-yn[i-2][0])/abs(ye(T)-yn[i-1][0])
#          
        #print("%5d %24.16e %24.16e %12.5e" % (n,h[i-1],abs(ye(T)-yn[i-1][0]),r)); 

if __name__ == "__main__":
    main()

        