import matplotlib.pyplot as plt
from math import exp
from numpy import arange

r = 2
K = 10**5
T = 20.0


def ye(t):
    # exact solution 
    return K/(1+(K-1)*exp(-r*t))




def f(t, y):  
   return r*(1 - y/K)*y

TOL = 0.01
hmax = 2
hmin = 0.001

h = 1

FLAG = 1

t = 0
tn = [t]

a = 0
b = 20

w = 1
w_list = [w]
t_list = [t]
h_list = []
erro_list = []

while (FLAG == 1):
    k1 = h*f(t,w)
    k2 = h*f(t+1/4*h, w+1/4*k1)
    k3 = h*f(t+3/8*h, w+3/32*k1+9/32*k2)
    k4 = h*f(t+12/13*h,w+1932/2197*k1-7200/2197*k2+7296/2197*k3)
    k5 = h*f(t+h,w+439/216*k1-8*k2+3680/513*k3-845/4104*k4)
    k6 = h*f(t+1/2*h, w-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4 - 11/40*k5)

    R = abs(1/360*k1 - 128/4275*k3 - 2197/75240*k4 +1/50*k5 + 2/55*k6)
    if(R <= TOL):
        erro_list.append(R)
        h_list.append(h)
        w = w + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 -1/5*k5
        w_list.append(w)
        t_list.append(t)
        t = t+h


    delta = 0.84*(TOL/R)**(1/4)
    
    if(delta <= 0.1):
        h = 0.1*h
    elif(delta > 4):
        h = 4*h
    else:
        h = delta*h
       
    if (h > hmax):
        h = hmax
    
    if (t >= b ):
        FLAG = 0
    elif (t+h > b):
        h = b-t
    elif (h < hmin):
        FLAG = 0
        print("h minimo excedido")
        break;

print("")
print("h mínimo: "+ str(min(h_list)))
print("h maximo: "+ str(max(h_list)))
print("erro minimo" + " "+ str(min(erro_list)))
print("erro máximo" +" "+ str(max(erro_list)))


plt.xlabel('time t   (in units)')
plt.ylabel('y  state variable')
plt.plot(t_list, w_list, linestyle=(':'),c='r')
plt.plot(list(arange(0,T,0.25)),list(map(ye,list(arange(0,T,0.25)))),color='black',linestyle='-',linewidth=1.0)
plt.show()
