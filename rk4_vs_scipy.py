import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import pickle

def f(x,y):  
  return x*y**0.5
  
def rk4(f,x0,y0,h,xn):
  n=(xn-x0)/h
  x=np.arange(x0,xn+h,h)
  y=[y0]
  for i in range(int(n)):
    k1=f(x[i],y[i])
    k2=f(x[i]+h/2,y[i]+k1*h/2)
    k3=f(x[i]+h/2,y[i]+k2*h/2)
    k4=f(x[i]+h,y[i]+k3*h)
    c=y[i]+((k1+2*k2+2*k3+k4)*h)/6
    y.append(c)
  return [x,y] #a list which stores sublist for x and y 

f1=open('rk4_ic_manav.txt','r')
a=f1.readline() #'0, 1/n'
y0=int(a[3])   
f1.close()

def rk4extractor():
    ylist=[]
    hlist=[0.01,0.05,0.1,0.2,0.5]
    for i in hlist:
        ylist.append(rk4(f,0,y0,i,1))
    return ylist

ylist=rk4extractor()

def out_creator(ylist): 
    f2=open('rkr_out_manav.txt','w')
    for i in range(len(ylist)):
        c=list(ylist[i][1])
        c=[str('%.5f'%i) for i in c]
        row=','.join(c)
        f2.write(f"{row}\n")
    f2.close()

out_creator(ylist)

def comparer(ylist):
    yscipy=[]
    for i in ylist:
        sol=si.odeint(f,y0,i[0])
        yscipy.append(sol)
    return yscipy

yscipy=comparer(ylist) #list with sublists of y solutions

for i in range(len(ylist)):
    print('xlist is:',ylist[i][0])
    print('ylist from rk4 is:',ylist[i][1])
    print('ylist from scipy is:',yscipy[i])
    print('*'*30)
    
def solution(x):
    return (4+x[i]*x[i])**2/16.0
    
    
def plotter(ylist,yscipy):
    hlist=[0.01,0.05,0.1,0.2,0.5]
    for i in range(len(ylist)):
        plt.plot(ylist[i][0],ylist[i][1],'.-b',label=f"rk4 for h={hlist[i]}")
        plt.plot(ylist[i][0],yscipy[i],'.-r',label=f"rk4 scipy for h={hlist[i]}")
    plt.xlabel(r'x values')
    plt.ylabel(r'y solutions')
    plt.grid(True)
    plt.legend()
    plt.title('RK4 Method for ODE')
    plt.show()

plotter(ylist,yscipy)

        
        
        
        
        
    
