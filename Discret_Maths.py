from math import exp, log 
import numpy
import rkf45
import matplotlib.pyplot as plt
from numpy import zeros, matrix
from scipy import linalg
from scipy.optimize import fsolve
from math import tan,tanh
import scipy.integrate as integrate
import scipy.special as special
__author__ = 'iSecret'


p = [1, 0.1, 0.01, 0.0001, 0.000001]

E = numpy.matrix([[1,0,0],
                [0,1,0],
                  [0,0,1]])

A = numpy.matrix([[441,378,162],
                [378,414,216],
                  [162,216,144]])
B= [4203,4194,2106]
C= [0,0 ,0 ]

#differential
def calc(t,y0):
    global k
    global q
    y= zeros(2)
    #y[0]=y0[0]
    y[0]=y0[1]
    y[1]=(-q/(4*k**2)*((1-k**2)*y0[0]+(2*k**2*y0[0]**3)))
    
    return y

#function to fin max solution
def func(x):
    return tan(x)+tanh(x)
# Function
def f(x) :
    ff=tan(x)+tanh(x)
    return ff
# integral Parameter
def Fa(x):
    Ffa=(1-numpy.exp(-0.6*x))/(x*(x+1))
    return Ffa
    
# A parameter
def Ax(value):
    Aa=pow((value-0.36691682),4);
    return Aa
#rkf45
def calcF(f, eps, n, xp, t, y0, step):

    flag = 1
    t_start = t[0]
    t_out = t[1]
    plt.style.use('ggplot')
    for i_step in interval(t_start, t_out + step, step):
        x_end, xp_end, t_end, flag_end = rkf45.r8_rkf45(f=calc, neqn=n, y=y0, yp=xp, t=t_start, tout=i_step, relerr=eps, abserr=0, flag=flag)
        if 0 == 0:
            plt.plot(i_step, [x_end[0]], 'ro-')
            plt.plot(i_step, [x_end[1]],'bo-')
                     
          
            print ('t =',str(round(i_step, 2)), ' x1 =', str(round(x_end[0],5)), '  x2 =', str(round(x_end[1],5)), ' flag =', flag_end, )
    plt.axis([t[0], t[1]+0.5, -15, 15])
    plt.show()
# interval
def interval(fro, to, step):
    while fro < to:
        yield fro
        fro += step

T1 = numpy.linalg.inv(A)
C=T1.dot(B)
k=C[0,0]
q=C[0,1]
T=C[0,2]
x0=-5
k=k+pow(10,-5)
print ('\k =', str(round(k,5)), '1-- -q-- =', str(round(q,5)),'-- -C-- =', str(round(T,5)))
my_list=range(-10,10)
for i in my_list:
    if x0>0:
        break
    else:
        x0 = fsolve(func,i)
i=0
R=0
B=0.01691317*x0
while i<=1:
    result=integrate.quadrature(Fa, i, 1) 
    
    R=R+result[0]
    i=i+0.1
print ("%.2f"% R,"  || ")

A=Ax(R)
print ("%.2f"% A,"  || ")
print ('\B =',B)
j=0
t = [0, T]
y0 = [A, B]
xp = calc(t[0],y0)
step = 0.2

calcF(calc, 1e-3, 2, xp, t, y0, step)


 

