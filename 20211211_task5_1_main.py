import math
import matplotlib.pyplot as plt
import numpy as np

def f ( x ):
    L = (x**2 - 6*x +12) / (x**2 + 6*x +12)
    return L

def golden(f ,*args):
    if not('a' in dir()):
        a=[]
        b=[]
    a.append(args[0])
    b.append(args[1])
    L=1e-3
    n=80

    leftpoint = a[0]+0.382*(b[0]-a[0])
    rightpoint = a[0]+0.618*(b[0]-a[0])

    for k in range(0,n):
        if abs(b[k]-a[k])<=L:
            solve=(a[k]+b[k])/2
            print('Number of iteration: ', k)
            break
        f_left=f(leftpoint)
        f_right=f(rightpoint)
        if f_left>f_right:
            a.append(leftpoint)
            b.append(b[k])
            temp1=rightpoint
            temp2=a[k+1]+0.618*(b[k+1]-a[k+1])
        else:
            a.append(a[k])
            b.append(rightpoint)
            temp2=leftpoint
            temp1=a[k+1]+0.382*(b[k+1]-a[k+1])
        leftpoint=temp1
        rightpoint=temp2

    return solve


print('Optimum point  is: ', golden(f,0,20))


x = np.arange(0, 20, 0.1)
y = (x**2 - 6*x +12) / (x**2 + 6*x +12)
plt.plot(x, y)
plt.show()
