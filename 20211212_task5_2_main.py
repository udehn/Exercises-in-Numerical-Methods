import math
import matplotlib.pyplot as plt
import notebook as notebook
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

def f (x, y):
    return 10*(y - x**2)**2 + (1 - x)**2

def f_dx (x, y):
    return 20*(y-x**2) * (-2*x) + 2*(1-x)*(-1)

def f_dy (x, y):
    return 20*(y-x**2)

def golden(f ,*args):
    a = args[0]
    b = args[1]
    x = args[2]
    y = args[3]
    h_dx = args[4]
    h_dy = args[5]

    epsilon = 1e-2
    k = 0
    leftpoint = a+0.382*(b-a)
    rightpoint = a+0.618*(b-a)

    while True:
        k += 1
        if abs(b-a)<=epsilon:
            solve=(a+b)/2
            #print('Number of iteration: ', k)
            break
        f_left=f(x + leftpoint*h_dx, y+leftpoint*h_dy)
        f_right=f(x + rightpoint*h_dx, y+rightpoint*h_dy)
        if f_left>f_right:
            a = leftpoint
            temp1=rightpoint
            temp2=a+0.618*(b-a)
        else:
            b = rightpoint
            temp2=leftpoint
            temp1=a+0.382*(b-a)
        leftpoint=temp1
        rightpoint=temp2
    return solve

def main():
    epsilon = 1e-2
    x = 0
    y = 0
    while True:
        h_dx = -f_dx(x, y)
        h_dy = -f_dy(x, y)
        a = golden(f, 0, 10, x, y, h_dx, h_dy)
        x += a*h_dx
        y += a*h_dy
        if a < epsilon:
            print('Optimum point  is: (', x, ', ', y, ')')
            break

if __name__ == '__main__':
    main()


fig = plt.figure()
ax3 = plt.axes(projection='3d')

xx = np.arange(-5,5,0.5)
yy = np.arange(-5,5,0.5)
X, Y = np.meshgrid(xx, yy)
Z = 10*(Y - X**2)**2 + (1 - X**2)**2

ax3.plot_surface(X,Y,Z,cmap='rainbow')
#ax3.contour(X,Y,Z, zdim='z',offset=0，cmap='rainbow)
plt.show()
