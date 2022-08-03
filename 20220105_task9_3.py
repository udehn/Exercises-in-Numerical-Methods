import math
import matplotlib.pyplot as plt
import notebook as notebook
#import numpy
import numpy
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

def v(x):
    if (x > 0.1 and x < 0.3):
        return 1
    else:
        return 0


#Init_data
x_low_bound = 0
x_up_bound = 2
t_low_bound = 0
t_up_bound = 2
"""
step = 5e-2
Nx = int( (x_up_bound - x_low_bound) / step )
Nt = int( (t_up_bound - t_low_bound) / step )
"""
Nx = 100
Nt = 100
stepx = (x_up_bound - x_low_bound) / Nx
stept = (t_up_bound - t_low_bound) / Nt

u = numpy.zeros((Nx, Nt))
"""
for i in range(0,Nt):
    u[0][i] = 0
    u[1][i] = 0
"""
for i in range(0,Nx):
    u[i][0] = v(x_low_bound+i*stepx)
    #print(u[i][0])
#start iteration


for j in range(0, Nt-1):
    for i in range(0, Nx):
        u[i, j+1] = (1.-stept/stepx)*u[i, j] + stept/stepx * u[((i-1) + Nx) % Nx, j]

"""
for j in range(1, Nt):
    for i in range(1, Nx-1):
        u[i][j] = (stept*u[i-1][j] + stepx*u[i][j-1]) / (stepx + stept)

"""
"""
        (ui,j - ui-1,j) / dx + (ui,j - ui,j-1) / dy = 0
        
        dy*(ui,j - ui-1,j) + dx*(ui,j - ui,j-1) = 0

        (dy+dx)*ui,j = dy*ui-1,j + dx*ui,j-1
        
        ui,j = (dy*ui-1,j + dx*ui,j-1) / (dx+dy)
"""


for i in range(0, Nx):
    for j in range(0, Nt):
        print(u[i][j] , end=" ")
    print('\n')


"""
if __name__ == '__main__':
    main()
"""



fig = plt.figure()
ax3 = plt.axes(projection='3d')

xx = np.arange(x_low_bound, x_up_bound, stepx)
yy = np.arange(t_low_bound, t_up_bound, stept)
X, Y = np.meshgrid(xx, yy)
#Z = u

ax3.plot_surface(X,Y,u,cmap='rainbow')
#ax3.contour(X,Y,Z, zdim='z',offset=0, cmap='rainbow')
plt.show()

