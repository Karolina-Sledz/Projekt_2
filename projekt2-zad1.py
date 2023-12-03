import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

data = np.array([[1.0, 3.0], [2.0, 1.0], [3.5, 4.0], [5.0, 0.0], [6.0, 0.5], [9.0 , -2.0], [9.5, -3.0]])
x = np.linspace(0, 10, 100)

N = data.shape[0]-1
h = np.zeros(N)
b = np.zeros(N)
u = np.zeros(N)
v = np.zeros(N)
z = np.zeros(data.shape[0]+1)

for i in range(N):
    h[i] = data[i+1,0] - data[i,0]
    b[i] = (6/h[i]) * (data[i+1,1] - data[i,1])
    
u[1] = 2*(h[0] + h[1])
   
for i in range(N):
    u[i] = (2*(h[i-1] + h[i])) - ((h[i-1]**2)/u[i-1])
    
    
for i in range(N):    
    v[1] = b[1] - b[0]
    v[i] = (b[i] - b[i-1]) - ((h[i-1]*v[i-1])/u[i-1]) 
    
for i in range(data.shape[0]-1, 0, -1):
    z[i] = (v[i-1] - h[i-1] * z[i+1]) / u[i-1]

z[0] = 0
    

A = np.zeros(data.shape[0]+1)
B = np.zeros(data.shape[0]+1)
C = np.zeros(data.shape[0]+1)

for i in range(data.shape[0]-1):
    A[i] = (z[i+1]-z[i])/(6*h[i])
    B[i] = z[i]/2
    C[i] = (-((h[i]*(z[i+1]+2*z[i]))/6) + (data[i+1,1]-data[i,1])/h[i])


Si = np.zeros(x.shape)
for i in range(data.shape[0]-1):
    Si += (data[i,1] + (x-data[i,0])*(C[i] + (x-data[i,0])*(B[i] + ((x-data[i,0])*A[i])))) * ((x > data[i,0]) & (x <= data[i+1,0]))


S = np.zeros(x.shape)

for i in range(data.shape[0]-1):
    if i == 0:
        S += Si * (x <= data[i+1, 0])
    elif i == data.shape[0]-2:
        S += Si * (x > data[i, 0])
    else:
        S += Si * ((x > data[i, 0]) & (x <= data[i+1, 0]))


fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
    
axes.plot(data[:,0], data[:,1], 'ko', label = "data")
axes.plot(x, S, 'r', label = "S(x)")

axes.set_xlabel("x")
axes.set_ylabel("y")
axes.legend(loc=2)

plt.show()






