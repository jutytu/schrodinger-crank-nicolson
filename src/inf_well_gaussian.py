# infinite well - gaussian packet

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, electron_mass


m = electron_mass
sigma = 0.5*1e-9
kappa = 5e10

P = 1000000 # num of time steps
N = 300 # num of spatial steps
L = 1e-8 # total length of the infinite well
del_t = 1e-19 # time step size
x = np.linspace(0, L, N) # space
del_x = x[1] - x[0]
x_0 = L/2 # well center
t = np.arange(0, P*del_t, del_t) # time


psi = np. zeros((P, N), dtype=np.complex_) # wave function of dimensions time x space
psi[0, :] = np.exp(-(x-x_0)**2/(2*sigma**2))*np.exp(1j*kappa*x)/(sigma*np.sqrt(2*np.pi)) # defining the function at time=0

V = np.zeros(N)
V[0] = 1e40  # using a big value as infinity approximation for the potential
V[N-1] = 1e40


# defining the A, B matrices
A = np.zeros((N, N) ,dtype=np.complex_)
B = np.zeros((N, N), dtype=np.complex_)

for i in range(N):
    if i==0:
        A[i, i] = 1 + 1j*del_t*hbar/(2*m*del_x**2) + 1j*del_t*V[i]/(2*hbar)
        A[i, i+1] = -1j*del_t*hbar/(4*m*del_x**2)
        B[i, i] = 1 - 1j*del_t*hbar/(2*m*del_x**2) - 1j*del_t*V[i]/(2*hbar)
        B[i, i+1] = 1j*del_t*hbar/(4*m*del_x**2)
    elif i==N-1:
        A[i, i] = 1 + 1j*del_t*hbar/(2*m*del_x**2) + 1j*del_t*V[i]/(2*hbar)
        A[i, i-1] = -1j*del_t*hbar/(4*m*del_x**2)
        B[i, i] = 1 - 1j*del_t*hbar/(2*m*del_x**2) - 1j*del_t*V[i]/(2*hbar)
        B[i, i-1] = 1j*del_t*hbar/(4*m*del_x**2)
    else:
        A[i, i] = 1 + 1j*del_t*hbar/(2*m*del_x**2) + 1j*del_t*V[i]/(2*hbar)
        A[i, i+1] = -1j*del_t*hbar/(4*m*del_x**2)
        A[i, i-1] = -1j*del_t*hbar/(4*m*del_x**2)
        B[i, i] = 1 - 1j*del_t*hbar/(2*m*del_x**2) - 1j*del_t*V[i]/(2*hbar)
        B[i, i+1] = 1j*del_t*hbar/(4*m*del_x**2)
        B[i, i-1] = 1j*del_t*hbar/(4*m*del_x**2)

# solving the matrix equation for consecutive time steps
for p in range(P-1):
    psi[p+1] = np.linalg.solve(A, B @ psi[p])
    if p%50==0:
        print(p)
        
        
# plots
plt.figure()
plt.subplot(5, 1, 1)
plt.plot(x, abs(psi[0]))
plt.plot(x, psi[0].imag)
plt.plot(x, psi[0].real)
plt.subplot(5, 1, 2)
plt.plot(x, abs(psi[15000]))
plt.plot(x, psi[15000].imag)
plt.plot(x, psi[15000].real)
plt.subplot(5, 1, 3)
plt.plot(x, abs(psi[30000]))
plt.plot(x, psi[30000].imag)
plt.plot(x, psi[30000].real)
plt.subplot(5, 1, 4)
plt.plot(x, abs(psi[60000]))
plt.plot(x, psi[60000].imag)
plt.plot(x, psi[60000].real)
plt.subplot(5, 1, 5)
plt.plot(x, abs(psi[80000]))
plt.plot(x, psi[80000].imag)
plt.plot(x, psi[80000].real)
plt.xlabel('Position')
plt.show()
        
        
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(x, abs(psi[500000]))
plt.plot(x, psi[500000].imag)
plt.plot(x, psi[500000].real)
plt.subplot(2, 1, 2)
plt.plot(x, abs(psi[900000]))
plt.plot(x, psi[900000].imag)
plt.plot(x, psi[900000].real)
plt.xlabel('Position')
plt.show()  
        
        
        