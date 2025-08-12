# stationary solution for infinite well

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, electron_mass


m = electron_mass

P = 60000  # number of timesteps
N = 300 # number of spacesteps
L = 1e-8 # size of the well
del_t = 1e-18 # timestep
x = np.linspace(0, L, N)  # space array
del_x = x[1] - x[0]
x_0 = L/2  # center of the well
t = np.arange(0, P*del_t, del_t)  # time array


n = 10  # mode number

psi = np. zeros((P, N), dtype=np.complex_)  # wave function matrix time x space
psi[0, :] = np.sqrt(2/L)*np.sin(n*np.pi*x/L)  # infinite well eigenfunction

V = np.zeros(N)
V[0] = 1e40  # approximation of infinite potential
V[N-1] = 1e40

# calculating the time needed to reach phase pi/2
E = (hbar*n*np.pi/L)**2/(2*m)
print('E = ', E)
t = np.pi*hbar/(2*E)  # imaginary part for phase pi/2
print('pi/2 timestep: ', int(t/del_t))


# defininig the matrices 
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

# solving the matrix equation for consecutive timesteps
for p in range(P-1):
    psi[p+1] = np.linalg.solve(A, B @ psi[p])
    if p%50==0:
        print(p)
        
        
        
# plots - modulus and imaginary part for different times 

plt.figure()
plt.title('abs(psi)')
plt.plot(x, abs(psi[0]), label='timestep = 0')
plt.plot(x, abs(psi[5000]), label='timestep = 5000')
plt.plot(x, abs(psi[10000]), label='timestep = 10000')
plt.plot(x, abs(psi[30000]), label='timestep = 30000')
plt.plot(x, abs(psi[50000]), label='timestep = 50000')
plt.legend()
plt.xlabel('Position')
plt.show()
        
plt.figure()
plt.title('psi.imag')
plt.plot(x, psi[0].imag, label='timestep = 0')
plt.plot(x, psi[5000].imag, label='timestep = 5000')
plt.plot(x, psi[10000].imag, label='timestep = 10000')
plt.plot(x, psi[30000].imag, label='timestep = 30000')
plt.plot(x, psi[50000].imag, label='timestep = 50000')
plt.xlabel('Position')
plt.legend()
plt.show()

plt.figure()
plt.title('psi.real')
plt.plot(x, psi[0].real, label='timestep = 0')
plt.plot(x, psi[5000].real, label='timestep = 5000')
plt.plot(x, psi[10000].real, label='timestep = 10000')
plt.plot(x, psi[30000].real, label='timestep = 30000')
plt.plot(x, psi[50000].real, label='timestep = 50000')
plt.xlabel('Position')
plt.legend()
plt.show()


# plots - imaginary part for phase pi/2

plt.figure()
plt.title('Imaginary part for phase pi/2')
plt.plot(x, abs(psi[int(t/del_t)]), color='orange')
plt.plot(x, -abs(psi[int(t/del_t)]), color='orange')
plt.plot(x, psi[int(t/del_t)].imag, label = 'pi/2 psi.imag')
plt.xlabel('Position')
plt.legend()
plt.show()

plt.figure()
plt.title('Imaginary part for phase 21*pi/2')
plt.plot(x, abs(psi[int((4*5+1)*t/del_t)]), color='orange')
plt.plot(x, -abs(psi[int((4*5+1)*t/del_t)]), color='orange')
plt.plot(x, psi[int((4*5+1)*t/del_t)].imag, label = 'pi/2 psi.imag')
plt.xlabel('Position')
plt.legend()
plt.show()

