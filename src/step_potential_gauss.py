# step potential

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.constants import hbar, electron_mass
from matplotlib import animation


m = electron_mass
sigma = 1e-9
kappa = 5e10   
E = hbar**2*kappa**2/(2*m) # energy
print(E)

P = 5000 # number of time steps
N = 1000 # number of space steps
L = 1e-8 # total space size
del_t = 1e-18
x = np.linspace(0, L, N)
x_0 = x[int(N/4)]
del_x = x[1] - x[0]
t = np.arange(0, P*del_t, del_t)

# Sinusoidal wave function
# psi = np. zeros((P, N), dtype=np.complex_) # dimensions of time x space
# psi[0, 0:int(N/4)] = np.exp(1j*kappa*x[0:int(N/4)])

psi = np. zeros((P, N), dtype=np.complex_) # wave function of dimensions time x space
psi[0, :] = 1e-8*np.exp(-(x-x_0)**2/(2*sigma**2))*np.exp(1j*kappa*x)/(sigma*np.sqrt(2*np.pi)) # defining the function at time=0

# Step potential
V = np.zeros(N)
V[int(N/2):N] = 1e-17 # potential step starts in the middle of the space

# defining the matrices
A = np.zeros((N, N), dtype=np.complex_)
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

# calculating the wave function for consecutive time steps
for p in range(P-1):
    psi[p+1] = np.linalg.solve(A, B @ psi[p])
    if p%50==0:
        print(p) 


# reflection and transmission factors - only calculated for E > V
if E>V[N-1]:
    kappa2 = np.sqrt(2*m*(E-V[N-1])/hbar**2)
    A1 = np.max(abs(psi[0, 0:400]))
    A2 = np.max(abs(psi[int(0.5*P), 0:175]))
    B = np.max(abs(psi[int(0.5*P), 225:N]))
    R = (A2/A1)**2
    T = (B/A1)**2*kappa2/kappa
    print('Simulated R: ', R)
    print('Simulated T: ', T)
    Rtheo = ((kappa-kappa2)/(kappa+kappa2))**2
    Ttheo = (2*np.sqrt(kappa*kappa2)/(kappa+kappa2))**2
    print('Theoretical R: ', Rtheo)
    print('Theoretical T: ', Ttheo)



# Interactive plot

p = 0
def update_abs(val):
    p = int(slider.val)
    line1.set_ydata(psi[p].real + E*1e17)
    fig.canvas.draw_idle()

global_min = psi.real.min() + E*1e17
global_max = psi.real.max() + E*1e17

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)

line1, = ax.plot(range(N), psi[p].real + E*1e17, label="Simulated |psi|", color='b')
lineV, = ax.plot(range(N), V * 1e17, label="Step potential", color='r', linestyle='--')

ax.set_ylim(global_min * 1.1, global_max * 1.1)
ax.set_xlabel("x")
ax.set_ylabel("|psi|")
ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, 'Time Step', 0, P-1, valinit=0, valstep=1)
slider.on_changed(update_abs)

ax.legend(loc='upper right')
plt.show()


# def update(p):
#     line1.set_ydata(psi[p].real + E*1e17)
#     return line1,

# ani = animation.FuncAnimation(
#     fig, update, frames=range(P), interval=200, blit=True
# )
# ani.save("psi_slider_like_animation.gif", writer="pillow")





        
        
        
        
        
        