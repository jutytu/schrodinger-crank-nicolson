# Schrödinger equation solved using the Crank-Nicolson method

This repository contains my part of the group project for a university Computational Physics course.  
The focus is on numerical solutions of the time-dependent Schrödinger equation using the Crank–Nicolson method, with applications to:

- Infinite potential well (Gaussian packets & eigenstates)
- Step potential (transmission and reflection probabilities)

![](figures/slider.gif)


## Repository structure

```bash
schrodinger-crank-nicolson/
│
├── README.md
├── LICENSE
├── requirements.txt       
│
├── src/                   # Source code
│   ├── inf_well_gaussian.py
│   ├── inf_well_stationary.py
│   ├── step_potential_gauss.py
│
├── figures/               # All generated plots & gifs
│
├── docs/                  # Additional documentation
    └── compphys_tekst.pdf

```


## 1. Theory

### 1.1 Schrödinger Equation

We study the time-dependent Schrödinger equation in 1D:

\[
i \hbar \frac{\partial \psi(x,t)}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \psi(x,t)}{\partial x^2} + V(x) \psi(x,t)
\]

---

### 1.2 Crank–Nicolson Method

The Crank–Nicolson scheme is obtained by averaging the explicit and implicit finite-difference schemes:

\[
\frac{\psi^{n+1} - \psi^n}{\Delta t} = \frac{1}{2} \hat{H} (\psi^{n+1} + \psi^n)
\]

which leads to the matrix equation:

\[
A \psi^{n+1} = B \psi^n
\]

- **Unconditionally stable**
- **Second-order accurate in time and space**
- **Low numerical diffusion**

---

### 1.3 Infinite Potential Well – Gaussian Packet

Initial wave function:

\[
\psi(x,0) = e^{-\frac{(x-x_0)^2}{2\sigma^2}} e^{ikx}
\]

Observations:
- Wave packet bounces between walls
- Modulus conserved for short times
- Long-term: numerical diffusion causes packet broadening

---

### 1.4 Infinite Well – Eigenstates

Eigenfunctions:

\[
\psi_n(x) = \sqrt{\frac{2}{L}} \sin \left( \frac{n \pi x}{L} \right)
\]
\[
E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}
\]

We verify:
- Modulus remains constant
- Real & imaginary parts oscillate with expected frequency:
\[
\omega_n = \frac{E_n}{\hbar}
\]

---

### 1.5 Step Potential

Potential:
\[
V(x) =
\begin{cases}
0 & x < x_0 \\
V_0 & x \ge x_0
\end{cases}
\]

Classical transmission and reflection:
For \(E > V_0\):
\[
T = \frac{4k_1 k_2}{(k_1+k_2)^2}, \quad R = \frac{(k_1-k_2)^2}{(k_1+k_2)^2}
\]
with \( k_i = \sqrt{\frac{2m(E-V_i)}{\hbar^2}} \).

---

## 2. Results & Plots

All results are stored in `figures/`:

- **Infinite Well Gaussian Packet**
  - `inf_gauss_8k_steps.png`, `inf_gauss_50k_steps.png`, `inf_gauss_500k_steps.png`
- **Eigenstates**
  - `inf_eigen_real.png`, `inf_eigen_imag.png`, etc.
- **Step Potential**
  - `step_pot_ref1.png`, `step_pot_ref2.png`, `step_pot_trans1.png`, `step_pot_trans2.png`, `step_pot_theory.png`

Key observations:
- Gaussian packet broadens over time due to numerical effects
- Eigenstate oscillations match theoretical predictions
- Step potential reflection/transmission agree well with theory for small time steps

---
