# ASTR406
data analysis with astropy tools

The whole repo is a collection of assignments in UIUC ASTR 406. 

## Homework 1

### 1. Generate three figures that demonstrate the effect on the Schechter function for a range of α, M ∗, and n∗, varied separately. Plot the Schechter function with absolute magnitude on the x-axis, over a galaxy absolute magnitude range of −16 > M > −24. Choose parameters and plot parameters such that the effect of changing each parameter can be clearly seen. Be sure to label your axes and use a legend. Plot the results such that the x-axis shows luminosity increasing to the right.

```python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit

# Schechter function parameters
def sf(m_arr, mstar, alpha, nstar):
phi = np.log(10)/2.5*nstar*(10**(-0.4*(m_arr-mstar)))**(alpha+1)*np.
↪exp(-1*10**(-0.4*(m_arr-mstar)))
return(np.log10(phi))
# create an array of luminosities over which to plot
m_arr = np.arange(-24,-16,0.1)
#vary alpha and plot the results
alpha = np.arange(-3,3,0.5)
colors = plt.cm.jet(np.linspace(0,1,12))
mstar = -20
n= 0.02
for j,a in enumerate(alpha):
plt.plot(m_arr,sf(m_arr,mstar,a,n),label=r'$\alpha=$'+str(a),
color=colors[j])
plt.legend()
plt.ylabel(r'log $\phi(M)$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.xlabel('M (mag)')
plt.xlim(-16,-24)
plt.show()

```
<img width="571" alt="Screenshot 2025-06-29 at 21 41 28" src="https://github.com/user-attachments/assets/43bad1e0-9609-4810-b92a-e566d037e72d" />

```python
alpha = -1
n = np.arange(-4,0,0.5)
colors = plt.cm.jet(np.linspace(0,1,8))
mstar = -20
for j,nn in enumerate(n):
plt.plot(m_arr,sf(m_arr,mstar,alpha,10**nn),
label=r'$\log n^* /Mpc^{-3}=$'+str(nn),
color=colors[j])
plt.legend()
plt.ylabel(r'log $\phi(M)$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.xlabel('M (mag)')
plt.xlim(-16,-24)
plt.show()
# vary M star and plot the results
alpha = -1
n = 0.02
colors = plt.cm.jet(np.linspace(0,1,10))
```
<img width="577" alt="Screenshot 2025-06-29 at 22 12 12" src="https://github.com/user-attachments/assets/52a462fa-3c9e-40fc-9108-931f8dc25921" />

```python
#vary M star and plot the results
alpha = -1
n = 0.02
colors = plt.cm.jet(np.linspace(0,1,10))
mstar = np.arange(-18,-23,-0.5)
for j,l in enumerate(mstar):
plt.plot(m_arr,sf(m_arr,l,alpha,n),
label=r'$M^*=$ '+str(round(l,2)),
color=colors[j])
plt.legend()
plt.ylabel(r'log $\phi(M)$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.xlabel('M (mag)')
plt.xlim(-16,-24)
plt.ylim(-30,5)
plt.show()
```
<img width="570" alt="Screenshot 2025-06-29 at 22 13 36" src="https://github.com/user-attachments/assets/5d913d43-3902-4993-958d-1dbef7b39b0c" />






### 2. The file schechter data.csv contains some simulated data for measurements of log φ (in units of Mpc−3 mag−1) at a range of absolute magnitudes. Read in this data and plot the results such that the x-axis shows luminosity increasing to the right.

### 3.Estimate* the parameters α, M ∗, or n∗ for the best-fitting Schechter function associated with this dataset. Include a plot comparing your best-fit model to the data. (*Estimate here means, you don’t need to report uncertainties on your final values.)

```python
plt.plot(tab['mag'],tab['logphi'],'ko')
plt.ylabel(r'log $\phi(M)$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.xlabel('M (mag)')
plt.xlim(-16,-24)
plt.show()

M_valid = np.array(data['mag'])
logphi_valid = np.array(data['logphi'])
def f_log(M, alpha, M_star, n_star):
return np.log10(f(alpha, n_star, M_star, M))
initial_guess = [-1.0,
-21.0, 1.0]
params, _ = sp.optimize.curve_fit(f_log, M_valid, logphi_valid,
p0=initial_guess)
alpha_fit, M_star_fit, n_star_fit = params
M_fit = np.linspace(-24,
-16, 400)
logphi_fit = f_log(M_fit, alpha_fit, M_star_fit, n_star_fit)
plt.figure(figsize=(8, 6))
plt.plot(M_valid, logphi_valid, 'bo', label='Simulated Data')
plt.plot(M_fit, logphi_fit, 'r-', label=f'Best-Fit Model
(alpha={alpha_fit:.2f}, M*={M_star_fit:.2f}, n*={n_star_fit:.2f})')
plt.gca().invert_xaxis()
plt.xlabel('Absolute Magnitude (M)')
plt.ylabel('log(φ) (Mpc$^{-3}$ mag$^{-1}$)')
plt.title('Best-Fit Schechter Function using User\'s Code')
plt.legend()
```
<img width="874" alt="Screenshot 2025-06-29 at 22 11 00" src="https://github.com/user-attachments/assets/cdb0af12-349e-4f41-9107-2ea34f4694a0" />


