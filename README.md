# ASTR406
data analysis with astropy tools

The whole repo is a collection of assignments in UIUC ASTR 406. This course mainly focuses on galaxy evolution. First, there would be an introduction of what asstropy is, then the description of each homework and code display.

## Astropy

Astropy is a free and open-source Python package for astronomy research, first released in 2011. Developed and maintained by a global community, it provides tools for coordinate transformations, file reading, unit manipulation, and data modeling. Managed by a Coordination Committee, Astropy has over 280,000 lines of code and is supported by donations and institutions like the Space Telescope Science Institute. It is included in major Python distributions.


## Homework 1

### 1. Generate three figures that demonstrate the effect on the Schechter function for a range of α, M ∗, and n∗, varied separately. Plot the Schechter function with absolute magnitude on the x-axis, over a galaxy absolute magnitude range of −16 > M > −24. Choose parameters and plot parameters such that the effect of changing each parameter can be clearly seen. Be sure to label your axes and use a legend. Plot the results such that the x-axis shows luminosity increasing to the right.

The Schechter luminosity function ϕ provides an approximation of the abundance of galaxies in a luminosity interval [L + dL]. The luminosity function has units of a number density n per unit luminosity and is given by a power law with an exponential cut-off at high luminosity.

<img width="604" alt="Screenshot 2025-07-01 at 15 47 18" src="https://github.com/user-attachments/assets/e777f778-d7a4-47f4-9075-b24fdc611ee3" />

Note that because the magnitude system is logarithmic, the power law has logarithmic slope α+1
 This is why a Schechter function with α = −1 is said to be flat.


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


## Homework 2

In order to measure how bright these stars are, subtract off an estimate of the background
noise using np.nanmedian for each image. Print the value of the background noise for
each image.

```python
with fits.open('iaby02x7q_drc.fits') as hdu:
f606 = hdu[1].data
f606_hdr = hdu[1].header
with fits.open('iaby02x4q_drc.fits') as hdu:
f336 = hdu[1].data
f336_hdr = hdu[1].header
[3]: print(np.nanmedian(f606))
print(np.nanmedian(f336))
f606 -= np.nanmedian(f606)
f336 -= np.nanmedian(f336)


0.16780604
0.0048121833
```
Use DAOStarFinder from the photutils package to generate a list of sources the F606
filter image (filename iaby02x7q drc.fits). Print a subset of your source table.

```python
image = f606
bkg_sigma = np.nanstd(image)
daofind = DAOStarFinder(fwhm=2.0, threshold=2.0 * bkg_sigma)
sources = daofind(image)
for col in sources.colnames:
sources[col].info.format = '%.8g' # for consistent table output (optional)
print(sources)

id xcentroid ycentroid sharpness roundness1 roundness2 npix sky peak
flux mag
---- --------- --------- ---------- ------------- ------------- ---- ---
--------- --------- ------------
1 3.0406029 18.694236 0.72781123 0.52785188 0.4277751 25 0
183.37491 1.1115466 -0.11481922
2 226.13532 18.881369 0.76138002 0.20785956 0.18362588 25 0
263.76447 1.7093722 -0.58209157
3 151.70402 29.103142 0.66134122 0.090862155 0.065903529 25 0
188.86624 1.2204299 -0.2162821
4 95.540585 32.573552 0.74137947 0.69706047 0.19501715 25 0
1325.764 6.2833886 -1.9954848
5 413.48356 38.77244 0.68142404 0.35486788 0.046435209 25 0
189.00046 1.1868695 -0.18600743
6 194.44305 43.49156 0.79763919 0.77449512 0.53015776 25 0
1158.8774 5.2824454 -1.8070875
7 131.0649 47.776011 0.72162032 0.28290766 0.23585871 25 0
297.05051 1.9370632 -0.71785949
8 344.51341 48.94181 0.67204465 0.16990086 0.14771862 25 0
222.64236 1.4626118 -0.41282268
9 277.00845 49.778098 0.74608099 0.17947522 0.11080877 25 0
212.85767 1.3805077 -0.35009709
10 493.85381 52.179783 0.7147804 -0.011016716 -0.066366939 25 0
2
```

