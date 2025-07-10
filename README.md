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
Make a new plot of the apertures overlaid on the F336 image to verify whether these
images are aligned, or if you will need to apply a small offset. If the final goal is to make
a color-magnitude diagram, why can’t we just use DAOStarFinder on the second image?

```python
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
# define apertures
apertures = CircularAperture(positions, r=8.0) #again, full credit for␣
↪something reasonable here
# get aperture photometry
phot_table = aperture_photometry(image, apertures)
# show preview of table for inspection
for col in phot_table.colnames:
phot_table[col].info.format = '%.8g' # for consistent table output
print(phot_table)
id xcenter ycenter aperture_sum
pix pix
---- --------- --------- ------------
1 3.0406029 18.694236 2255.0242
2 226.13532 18.881369 nan
3 151.70402 29.103142 1906.4741
4 95.540585 32.573552 19863.044
5 413.48356 38.77244 2185.7362
6 194.44305 43.49156 24007.486
7 131.0649 47.776011 2476.4301
8 344.51341 48.94181 2281.5525
9 277.00845 49.778098 1609.5245
10 493.85381 52.179783 1228.201
11 56.133717 54.155248 2395.3979
12 906.3341 55.735058 nan
13 510.52482 58.938272 2661.3064

```
```python
plt.imshow(image, cmap='gray_r', origin='lower',vmin=0,vmax=50)
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.colorbar()
plt.xlim(500,1000)
plt.ylim(500,1000)
plt.show()
```

<img width="487" alt="Screenshot 2025-07-03 at 17 24 59" src="https://github.com/user-attachments/assets/c3520316-6710-4933-9b41-36c96c21ad09" />

Use the aperture sums from your tables above, as well as your answer to 1c above to
calculate the F336-F606 color and F606 apparent magnitude for your list of stars. Make
a plot of these values following the typical conventions of color-magnitude diagrams
(color becoming redder to the right on the x-axis, apparent magnitude decreasing as the
y-axis increases).
Label the main sequence, horizontal branch, and red giant stars on your HR diagram.

```python
color = -2.5*np.log10(phot_table2['aperture_sum']/phot_table['aperture_sum'])
flam = phot_table['aperture_sum']*f606_hdr['PHOTFLAM']
mag = f606_hdr['PHOTZPT']-2.5*np.log10(flam)
plt.plot(color, mag,linestyle='None',alpha=0.5,marker='.',markersize=1)
9
plt.xlim(1,5)
plt.ylim(20,12)
plt.xlabel('F336-F606')
plt.ylabel('m (F606)')
plt.text(3,19.5,'Main Sequence')
plt.text(2.5,15,'Horizontal Branch')
plt.text(4,14,'Red Giants')
plt.show()
```
<img width="576" alt="Screenshot 2025-07-09 at 16 50 46" src="https://github.com/user-attachments/assets/d1910d05-d813-4f93-94b1-5a87ddf24f21" />

## Homework 3

### 1. (a) Determine the maximum absolute magnitude MV,max of a star that can be observed with the unaided eye, as a function of its distance d in parsecs. Graph this function on a plot with MV on the y-axis and log d on the x-axis, and overplot (for comparison) points corresponding the 25 brightest naked-eye stars (excluding the Sun) from the file provided, which includes the names of the stars, the distance (in lightyears), and apparent and absolute magnitude.

```python
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import constants as const
[3]: # plot this line alongside the stars from the file
dat = ascii.read('brightest.csv')
x = np.logspace(0,3,20)
y = 10 - 5*np.log10(x)
plt.semilogx(x,y,label='obs limit')
dpc = dat['Distance']*u.lyr.to(u.pc)
plt.semilogx(dpc,dat['AbsMag'],linestyle='None',marker='o',label='25 brightest')
plt.gca().invert_yaxis()
plt.xlabel('Distance [pc]')
plt.ylabel('$M_V$')
plt.legend()
plt.show()
```
<img width="581" alt="Screenshot 2025-07-09 at 17 45 02" src="https://github.com/user-attachments/assets/4b1d17d6-ef4b-44f1-8441-728abe87a06b" />

### 2. 
The more distant of the bright stars are more luminous for two reasons: 1) they have to be more
luminous to produce the same flux as the nearby bright stars and 2) these stars are extremely rare,
so they tend to be very far apart and would be unlikely to be found at close distances.

### 3.
Create a model longitude-velocity diagram for the Galaxy by assuming this model galaxy
is a series of thin rings at radii R ranging from 1 to 15 kpc spaced by 1 kpc. Plot the
curves
vr(l) = R0[Ω(R) − Ω(R0)] sin l
for each ring on one plot to construct your model diagram. Be sure you limit your range
of l if needed using your results from parts a and b (Hint: while you could iterate over
each ring azimuthally, it is much easier to iterate over l). Use colors and/or line styles
to indicate the radius of each galaxy component and whether it is interior to or exterior
to the Sun.

```python
plt.figure(figsize=(10,6))
R0=8 # kpc
V0=220 # km/s
radii = np.linspace(1,15, num=15) # kpc
color = plt.cm.rainbow(np.linspace(0, 1, 15))
for j, R in enumerate(radii):
omega = V0/R
omega0 = V0/R0
if R > R0:
long = np.linspace(-180, 180, num=100)
vr = R0 * (omega-omega0) * np.sin(np.radians(long))
plt.plot(long, vr,linestyle=':',c=color[j],label=str(R)+' kpc')
else: # inner rings have narrower longitude range
lmax = np.degrees(np.arcsin(R/R0))
long = np.linspace(-lmax, lmax, num=100)
vr = R0 * (omega-omega0) * np.sin(np.radians(long))
plt.plot(long, vr,c=color[j],label=str(R)+' kpc')
plt.gca().invert_xaxis() #reversed to match S&G, but not needed for full credit
plt.xlabel('Galactic Longitude [deg]', fontsize=14)
plt.ylabel('$v_R$ [km/s]', fontsize=14)
plt.xlim(180,-180)
plt.legend()
plt.show()
```
<img width="577" alt="Screenshot 2025-07-09 at 17 58 26" src="https://github.com/user-attachments/assets/9b1bc9e7-32d5-4446-ba8c-c0bfe324c9d0" />

## Homework 4

### 1. For a star moving on a radial orbit in this potential (i.e., in a straight line through the center), it will oscillate harmonically with a period P ∼ 3tf f , where the free-fall timescale is tf f = (Gρ)−1/2. Using the circular velocity of the Sun’s orbit of v = 200 km/s, find the mass (in units of solar masses) and average density inside the sphere (in units of solar masses per cubic parsec) with radius r = 8 kpc. What is the free-fall time for gas at this density, in units of years? Despite the fact that we are assuming aspherical constant density, the value you obtain is close to the free-fall time we observe in the galaxy, which limits the timescales of processes such as star formation bursts, because gravitational forces cannot move material any faster through the galaxy.

```python
import astropy.units as u
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt
v = 200*u.km/u.s
r = 8*u.kpc
G = const.G
rho = v**2*3/(4*np.pi*G*r**2)
print('Average density is ',rho.to(u.Msun/u.pc**3))
mass = rho*4*np.pi*r**3/3
print('Mass is {:.2e}'.format(mass.to(u.Msun)))
t_ff = (G*rho)**(-1/2)
print('Free fall time is {:.2e}'.format(t_ff.to(u.year)))
Average density is 0.03469207840572885 solMass / pc3


#Mass is 7.44e+10 solMass
#Free fall time is 8.00e+07 yr
```
### Derive an expression for the the rotation curve vc(r) in terms of a and M∞. For a=3 kpc and M∞ = 2 × 1011 M , make a plot of the rotation curve. Calculate the maximum value for vc and the radius where vc is at its maximum.

```python
mtot = 2.e+11 * const.M_sun
r = np.linspace(0,20,num=100) * u.kpc
a = 3 * u.kpc
vcirc = np.sqrt((const.G*mtot*r/(r+a)**2)).to(u.km/u.s)
# Plot the rotation curve
plt.plot(r, vcirc)
plt.xlabel('Radius (kpc)', fontsize=12)
plt.ylabel('Circular velocity (km/s)', fontsize=12)
plt.show()
print('The maximum value of vcirc is {0} at {1}'.format(np.max(vcirc),r[np.
↪argmax(vcirc)]))
```
<img width="549" alt="Screenshot 2025-07-09 at 18 26 33" src="https://github.com/user-attachments/assets/71d15c2e-c873-45df-9aa3-b224587b4f32" />

### Graph Ω, Ω − κ/2, and Ω + κ/2 in units of Ω0 as functions of r/a0, and show that the maximum pattern speed for which an inner m=2 Lindblad resonance is possible is given by Ωp ≈ Ω0/8.

```python
def omega(r, a_0, omega_0):
return omega_0/((r/a_0)**2+1)**0.75
def kappa(r, a_0, omega_0):
return omega_0*((r/a_0)**2+4)**0.5/((r/a_0)**2+1)**1.25
def dvdr(r, a_0, omega_0):
return omega_0/((r/a_0)**2+1)**0.75 - 3*r**2*omega_0/(
2*a_0**2*((r/a_0)**2+1)**1.75)
# Orbital frequencies vs. radius
rad = np.linspace(0.1,10,num=100) * u.kpc
a_0 = 1 * u.kpc
omega_min_kappa = omega(rad,a_0,1)-0.5*kappa(rad,a_0,1)
omega_plu_kappa = omega(rad,a_0,1)+0.5*kappa(rad,a_0,1)
plt.plot(rad/a_0,omega(rad,a_0,1),linestyle='-',label=r'$\Omega$')
plt.plot(rad/a_0,omega_min_kappa,linestyle='--',label=r'$\Omega-\kappa/2$')
plt.plot(rad/a_0,omega_plu_kappa,linestyle='--',label=r'$\Omega+\kappa/2$')
plt.xlabel('r / $a_0$', size=14)
plt.ylabel(r'frequency / $\Omega_0$', size=14)
plt.legend(loc='upper right', fontsize=12)
plt.show()
print("The maximum frequency for an ILR is {0}".format(np.amax(
omega_min_kappa)))
```

<img width="568" alt="Screenshot 2025-07-10 at 10 03 03" src="https://github.com/user-attachments/assets/ae98d041-f21e-4ba3-9c01-1bfd9eca52df" />

### Consider a bar (m = 2) with pattern speed Ωp = Ω0/10. Where would the Outer Lindblad Resonance be found, in units of r/a0?

```python
rad[np.argmin(np.abs(0.1-omega_plu_kappa))].value
[6]: 6.1
```

### Make a graph showing both the potential and the effective potential Φeff (r) for an orbit in a Plummer potential with a0 = 1 kpc and M = 1010M , with guiding center rg = 8 kpc. You will first need to determine the angular momentum L of the guiding center orbit. Confirm numerically, using the np.argmin function, that the minimum of the effective potential actually occurs at r = rg.


```python
def phi(mtot, r, a_0):
return (-const.G*mtot/np.sqrt(r**2+a_0**2)).to(u.km**2/u.s**2)
a_0 = 1*u.kpc
mtot = 1e10*u.solMass
rg=8*u.kpc
def phi_eff(mtot, r, a_0, rg):
omega_0 = np.sqrt(const.G*mtot/a_0**3)
Lz = rg**2*omega(rg, a_0, omega_0)
return phi(mtot, r, a_0) + Lz**2/(2*r**2)
rad = np.linspace(2,16,num=100) * u.kpc
# Plot the actual galaxy potential
plt.plot(rad,phi(mtot,rad,a_0),linestyle='-',label='$\Phi$')
plt.plot(rad,phi_eff(mtot,rad,a_0,rg),linestyle='-',label='$\Phi_{eff}$')
plt.legend()
plt.xlabel('Radius (kpc)', size=14)
plt.ylabel('Potential $(km/s)^2$', size=14)
plt.show()
rmin = rad[np.argmin(phi_eff(mtot,rad,a_0,rg))]
print('The minimum of the effective potential is at {}'.format(rmin))
```

<img width="593" alt="Screenshot 2025-07-10 at 10 05 34" src="https://github.com/user-attachments/assets/2d802f12-3a31-4df7-b9fb-e2c2b631a97c" />







