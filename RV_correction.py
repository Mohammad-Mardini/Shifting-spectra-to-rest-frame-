#Import necessary packages
from __future__ import print_function, division
from matplotlib import gridspec
import matplotlib.pylab as plt
from PyAstronomy import pyasl
from pylab import *
import numpy as np

# Read template
Template = np.loadtxt("J2114-0616_shifted.txt" , delimiter="  ", skiprows=2)
# Read data, which are not RV corrected
Spectrum = np.loadtxt("J1645+4357.txt", delimiter="  ", skiprows=2)


# Locate portion of the spectrum with strong lines to be used in cross-correlation, here we used Halpha at 6562.797
ind = np.where((6553.8 < Spectrum[:,0]) & (Spectrum[:,0]< 6572.0)) 
indx= np.where((6540.0 < Template[:,0]) & (Template[:,0]< 6590.8)) 
waveTemplate = Template[:,0][indx[0]] 
fluxTemplate = Template[:,1][indx[0]] 
waveSpectrum = Spectrum[:,0][ind[0]] 
fluxSpectrum = Spectrum[:,1][ind[0]] 




# Check the selected lines
plt.plot(waveTemplate, fluxTemplate, 'b.-')
plt.plot(waveSpectrum, fluxSpectrum, 'r.-')
plt.show()



# Carry out the cross-correlation.
# The RV-range is -300 - +300 km/s in steps of 0.1 km/s.
# The first and last 50 points of the data are skipped.
rv, cc = pyasl.crosscorrRV(waveSpectrum, fluxSpectrum,waveTemplate, fluxTemplate, -300., 300., 0.10, skipedge=50)


# Find the index of maximum cross-correlation function
maxind = np.argmax(cc)

print("RV = ", rv[maxind], " km/s")
if rv[maxind] > 0.0:
  print("  A red-shift with respect to the template")
else:
  print("  A blue-shift with respect to the template")

# Check the cross-correlation peak
plt.plot(rv, cc, 'b-')
plt.plot(rv[maxind], cc[maxind], 'ro')
plt.show()

# Shift wavelength to rest frame
indxx= np.where((4028.9 < Spectrum[:,0]) & (Spectrum[:,0]< 6737.0)) 
Wave_rest=Spectrum[:,0][indxx[0]]  * (1- (rv[maxind]/299792))
flux_rest=Spectrum[:,1][indxx[0]]

# Save the RV-corrected spectrum
np.savetxt('RV-corrected spectrum.txt', (np.c_[Wave_rest, flux_rest]), delimiter='   ', fmt='%1.4f',newline='\n ')
