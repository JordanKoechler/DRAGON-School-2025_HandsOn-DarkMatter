import healpy
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.ticker as tick
import matplotlib.colors as colors
from matplotlib import cm
from astropy.io import fits as pyfits
import numpy as np

### MODELS
print('Opening Hermes maps')
HermesRuns = 'HERMES_output/'

GeV_units = 1/(1.6021766339999998e-10)  # From Joules to GeV
SI_units = 1e4 # m2 to cm2

Comp_IC_isrf_ = HermesRuns + 'hermes_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann-IC_isrf_nside64.fits.gz'
Comp_IC_cmb_  = HermesRuns + 'hermes_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann-IC_cmb_nside64.fits.gz'

hdulisrf = pyfits.open(Comp_IC_isrf_)
hdulcmb  = pyfits.open(Comp_IC_cmb_)

Hermes_IC_isrf, Hermes_IC_cmb, E_val = [], [], []
for iE in range(len(hdulisrf))[1:]:
    E_val += [hdulisrf[iE].header['ENERGY']*GeV_units]
    Hermes_IC_isrf += [np.array(hdulisrf[iE].data.astype(float)/SI_units)]
    Hermes_IC_cmb  += [np.array(hdulcmb [iE].data.astype(float)/SI_units)]

nside = hdulisrf[1].header['NSIDE']

hdulisrf.close()
hdulcmb.close()

MapISRF = np.array(Hermes_IC_isrf)
MapCMB  = np.array(Hermes_IC_cmb)

iE = 2
print(E_val[iE]*1e3, 'MeV')

plt.rcParams.update({'font.size':14.5})
healpy.mollview(MapISRF[iE], title = 'amb photons: ISRF', unit=r'$cm^{-2}s^{-1}GeV^{-1}sr^{-1}$', norm = "hist", cmap = plt.cm.rainbow) 
healpy.graticule()
plt.tight_layout(pad=0.5)
plt.savefig('IC_flux_map_ISRF.pdf')
plt.show()

plt.rcParams.update({'font.size':14.5})
healpy.mollview(MapCMB[iE], title = 'amb photons: CMB', unit=r'$cm^{-2}s^{-1}GeV^{-1}sr^{-1}$', norm = "hist", cmap = plt.cm.rainbow) 
healpy.graticule()
plt.tight_layout(pad=0.5)
plt.savefig('IC_flux_map_CMB.pdf')
plt.show()