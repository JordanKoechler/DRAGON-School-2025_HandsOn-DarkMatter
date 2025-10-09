import healpy
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.ticker as tick
import matplotlib.colors as colors
from matplotlib import cm
import astropy.io.fits as pyfits
import numpy as np

eBand1 = np.loadtxt('./data/INTEGRAL_27-49keV.txt',    skiprows=2)
eBand2 = np.loadtxt('./data/INTEGRAL_49-90keV.txt',    skiprows=2)
eBand3 = np.loadtxt('./data/INTEGRAL_100-200keV.txt',  skiprows=2)
eBand4 = np.loadtxt('./data/INTEGRAL_200-600keV.txt',  skiprows=2)
eBand5 = np.loadtxt('./data/INTEGRAL_600-1800keV.txt', skiprows=2)

eBands = [eBand1, eBand2, eBand3, eBand4, eBand5]
colours = ['red', 'orange', 'green', 'blue', 'purple']

sigmav = 1e-27

def pixbl(nside, bmin, bmax, lmin, lmax):
    npix=12*nside**2
    listpix = np.arange(npix)
    l, b = healpy.pixelfunc.pix2ang(nside, listpix, nest=False, lonlat=True)
    mask = []
    for i in np.arange(npix):
        if(l[i]>180):
            l[i]-=360
        if(bmax >= b[i] >= bmin and lmax >= l[i] >= lmin):
            mask.append(1)
        else:
            mask.append(0)
    return mask

npix = healpy.nside2npix(64)

### MODELS
print('Opening Hermes maps')
HermesRuns = 'HERMES_output/'

GeV_units = 1/(1.6021766339999998e-10)  # From Joules to GeV
SI_units = 1e4 # m2 to cm2

Comp_IC_isrf_ = HermesRuns + 'hermes_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann-IC_isrf_nside64.fits.gz'
Comp_IC_cmb_  = HermesRuns + 'hermes_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann-IC_cmb_nside64.fits.gz'

lmins = [-23.1, -23.1, -23.1, -23.1, -60]
lmaxs = [ 23.1,  23.1,  23.1,  23.1,  60]

Emins = [27, 49, 100, 200, 600]
Emaxs = [49, 90, 200, 600, 1800]

nbins = [21, 21, 21, 21, 15]

hdulisrf = pyfits.open(Comp_IC_isrf_)
hdulcmb  = pyfits.open(Comp_IC_cmb_)

plt.figure(figsize=(15, 8))

ax1 = plt.subplot2grid((2,6), (0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

plt.subplots_adjust(wspace=1., hspace=0.3)

axs = [ax1, ax2, ax3, ax4, ax5]

ax1.set_ylim(-0.025, 0.305)
ax2.set_ylim(-0.01, 0.1)
ax3.set_ylim(-0.01, 0.1)
ax4.set_ylim(-0.015, 0.14)
ax5.set_ylim(-0.002, 0.02)

ax1.set_xlim(-90, 90)
ax2.set_xlim(-90, 90)
ax3.set_xlim(-90, 90)
ax4.set_xlim(-90, 90)
ax5.set_xlim(-90, 90)

ax1.set_xlabel(r'Latitude $b$ [degree]')
ax2.set_xlabel(r'Latitude $b$ [degree]')
ax3.set_xlabel(r'Latitude $b$ [degree]')
ax4.set_xlabel(r'Latitude $b$ [degree]')
ax5.set_xlabel(r'Latitude $b$ [degree]')

ax1.set_ylabel(r'Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
ax2.set_ylabel(r'Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
ax3.set_ylabel(r'Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
ax4.set_ylabel(r'Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
ax5.set_ylabel(r'Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')

ax1.text(0.15, 0.9, r'27-49 keV', color='red', horizontalalignment='center',
        verticalalignment='center', transform=ax1.transAxes)

ax2.text(0.15, 0.9, r'49-90 keV', color='orange', horizontalalignment='center',
        verticalalignment='center', transform=ax2.transAxes)

ax3.text(0.15, 0.9, r'100-200 keV', color='green', horizontalalignment='center',
        verticalalignment='center', transform=ax3.transAxes)

ax4.text(0.15, 0.9, r'200-600 keV', color='blue', horizontalalignment='center',
        verticalalignment='center', transform=ax4.transAxes)

ax5.text(0.15, 0.9, r'600-1800 keV', color='purple', horizontalalignment='center',
        verticalalignment='center', transform=ax5.transAxes)

for j in range(5):
    
    DMflux = []

    for i in range(nbins[j]):

        bmin = eBands[j][i, 0]-eBands[j][i, 1]
        bmax = eBands[j][i, 0]+eBands[j][i, 1]

        ipix = pixbl(64, bmin, bmax, lmins[j], lmaxs[j])

        solidangle = 4 * np.pi * np.sum(ipix) / npix

        Hermes_IC_isrf, Hermes_IC_cmb, E_val = [], [], []
        for iE in range(len(hdulisrf))[1:]:
            Egamma = hdulisrf[iE].header['ENERGY']*GeV_units*1e6
            if Emins[j] < Egamma < Emaxs[j]:
                E_val += [Egamma]
                Hermes_IC_isrf += [np.array(hdulisrf[iE].data.astype(float)/SI_units)]
                Hermes_IC_cmb  += [np.array(hdulcmb [iE].data.astype(float)/SI_units)]

        FluxHermes_IC_isrf = np.array(Hermes_IC_isrf) @ ipix / npix * 4*np.pi / 1e-26 / 1e6 / solidangle / 0.07
        FluxHermes_IC_cmb  = np.array(Hermes_IC_cmb)  @ ipix / npix * 4*np.pi / 1e-26 / 1e6 / solidangle / 0.07
        FluxHermes_IC_tot  = FluxHermes_IC_isrf + FluxHermes_IC_cmb

        DMflux += [[eBands[j][i, 0], np.trapz(FluxHermes_IC_tot, E_val)]]

    DMflux = np.array(DMflux)

    axs[j].plot(DMflux[:, 0], sigmav*DMflux[:, 1], ls='-', color='k')
    axs[j].errorbar(eBands[j][:, 0], eBands[j][:, 2], xerr=eBands[j][:, 1], yerr=eBands[j][:, 3], ls='', color=colours[j])

hdulisrf.close()
hdulcmb.close()

plt.tight_layout(pad=0.5)
plt.savefig('plots/IC_flux_spec_INTEGRAL.pdf')
plt.show()