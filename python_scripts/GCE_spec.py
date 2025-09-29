import healpy
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.ticker as tick
import matplotlib.colors as colors
from matplotlib import cm
import astropy.io.fits as pyfits
import numpy as np

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1

sigmav = 1.3e-26

nside = 64

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

npix = healpy.nside2npix(nside)

### MODELS
print('Opening Hermes maps')
HermesRuns = 'HERMES_output/'

GeV_units = 1/(1.6021766339999998e-10)  # From Joules to GeV
SI_units = 1e4 # m2 to cm2

Comp_ph = HermesRuns + 'hermes_2D_DMph_bb_40GeV_NFW_ann_nside' + str(nside) + '.fits.gz'

hdul = pyfits.open(Comp_ph)

ipixNorth = pixbl(nside,   1, 20, -20, 20)
ipixSouth = pixbl(nside, -20, -1, -20, 20)

solidangleNorth = 4 * np.pi * np.sum(ipixNorth) / npix
solidangleSouth = 4 * np.pi * np.sum(ipixSouth) / npix

Hermes_phNorth, Hermes_phSouth, E_val = [], [], []
for iE in range(len(hdul))[1:]:
    Egamma = hdul[iE].header['ENERGY']*GeV_units
    E_val += [Egamma]
    Hermes_phNorth += [np.array(hdul[iE].data.astype(float)/SI_units)]
    Hermes_phSouth += [np.array(hdul[iE].data.astype(float)/SI_units)]

FluxHermes_North = np.array(Hermes_phNorth) @ ipixNorth / npix * 4*np.pi / 3e-26 / solidangleNorth
FluxHermes_South = np.array(Hermes_phSouth) @ ipixSouth / npix * 4*np.pi / 3e-26 / solidangleSouth
FluxHermes_Tot   = FluxHermes_North + FluxHermes_South

E_val = np.array(E_val)

hdul.close()

data = np.loadtxt('data/Fermi_40x40.txt', skiprows=2)

fig, ax = plt.subplots(figsize=(8,6))

ax.set_xscale('log')

ax.set_xlim(0.3, 50)
ax.set_ylim(-1e-6, 4.5e-6)

ax.set_xlabel(r'$E_\gamma$ [GeV]', fontsize=18)
ax.set_ylabel(r'$E_\gamma^2d\Phi_\gamma/dE_\gamma$ [GeV/cm$^2$/s/sr]', fontsize=18)

ax.plot(E_val, sigmav*E_val**2*FluxHermes_Tot, ls='-', color='k')
ax.plot(E_val, [0 for _ in E_val], ls='--', color='k')

ax.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], ls='', color='k', marker='o')

ax.text(0.15, 0.95, r'Fermi, 40°$\times$40°', color='k', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes)
ax.text(0.15, 0.87, r'gNFW, $\gamma=1.26$', color='k', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes)
ax.text(0.85, 0.95, r'DM DM $\to b\bar{b}$', color='k', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes)
ax.text(0.85, 0.87, r'$m_{DM}=40$ GeV', color='k', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes)
ax.text(0.78, 0.79, r'$\langle\sigma v\rangle = 1.3\times10^{-26}$ cm$^3$/s', color='k', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes)

fig.tight_layout(pad=0.5)
plt.savefig('GCE_flux_spec_Fermi.pdf')
plt.show()