import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1

fig, axs = plt.subplots(1, 2, figsize=(15,6))
rc('font', size=18)

DragonRuns = '/Users/jordankoechler/dragon-3.1.0/output/'

Voyager1 = np.loadtxt('./data/Voyager1_2016.txt', skiprows=2)
xdata    = Voyager1[:, 2]
xerrdata = [Voyager1[:, 1]-Voyager1[:, 2], Voyager1[:, 2]-Voyager1[:, 0]]
ydata    = Voyager1[:, 3]
yerrdata = Voyager1[:, 4]

run = pyfits.open(DragonRuns + 'run_2D_DMe_mm_105.7MeV_vA_0kms_NFW_ann.fits.gz')
run_header = run[0].header
run_ep = run[1].data
run_em = run[2].data

emin   = run_header['ekmin']
ek_fac = run_header['ekin_fac']
dimE   = run_header['dimE'] 
E = [emin*(ek_fac**i) for i in range(0,dimE)]
E = np.array(E)

iRsun = run_header['irsun']
izsun = run_header['izsun']

axs[0].plot(E * 1e3, (run_ep[izsun, iRsun]+run_em[izsun, iRsun]) * E**2 * 1e3, lw=2, color='red', ls="--", label=r'$v_A = 0$ km/s')
axs[0].errorbar(xdata, xdata**2 * ydata, xerr=xerrdata, yerr=np.sqrt((2*xdata*ydata*xerrdata)**2+(xdata**2*yerrdata)**2), ls='', color='blue')

run = pyfits.open(DragonRuns + 'run_2D_DMe_mm_105.7MeV_vA_13.4kms_NFW_ann.fits.gz')
run_header = run[0].header
run_ep = run[1].data
run_em = run[2].data

emin   = run_header['ekmin']
ek_fac = run_header['ekin_fac']
dimE   = run_header['dimE'] 
E = [emin*(ek_fac**i) for i in range(0,dimE)]
E = np.array(E)

iRsun = run_header['irsun']
izsun = run_header['izsun']

axs[0].plot(E * 1e3, (run_ep[izsun, iRsun]+run_em[izsun, iRsun]) * E**2 * 1e3, lw=2, color='red', ls='-', label=r'$v_A = 13.4$ km/s')

axs[0].set_xlabel(r"$E_e$ [MeV]", fontsize=18)
axs[0].set_ylabel(r"$E_e^2 \cdot \frac{d\Phi_e}{dE_e}$ at Earth [MeV s$^{-1}$ m$^{-2}$ sr$^{-1}$]", fontsize=18)

axs[0].legend(loc='upper right', fontsize=18, frameon=True)

axs[0].grid(color='lightgrey')

axs[0].set_xscale("log")
axs[0].set_yscale("log")

axs[0].set_xlim(1e0, 5e3)
axs[0].set_ylim(1e1, 1e7)

axs[0].text(0.25, 0.95, r'DM DM $\to \mu^+ \mu^- \to e^+e^-$', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[0].transAxes)

axs[0].text(0.25, 0.88, r'$\langle\sigma v\rangle = 2.3\times10^{-26}$ cm$^3$/s', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[0].transAxes)

axs[0].text(0.25, 0.81, r'$m_{DM} = 105.7$ MeV', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[0].transAxes)

axs[0].text(0.25, 0.45, r'Voyager1', color='blue', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[0].transAxes)



run = pyfits.open(DragonRuns + 'run_2D_DMe_mm_1GeV_vA_0kms_NFW_ann.fits.gz')
run_header = run[0].header
run_ep = run[1].data
run_em = run[2].data

emin   = run_header['ekmin']
ek_fac = run_header['ekin_fac']
dimE   = run_header['dimE'] 
E = [emin*(ek_fac**i) for i in range(0,dimE)]
E = np.array(E)

iRsun = run_header['irsun']
izsun = run_header['izsun']

axs[1].plot(E * 1e3, (run_ep[izsun, iRsun]+run_em[izsun, iRsun]) * E**2 * 1e3, lw=2, color='red', ls="--", label=r'$v_A = 0$ km/s')

run = pyfits.open(DragonRuns + 'run_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann.fits.gz')
run_header = run[0].header
run_ep = run[1].data
run_em = run[2].data

emin   = run_header['ekmin']
ek_fac = run_header['ekin_fac']
dimE   = run_header['dimE'] 
E = [emin*(ek_fac**i) for i in range(0,dimE)]
E = np.array(E)

iRsun = run_header['irsun']
izsun = run_header['izsun']

axs[1].plot(E * 1e3, (run_ep[izsun, iRsun]+run_em[izsun, iRsun]) * E**2 * 1e3, lw=2, color='red', ls='-', label=r'$v_A = 13.4$ km/s')
axs[1].errorbar(xdata, xdata**2 * ydata, xerr=xerrdata, yerr=np.sqrt((2*xdata*ydata*xerrdata)**2+(xdata**2*yerrdata)**2), ls='', color='blue')

axs[1].set_xlabel(r"$E_e$ [MeV]", fontsize=18)
axs[1].set_ylabel(r"$E_e^2 \cdot \frac{d\Phi_e}{dE_e}$ at Earth [MeV s$^{-1}$ m$^{-2}$ sr$^{-1}$]", fontsize=18)

axs[1].legend(loc='upper right', fontsize=18, frameon=True)

axs[1].grid(color='lightgrey')

axs[1].set_xscale("log")
axs[1].set_yscale("log")

axs[1].set_xlim(1e0, 5e3)
axs[1].set_ylim(1e1, 1e7)

axs[1].text(0.25, 0.95, r'DM DM $\to \mu^+ \mu^- \to e^+e^-$', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[1].transAxes)

axs[1].text(0.25, 0.88, r'$\langle\sigma v\rangle = 2.3\times10^{-26}$ cm$^3$/s', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[1].transAxes)

axs[1].text(0.25, 0.81, r'$m_{DM} = 1$ GeV', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[1].transAxes)

axs[1].text(0.25, 0.45, r'Voyager1', color='blue', fontsize=18, horizontalalignment='center',
        verticalalignment='center', transform=axs[1].transAxes)

fig.tight_layout(pad=0.5)
plt.savefig('DMe_Flux_vA.pdf')
plt.show()
