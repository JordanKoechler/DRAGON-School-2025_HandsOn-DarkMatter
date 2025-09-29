from pyhermes import *
from pyhermes.units import *
import numpy as np

Output_Dragon = '/Users/jordankoechler/dragon-3.1.0/output/'
nside = 64
sun_pos = Vector3QLength(8.33*kpc, 0*pc, 0*pc)

MinEnergy = 2.7e-5  # GeV
MaxEnergy = 1.8e-3   # GeV
E_points = 15

kn_crosssection = interactions.KleinNishina()
CMB_photons = photonfields.CMB()
ISRF_photons = photonfields.ISRF()

mask_edges = ([90*deg, -90*deg], [360*deg, 0*deg])
mask = RectangularWindow(*mask_edges)

def Calc_IC(FName, OutName, GalMod):
    print('Variable model: ' + FName.replace(Output_Dragon, ''))

    if GalMod == '2D':  dragon_leptons = cosmicrays.Dragon2D(FName, [Electron, Positron])
    elif GalMod == '3D':  dragon_leptons = cosmicrays.Dragon3D(FName, [Electron, Positron])
    else: dragon_leptons = cosmicrays.Dragon2D([Electron, Positron])

    integratorIC_cmb  = InverseComptonIntegrator(dragon_leptons, CMB_photons, kn_crosssection)
    integratorIC_isrf = InverseComptonIntegrator(dragon_leptons, ISRF_photons, kn_crosssection)

    integratorIC_cmb.setupCacheTable(30, 30, 12)
    integratorIC_cmb.setObsPosition(sun_pos)
    integratorIC_isrf.setupCacheTable(30, 30, 12)
    integratorIC_isrf.setObsPosition(sun_pos)

    skymapIC_cmb_range = GammaSkymapRange(nside, MinEnergy*GeV, MaxEnergy*GeV, E_points)
    skymapIC_isrf_range = GammaSkymapRange(nside, MinEnergy*GeV, MaxEnergy*GeV, E_points)

    skymapIC_cmb_range.setMask(mask)
    skymapIC_isrf_range.setMask(mask)
    
    skymapIC_cmb_range.setIntegrator(integratorIC_cmb)
    skymapIC_cmb_range.compute()
    nameIC1 = "{}-IC_cmb_nside{}".format(OutName, nside)
    skymapIC_cmb_range.save(outputs.HEALPixFormat("!{}.fits.gz".format(nameIC1)))
    print("IC on CMB done")
    
    skymapIC_isrf_range.setIntegrator(integratorIC_isrf)
    skymapIC_isrf_range.compute()
    nameIC2 = "{}-IC_isrf_nside{}".format(OutName, nside)
    skymapIC_isrf_range.save(outputs.HEALPixFormat("!{}.fits.gz".format(nameIC2)))
    print("IC on ISRF done")


FN = 'run_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann.fits.gz'
OutName = 'hermes_2D_DMe_mm_1GeV_vA_13.4kms_NFW_ann'
Calc_IC(Output_Dragon + FN, OutName, '2D')

print("\nDone!\n")
