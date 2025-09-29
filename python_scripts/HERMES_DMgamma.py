from pyhermes import *
from pyhermes.units import *
import numpy as np

MinEnergy = 0.1  # GeV
MaxEnergy = 50   # GeV
E_points = 50

nside = 64

gamma_slope = 1.26
concentration = 11
M_200 = 1.14e12 * sun_mass

# rho_s = 0.243 GeV/cm3 and r_s = 20 kpc
# rho_sun = 0.4 GeV/cm3

dmProfile = darkmatter.NFWGProfile(gamma_slope, concentration, M_200)
dmSpectrum = darkmatter.PPPC4DMIDSpectrum(darkmatter.Channel.b, darkmatter.Mass.m40GeV)
# sigmav always 3e-26 cm3/s

mask_edges = ([90*deg, -90*deg], [360*deg, 0*deg])
mask = RectangularWindow(*mask_edges)

integratorDM = DarkMatterIntegrator(dmSpectrum, dmProfile)

sun_pos = Vector3QLength(8.33*kpc, 0*pc, 0*pc)
integratorDM.setObsPosition(sun_pos)

skymap_range = GammaSkymapRange(nside, MinEnergy*GeV, MaxEnergy*GeV, E_points)

skymap_range.setMask(mask)
    
skymap_range.setIntegrator(integratorDM)
skymap_range.compute()
nameDM = "{}_nside{}".format('hermes_2D_DMph_bb_40GeV_NFW_ann', nside)
skymap_range.save(outputs.HEALPixFormat("!{}.fits.gz".format(nameDM)))

print("\nDone!\n")
