#!/usr/bin/python3
"""
 Script to plot STIX preview images  with matplotlib

 Created at 2025-07-30T10:06:33.049142 by STIX data center online image reconstruction tool
 Contact:  hualin.xiao@fhnw.ch
 If this tool is helpful for your work, please cite our paper: 
 * Hualin Xiao, et al., The data center for the Spectrometer and Telescope for Imaging X-rays (STIX) on board Solar Orbiter, 2023, A&A, 673, A142, DOI: https://doi.org/10.1051/0004-6361/202346031

"""
import wget
import sunpy
import sunpy.map
from astropy import units as u
from matplotlib import pyplot as plt
#%matplotlib notebook # for jupyter notebook

# define more units
u.add_enabled_units(
    [u.def_unit("arcsecs", 1 * u.arcsec),
     u.def_unit("meters", 1 * u.m)])

download_location = '.'

#download image fits files from STIX data center
fits_urls = ['https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_bp_hpc.fits', 'https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_vff_hpc.fits', 'https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_mem_hpc.fits', 'https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_clean_hpc.fits', 'https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_full_hpc.fits', 'https://datacenter.stix.i4ds.net/image-archive/2024/11/13/16/stix_L3A_quicklook_20241113T165600-20241113T165700-7.0-50.0keV_2411130279_989_user_spec_fitting.fits']
filenames = [wget.download(url, out=download_location) for url in fits_urls]

maps = sunpy.map.Map(filenames)

print(f"total number of maps: {len(maps)}")

#plotting images using sunpy.map
for m in maps:
    plt.figure()
    m.plot(cmap="std_gamma_2")
    m.draw_grid(color='w', ls='--', grid_spacing=10 * u.deg)
    m.draw_limb(color='w')
plt.show()

