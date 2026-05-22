import xspec
from stix2xspec.xspec_utils import *
#from sunxspex import xspec_models

#mod_th = sunxspex.xspec_models.ThermalModel()
#xspec.AllModels.addPyMod(mod_th.model, mod_th.ParInfo, 'add')
#mod_th.print_ParInfo() # see the initial configuration of parameters

#mod_nt = sunxspex.xspec_models.ThickTargetModel()
#xspec.AllModels.addPyMod(mod_nt.model, mod_nt.ParInfo, 'add')
#mod_nt.print_ParInfo() # see the initial configuration of parameters

xspec.AllData.clear() # get rid of any data that is still loaded from previous runs
xspec.AllData(f"stix_215_3s.pha") # fit the 1140th data row in the converted spectrogram file. make sure the .srm file is in the same folder as the spectrogram file.

#spectime = fits_time_to_datetime('stx_spectrum_2410315184.fits,', idx=1140)
#plot_data(xspec, erange = [4,50],title = f'STIX spectrum at {spectime:%Y-%m-%d %H:%M:%S}').show()
