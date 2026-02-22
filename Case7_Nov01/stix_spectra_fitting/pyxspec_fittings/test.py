import xspec
from stix2xspec.xspec_utils import *
import sunkit_spex
# 1. Use the specific path you identified
from sunkit_spex.models.physical.thermal import ThermalEmission
from sunkit_spex.models.physical.nonthermal import ThickTarget
import astropy.units as u
import numpy as np
# =========================================================
# GLOBAL INSTANCE HOLDER
# =========================================================
# We use a global variable to store the model instance so we 
# only load the atomic data ONCE, not every time XSPEC calls the function.
_THERMAL_MODEL_INSTANCE = None

def get_thermal_model():
    """
    Singleton accessor for the ThermalEmission class.
    """
    global _THERMAL_MODEL_INSTANCE
    if _THERMAL_MODEL_INSTANCE is None:
        print("Initializing ThermalEmission model (loading atomic data)...")
        # You can pass arguments here if needed (e.g., specific abundance file)
        _THERMAL_MODEL_INSTANCE = ThermalEmission()
    return _THERMAL_MODEL_INSTANCE

# =========================================================
# XSPEC WRAPPER FUNCTION
# =========================================================
def sunkit_thermal_wrapper(energies, params, flux):
    """
    Wraps sunkit_spex.models.physical.thermal.ThermalEmission for PyXspec.
    
    XSPEC Parameters:
    1. Emission Measure (1e49 cm^-3)
    2. Temperature (MK)
    """
    # 1. Get the global model instance
    model = get_thermal_model()

    # 2. Prepare Inputs
    # XSPEC provides energy bin edges in keV
    energy_edges = np.array(energies) * u.keV
    
    # Map XSPEC params to physical quantities
    # param[0]: Emission Measure. We scale by 1e49 to keep XSPEC values roughly near 1.0
    emission_measure = params[0] * 1e49 * u.cm**(-3) 
    
    # param[1]: Temperature in MegaKelvin
    temperature = params[1] * u.MK

    try:
        # 3. Calculate Flux
        # We assume the model is callable and accepts edges, temp, and EM.
        # Note: If the model expects specific keyword arguments, adjust below:
        # e.g., model(energy_edges, temperature=temperature, emission_measure=emission_measure)
        
        calculated_flux = model(energy_edges, temperature=temperature, emission_measure=emission_measure)

        # 4. Fill XSPEC flux array
        # XSPEC 'flux' array size is N-1 (bins), while 'energies' size is N (edges).
        # sunkit-spex usually returns the integrated photon flux per bin, matching N-1.
        
        # We extract .value to get raw floats (removing units)
        flux_values = calculated_flux.value
        
        for i in range(len(flux)):
            flux[i] = flux_values[i]
            
    except Exception as e:
        # Catch errors so XSPEC doesn't crash completely during fit exploration
        print(f"DEBUG: Error in thermal calc: {e}")
        for i in range(len(flux)):
            flux[i] = 0.0

# =========================================================
# REGISTER WITH XSPEC
# =========================================================

# Define parameter descriptions
# Format: "Name Unit Initial Delta Min Bot Top Max"
thermal_par_info = (
    "EM \"1e49 cm^-3\" 1.0 0.01 0.0 0.0 1e5 1e5",  # Emission Measure
    "T \"MK\" 10.0 0.1 0.1 0.1 100.0 100.0"       # Temperature
)

# # from sunkit_spex.legacy import bremsstrahlung_thick_target

# def sunkit_nonthermal_wrapper(energies, params, flux):
#     """
#     Wraps sunkit-spex thick-target bremsstrahlung.
    
#     XSPEC Params:
#     1. Total Integrated Electron Flux (10^35 electrons/s)
#     2. Spectral Index (delta)
#     3. Low Energy Cutoff (keV)
#     """
#     # 1. Prepare Energy Grid (Bin Centers & Widths)
#     # XSPEC passes bin edges. We need centers for the calculation.
#     e_edges = np.array(energies)
#     e_centers = (e_edges[:-1] + e_edges[1:]) / 2.0
#     e_widths = (e_edges[1:] - e_edges[:-1])

#     # Convert to Astropy quantities
#     energy_centers = e_centers * u.keV
    
#     # 2. Prepare Parameters
#     # Param 0: Total Electron Flux. 
#     # We scale by 1e35 to keep XSPEC parameter values readable (around ~1.0)
#     total_particle_flux = params[0] * 1e35 * (u.electron / u.s)
    
#     # Param 1: Spectral Index (Delta)
#     spectral_index = params[1]
    
#     # Param 2: Low Energy Cutoff (keV)
#     low_cutoff = params[2] * u.keV
    
#     try:
#         # 3. Call sunkit-spex model
#         # Returns Flux Density: photons s^-1 cm^-2 keV^-1
#         # Note: If your version requires 'observer_distance', add it. 
#         # Default is usually 1 AU.
#         flux_density = bremsstrahlung_thick_target(
#             energy=energy_centers, 
#             total_particle_flux=total_particle_flux, 
#             spectral_index=spectral_index, 
#             low_energy_cutoff=low_cutoff
#         )
        
#         # 4. Integrate Flux (Density * Bin Width)
#         # Result: photons s^-1 cm^-2
#         # We strip units (.value) for XSPEC
#         integrated_flux = flux_density.value * e_widths
        
#         # 5. Fill XSPEC array
#         for i in range(len(flux)):
#             flux[i] = integrated_flux[i]
            
#     except Exception as e:
#         print(f"Non-thermal Calc Error: {e}")
#         for i in range(len(flux)):
#             flux[i] = 0.0

# # =========================================================
# # REGISTER WITH XSPEC
# # =========================================================

# # Define Parameter Info
# # Format: "Name Unit Initial Delta Min Bot Top Max"
# # 1. Flux: Initial 1.0 (representing 1e35), Range 1e-5 to 1e5
# # 2. Index: Initial 4.0, Range 1.1 to 10
# # 3. Cutoff: Initial 20 keV, Range 5 to 100
# nonthermal_par_info = (
#     "F_e \"1e35 e/s\" 1.0 0.01 1e-5 1e-5 1e5 1e5", 
#     "Index \"\" 4.0 0.1 1.1 1.1 10.0 10.0",
#     "E_c \"keV\" 20.0 1.0 5.0 5.0 100.0 100.0"
# )

# # Clear (Optional: Be careful if you have other models loaded)
# # AllModels.clear()

# # Add to XSPEC
# AllModels.addPyMod(sunkit_nonthermal_wrapper, nonthermal_par_info, 'add')

# print("Non-thermal model registered! Usage: Model('sunkit_nonthermal_wrapper')")
from xspec import AllModels
AllModels.clear()

# Add to XSPEC
AllModels.addPyMod(sunkit_thermal_wrapper, thermal_par_info, 'add')

print("Model registered! usage: Model('sunkit_thermal_wrapper')")

mod_th.print_ParInfo() 