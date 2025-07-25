from pyzdcf import pyzdcf

input = './input/'           # Path to the input data
output = './output/'         # Path to the directory for saving the results

# Light curve names
lc1 = 'lc_name1'
lc2 = 'lc_name2'

# Parameters are passed to the pyZDCF as a dictionary

params = dict(autocf            =  False, # Autocorrelation (T) or cross-correlation (F)
              prefix            = 'ccf',  # Output files prefix
              uniform_sampling  =  False, # Uniform sampling?
              omit_zero_lags    =  True,  # Omit zero lag points?
              minpts            =  0,     # Min. num. of points per bin (0 is a flag for default value of 11)
              num_MC            =  100,   # Num. of Monte Carlo simulations for error estimation
              lc1_name          =  lc1,   # Name of the first light curve file
              lc2_name          =  lc2    # Name of the second light curve file (required only if we do CCF)
             )

# Here we use non-interactive mode (intr=False)
dcf_df = pyzdcf(input_dir  = input,
                output_dir = output,
                intr       = False,
                parameters = params,
                sep        = ',',
                sparse     = 'auto',
                verbose    = True)