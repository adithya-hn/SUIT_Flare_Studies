from astropy.io import fits
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt 


def stix_genspec(stix_file, tstart, tstop, outfile):

    hdul = fits.open(stix_file)
    rate_hdu = hdul['RATE']
    data = rate_hdu.data
    hdr = rate_hdu.header

    # ----- reconstruct absolute time -----

    mjdref = hdr['MJDREF']
    timezero = hdr.get('TIMEZERO',0)

    time_sec = data['TIME']

    time_mjd = mjdref + timezero + time_sec/86400.0
    time = Time(time_mjd, format='mjd')

    tstart = Time(tstart)
    tstop  = Time(tstop)

    mask = (time >= tstart) & (time < tstop)

    if np.sum(mask) == 0:
        raise ValueError("No spectral bins in requested time range")

    # ----- extract arrays -----

    rate = data['RATE'][mask]          # counts/s
    exposure = data['EXPOSURE'][mask]  # seconds
    livetime = data['LIVETIME'][mask]  #not used since i am using exposure which is already corrrected
    sys_err = data['SYS_ERR'][mask]
    stat_err_rate = data['STAT_ERR'][mask]
    timedel = data['TIMEDEL'][mask]
    # # exposure without deadtime
    total_exposure = timedel.sum()
    # print(total_exposure)

    # ----- reconstruct counts ----- # exposure is already corrected for deadtime

    counts = rate * exposure[:,None]

    # ----- integrate spectrum -----

    spec_counts = np.sum(counts, axis=0)

    stat_err_ = np.sqrt(spec_counts) #cheking

    sys_err = np.mean(sys_err, axis=0)

    # propagate statistical error
    stat_err_counts = stat_err_rate * exposure[:,None]

    # quadrature propagation
    stat_err = np.sqrt(np.sum(stat_err_counts**2, axis=0))

    plt.plot(stat_err)
    plt.plot(stat_err_)
    plt.show()

    # ----- exposure -----

    exposure = np.sum(exposure)
    # print('exposure: ',exposure,livetime.sum()/len(livetime))

    # ----- channels -----

    channel = data['CHANNEL'][0]

    nchan = len(channel)

    # ----- build PHA -----

    col1 = fits.Column(name='CHANNEL', format='1J', array=channel)
    col2 = fits.Column(name='COUNTS', format='1E', array=spec_counts)
    col3 = fits.Column(name='STAT_ERR', format='1E', array=stat_err)
    col4 = fits.Column(name='SYS_ERR', format='1E', array=sys_err)

    cols = fits.ColDefs([col1,col2,col3,col4])

    pha = fits.BinTableHDU.from_columns(cols)
    pha.name = 'SPECTRUM'

    # ----- OGIP keywords -----

    pha.header['EXTNAME']  = 'SPECTRUM'
    pha.header['HDUCLASS'] = 'OGIP'
    pha.header['HDUCLAS1'] = 'SPECTRUM'
    pha.header['HDUCLAS2'] = 'TOTAL'
    pha.header['HDUCLAS3'] = 'COUNT'
    pha.header['HDUCLAS4'] = 'TYPE:I'
    pha.header['CHANTYPE'] = 'PI'
    pha.header['DETCHANS'] = nchan
    pha.header['EXPOSURE'] = float(exposure)
    pha.header['TSTART'] = tstart.isot
    pha.header['TSTOP']  = tstop.isot
    pha.header['DTCORR'] = True
    pha.header['POISSERR'] = False

    # ----- write file -----

    fits.HDUList([fits.PrimaryHDU(),pha]).writeto(outfile,overwrite=True)

    print("Type-I spectrum written:", outfile)

stix_genspec(
    "stx_spectrum_2410315184.fits",
    "2024-11-01T02:15:23",
    "2024-11-01T02:15:26",
    "stix_012230_012300.pha"
)