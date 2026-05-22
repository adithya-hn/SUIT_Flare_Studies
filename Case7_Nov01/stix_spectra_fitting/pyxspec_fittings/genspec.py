#####################################################
# @Author: SoLEXSPOC
# @Date:   2024-11-15 09:00:07 am
# @email: sarwade@ursc.gov.in
# @File Name: genspec.py
# @Project: solexs_tools
#
# @Last Modified time: 2026-02-05 04:31:50 pm
#####################################################

import argparse
import datetime
from astropy.io import fits
import numpy as np
import os
import warnings
from astropy.io.fits.verify import VerifyWarning
from astropy.time import Time

# from . import __version__
# from .time_utils import unix_time_to_utc
# from .caldb_utils import get_caldb_file


QUALITY_THRESHOLD_CHANNEL = 56 #2.74 - 2.79

def write_spec(channel, spec_data, stat_err, sys_err, tstart, tstop, exposure, filter_sdd, outfile, dtcorr=False, clobber=True):
    # writing file
    n_ch = len(channel)
    hdu_list = []
    primary_hdu = fits.PrimaryHDU()
                                    
    hdu_list.append(primary_hdu)

    quality = np.where(channel <= QUALITY_THRESHOLD_CHANNEL, 1, 0)

    fits_columns = []
    col1 = fits.Column(name='CHANNEL',format='1J',array=channel)
    col2 = fits.Column(name='COUNTS',format='1E',array=spec_data)
    col3 = fits.Column(name='STAT_ERR',format='1E',array=stat_err)
    col4 = fits.Column(name='SYS_ERR',format='1E',array=sys_err)
    col5 = fits.Column(name='QUALITY',format='1J',array=quality)

    fits_columns.append(col1)
    fits_columns.append(col2)
    fits_columns.append(col3)
    fits_columns.append(col4)
    fits_columns.append(col5)
    
    hdu_pha = fits.BinTableHDU.from_columns(fits.ColDefs(fits_columns))
    hdu_pha.name = 'SPECTRUM'
                                                                       
    hdu_list.append(hdu_pha)
                                                                       
    _hdu_list = fits.HDUList(hdus=hdu_list)

    # tstart_dt = datetime.datetime.fromtimestamp(tstart)
    # tstop_dt = datetime.datetime.fromtimestamp(tstop)

    tstart_utc_time_str = unix_time_to_utc(tstart)
    tstop_utc_time_str = unix_time_to_utc(tstop)
    
    # filter_sdd = hdu1[1].header['FILTER']
    # arf_file = get_caldb_file('arf',filter_sdd, obs_date) #os.path.join(CALDB_BASE_DIR,'arf',f'solexs_arf_{filter_sdd}_v{__caldb_version__}.arf')
    # rmf_file = get_caldb_file('rmf',filter_sdd, obs_date) #os.path.join(CALDB_BASE_DIR,'response','rmf',f'solexs_gaussian_{filter_sdd}_v{__caldb_version__}.rmf')

    # print(f'ARF: {arf_file}')
    # print(f'RMF: {rmf_file}')

    _hdu_list[1].header.set('TSTART',tstart_utc_time_str)
    _hdu_list[1].header.set('TSTOP',tstop_utc_time_str)
    _hdu_list[1].header.set('TIMESYS', 'UTC')
    _hdu_list[1].header.set('EXPOSURE',float(f'{exposure:.5f}'))

    
    _HEADER_KEYWORDS = (
        ("EXTNAME", "SPECTRUM", "Extension name"),
        ("CONTENT", "OGIP PHA data", "File content"),
        ("MISSION" , 'ADITYA L-1', 'Name of mission/satellite'),
        ("TELESCOP", 'AL1' , 'Name of mission/satellite'),
        ("INSTRUME", 'SoLEXS'      , 'Name of Instrument/detector'),        
        ("HDUCLASS", "OGIP    ", "format conforms to OGIP standard"),
        ("HDUVERS", "1.1.0   ", "Version of format (OGIP memo CAL/GEN/92-002a)"),
        (
            "HDUDOC",
            "OGIP memos CAL/GEN/92-002 & 92-002a",
            "Documents describing the format",
        ),
        ("HDUVERS1", "1.0.0   ", "Obsolete - included for backwards compatibility"),
        ("HDUVERS2", "1.1.0   ", "Obsolete - included for backwards compatibility"),
        ("HDUCLAS1", "SPECTRUM", "Extension contains spectral data  "),
        ("HDUCLAS2", "TOTAL ", ""),
        ("HDUCLAS3", "COUNT ", ""),
        ("HDUCLAS4", "TYPE:I ", ""),
        ("FILTER", filter_sdd, "Filter used"),
        # ('RESPFILE', rmf_file),
        # ('ANCRFILE', arf_file),
        ('BACKFILE','None'),        
        ("CHANTYPE", "PI", "Channel type"),
        ("POISSERR", False, "Are the rates Poisson distributed"),
        ("DETCHANS", n_ch, "Number of channels"),
        ("CORRSCAL", 1.0, ""),
        ("AREASCAL", 1.0, ""),
        ("BACKSCAL", 1.0, ""),
        ("SYS_ERR", 1, "Systematic error to be applied"),
#         ("QUALITY", 0, "Data quality flag"),
        ("GROUPING", 0, "Whether data is grouped"),
        ("TLMIN", 0, "Minimum legal value for 'CHANNEL' column"),
        ("TLMAX", n_ch-1, "Maximum legal value for 'CHANNEL' column"),
    )
    
    data_header = _hdu_list[1].header
    
    for k in _HEADER_KEYWORDS:
        data_header.append(k)
    
    if dtcorr:
        data_header['DTCORR'] = (True, 'Propagated from Input: Deadtime correction applied')
        data_header['HISTORY'] = 'Data derived from a deadtime-corrected type II PI file.'
    else:
        data_header['DTCORR'] = (False, 'No deadtime correction applied')

    _hdu_list[1].header = data_header
    



    _PRIMARY_HEADER_KEYWORDS = (
        ("MISSION" , 'ADITYA L-1', 'Name of mission/satellite'),
        ("TELESCOP", 'AL1' , 'Name of mission/satellite'),
        ("INSTRUME", 'SoLEXS'      , 'Name of Instrument/detector'),
        ("ORIGIN"  , 'SoLEXSPOC'       , 'Source of FITS file'),
        # ("CREATOR" , f'solexs_tools-{__version__}'  , 'Creator of file'),
        ("FILENAME", os.path.basename(outfile)            , 'Name of file'),
        ("CONTENT" , 'Type I PI file' , 'File content'),
        # ("VERSION" , __data_version__ , 'Data Product Version'),
        ("DATE", datetime.datetime.now().strftime("%Y-%m-%d"), 'Creation Date'),
    )
    
    primary_header = _hdu_list[0].header

    for k in _PRIMARY_HEADER_KEYWORDS:
        primary_header.append(k)

    _hdu_list[0].header = primary_header
        
    outfile = outfile[:-3] if outfile.endswith('.pi') else outfile
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=VerifyWarning)
        _hdu_list.writeto(f'{outfile}.pi',overwrite=clobber)

    return f'{outfile}.pi' 




def solexs_genspec(spec_file,tstart,tstop,outfile=None,clobber=True): # times in unix seconds

    hdu1 = fits.open(spec_file)
    hdu=fits.BinTableHDU.from_columns(hdu1[1].columns)
    data=hdu.data
    t0=Time(tstart)
    t1=Time(tstop)

    # time_solexs = data['TSTART']
    rate = hdu1['RATE'].data
    hdr  = hdu1['RATE'].header  
    mjdref = hdr['MJDREF']
    timezero = hdr['TIMEZERO']
    exposure = hdr['EXPOSURE']

    mjd_start = mjdref + timezero

    time_mjd = mjd_start + rate['TIME']/86400
    time_abs = Time(time_mjd, format='mjd')
    print(time_abs)
    dt=time_abs.to_value('datetime')
    # print(dt)

    

    
    # exposure=data['EXPOSURE']

    max_time = np.nanmax(dt)
    if t1 > max_time:
        warnings.warn(
            f"tstop {tstop}) is greater than the last available time in the L1 PI file ({max_time}). "
            f"Setting tstop to {max_time}.",
            UserWarning
            )
        t1 = max_time

    min_time = np.nanmin(dt)
    if t0 < min_time:
        warnings.warn(
            f"tstart {tstart}) is less than the first available time in the L1 PI file ({min_time}). "
            f"Setting tstart to {min_time}.",
            UserWarning
            )
        t0 = min_time        

    mask = (time_abs >= t0) & (time_abs < t1)
    print(f'Start Time: {t0}')
    print(f'Stop Time: {t1}')

    counts = rate * exposure[:, None]

    spec_counts = counts[mask].sum(axis=0)
    total_exposure = exposure[mask].sum()


    channel = data_f[0][3]
    n_ch = len(channel)
    spec_data = np.zeros(n_ch)
    stat_err = np.zeros(n_ch)
    sys_err = np.zeros(n_ch)
    exposure = 0

    for di in data_f:
        spec_data = spec_data + di[4]
        # stat_err = stat_err + np.sqrt(di[4])
        # sys_err = sys_err
        exposure = exposure + di[5]

    stat_err = np.sqrt(spec_data)
    
    # writing file
    tstart_dt = datetime.datetime.fromtimestamp(tstart, datetime.timezone.utc)
    tstop_dt = datetime.datetime.fromtimestamp(tstop, datetime.timezone.utc)

    if outfile == None:
        pi_file_basename = os.path.basename(spec_file)
        pi_file_basename = pi_file_basename.split('.')[0]
        outfile = pi_file_basename + '_' + tstart_dt.strftime('%H%M%S') + '_' + tstop_dt.strftime('%H%M%S')

    filter_sdd = hdu1[1].header['FILTER']

    is_dt_corr = hdu1[1].header.get('DTCORR',False)

    outfile = write_spec(channel, spec_data, stat_err, sys_err, tstart, tstop, exposure, filter_sdd, outfile, is_dt_corr, clobber)

    obs_date = hdu1[0].header.get('OBS_DATE',None)
    with fits.open(outfile,mode='update') as out_hdu:
        out_hdu[0].header['OBS_DATE'] = obs_date

        arf_file = get_caldb_file('arf',filter_sdd, obs_date) 
        rmf_file = get_caldb_file('rmf',filter_sdd, obs_date) 

        print(f'ARF: {arf_file}')
        print(f'RMF: {rmf_file}')

        out_hdu[1].header['ANCRFILE'] = arf_file
        out_hdu[1].header['RESPFILE'] = rmf_file

    return outfile


def solexs_genspec_cli():
    # Create the parser
    parser = argparse.ArgumentParser(description='Generate a type-I PI spectral file from PI spectrogram file (Type II) for a specified time range.')

    # Add arguments
    parser.add_argument('-i','--infile', type=str, help='Path to the PI spectrogram file (Type II)')
    parser.add_argument('-tstart', type=float, help='Start time in Unix seconds')
    parser.add_argument('-tstop', type=float, help='Stop time in Unix seconds')
    parser.add_argument('-gti', '--gti_file', type=str, help='Path to the Level 1 Good Time Interval File')
    parser.add_argument('-o','--outfile', type=str, help='Output file name (optional)', default=None)
    parser.add_argument('-c','--clobber', type=bool, default= False, help='Overwrite existing file if it exists')

    # Parse arguments
    args = parser.parse_args()


    try:
        outfile_name = solexs_genspec(args.infile, args.tstart, args.tstop, args.gti_file, outfile=args.outfile, clobber=args.clobber)
        print(f"Output written to {outfile_name}.")
    except Exception as e:
        print(f"Error: {e}")


def solexs_genmultispec(spec_file, tstart, tstop, time_bin, gti_file, output_dir='.', clobber=True):
    hdu1 = fits.open(spec_file)

    if hdu1[0].header['CONTENT'] != 'Type II PHA file':
        raise TypeError('Input File is not Type II PHA file.')

    hdu = fits.BinTableHDU.from_columns(hdu1[1].columns)
    data = hdu.data
    time_solexs = data['TSTART']

    hdu_gti = fits.open(gti_file)
    gti_data = hdu_gti[1].data
    gti_inds = np.array([False] * len(time_solexs))

    for i in range(len(gti_data)):
        row_gti_inds = (time_solexs >= gti_data['START'][i]) & (time_solexs <= gti_data['STOP'][i])
        gti_inds[row_gti_inds] = True

    max_time = np.nanmax(time_solexs)
    if tstop > max_time:
        warnings.warn(
            f"tstop ({tstop}) is greater than the last available time in the L1 PI file ({max_time}). "
            f"Setting tstop to {max_time}.",
            UserWarning
        )
        tstop = max_time

    pi_file_basename = os.path.basename(spec_file)
    pi_file_basename = pi_file_basename.split('.')[0]

    current_tstart = tstart
    while current_tstart < tstop:
        current_tstop = min(current_tstart + time_bin, tstop)

        # Filter data for the current time bin and GTI
        inds = (time_solexs >= current_tstart) & (time_solexs < current_tstop) & gti_inds
        data_f = data[inds]

        if len(data_f) == 0:
            warnings.warn(
                f"No valid data found for the time range ({current_tstart} to {current_tstop}). Skipping.",
                UserWarning
            )
            current_tstart += time_bin
            continue

        # Combine spectral data for the current time bin
        channel = data_f[0][3]
        n_ch = len(channel)
        spec_data = np.zeros(n_ch)
        stat_err = np.zeros(n_ch)
        sys_err = np.zeros(n_ch)
        exposure = 0

        for di in data_f:
            spec_data = spec_data + di[4]
            # stat_err = stat_err + np.sqrt(di[4])
            # sys_err = sys_err
            exposure = exposure + di[5]

        stat_err = np.sqrt(spec_data)

        filter_sdd = hdu1[1].header['FILTER']
        current_tstart_dt = datetime.datetime.fromtimestamp(current_tstart, datetime.timezone.utc)
        current_tstop_dt = datetime.datetime.fromtimestamp(current_tstop, datetime.timezone.utc)

        outfile_name = pi_file_basename + '_' + current_tstart_dt.strftime('%H%M%S') + '_' + current_tstop_dt.strftime('%H%M%S')
        outfile = os.path.join(output_dir,outfile_name)

        is_dt_corr = hdu1[1].header.get('DTCORR',False)

        tmp_outfile = write_spec(channel, spec_data, stat_err, sys_err, current_tstart, current_tstop, exposure, filter_sdd, outfile, is_dt_corr, clobber)

        obs_date = hdu1[0].header.get('OBS_DATE',None)
        with fits.open(tmp_outfile,mode='update') as out_hdu:
            out_hdu[0].header['OBS_DATE'] = obs_date

            arf_file = get_caldb_file('arf',filter_sdd, obs_date) 
            rmf_file = get_caldb_file('rmf',filter_sdd, obs_date) 

            print(f'ARF: {arf_file}')
            print(f'RMF: {rmf_file}')

            out_hdu[1].header['ANCRFILE'] = arf_file
            out_hdu[1].header['RESPFILE'] = rmf_file            

        print(f"Generated spectrum for time range {current_tstart_dt.isoformat()} to {current_tstop_dt.isoformat()}: {outfile}")

        current_tstart += time_bin


def solexs_genmultispec_cli():
    # Create the parser
    parser = argparse.ArgumentParser(description='Generate multiple type-I PI spectral files from PI spectrogram file (Type II) for a specified time range and time binning.')

    # Add arguments
    parser.add_argument('-i','--infile', type=str, help='Path to the PI spectrogram file (Type II)')
    parser.add_argument('-tstart', type=float, help='Start time in Unix seconds')
    parser.add_argument('-tstop', type=float, help='Stop time in Unix seconds')
    parser.add_argument('-tbin', '--time_bin', type=float, help='Time bin size in seconds')
    parser.add_argument('-gti', '--gti_file', type=str, help='Path to the Level 1 Good Time Interval File')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='Directory to store the generated spectra')
    parser.add_argument('-c','--clobber', type=bool, default= False, help='Overwrite existing file if it exists')
    # Parse arguments
    args = parser.parse_args()

    tstart_utc_time_str = unix_time_to_utc(args.tstart)
    tstop_utc_time_str = unix_time_to_utc(args.tstop)

    print(f'Start Time: {tstart_utc_time_str}')
    print(f'Stop Time: {tstop_utc_time_str}')
    print(f'Time Bin: {args.time_bin} seconds')
    print(f'Output Directory: {args.output_dir}')

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    try:
        solexs_genmultispec(
            spec_file=args.infile,
            tstart=args.tstart,
            tstop=args.tstop,
            time_bin=args.time_bin,
            gti_file=args.gti_file,
            output_dir=args.output_dir,
            clobber=args.clobber
        )
        print("Spectra generation completed successfully.")
    except Exception as e:
        print(f"Error: {e}")


stix_spec_file='stx_spectrum_2410315184.fits'
tstart='2024-10-31T22:58:00'
tstop ='2024-10-31T22:58:30'

outfile_name = solexs_genspec(stix_spec_file, tstart, tstop, outfile='test_time_bin.pi', clobber=True)
print(f"Output written to {outfile_name}.")