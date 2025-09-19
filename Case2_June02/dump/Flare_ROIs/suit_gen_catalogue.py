import os
import csv
import sunpy.map
import datetime
import timeit
from sunpy.time import parse_time

def create_catalogue(base_dir, output_csv):
    """Create a catalogue from FITS files in the base directory."""
    catalogue = []

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.fits'):
                file_path = os.path.join(root, file)
                
                try:
                    smap = sunpy.map.Map(file_path)
                    header = smap.meta
                    
                    image_name = file
                    time = header.get('DATE-OBS', '')
                    filter_name = header.get('FTR_NAME', '')
                    naxis1 = header.get('NAXIS1', '')
                    naxis2 = header.get('NAXIS2', '')
                    roi_id = header.get('ROI_ID', '')
                    crpix1 = header.get('CRPIX1', '')
                    crpix2 = header.get('CRPIX1', '')
                    crval1 = header.get('CRVAL1', '')
                    crval2 = header.get('CRVAL2', '')
                    rotVal = header.get('CROTA2', '')
                    solex_trigger1 = header.get('SOLX1TR', '')
                    solex_trigger2 = header.get('SOLX2TR', '')
                    helios_trigger = header.get('HELIOSTR', '')
                    
                    catalogue.append([image_name, time, filter_name, naxis1, naxis2, roi_id,
                                      crpix1, crpix2, crval1, crval2,rotVal ,solex_trigger1, 
                                      solex_trigger2, helios_trigger, file_path])
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
    
    catalogue.sort(key=lambda x: parse_time(x[1]).datetime if x[1] else None)
    # Write the catalogue to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['image_name', 'time', 'filter_name', 'NAXIS1', 'NAXIS2', 'roi_id',
                         'crpix1', 'crpix2', 'crval1', 'crval2','CrotVal', 'solex_trigger1', 
                         'solex_trigger2', 'helios_trigger', 'file_path'])
        writer.writerows(catalogue)

if __name__ == "__main__":
    base_directory = '/scratch/suit_data/level1.1fits/2024/07/'
    output_csv = 'SUIT_Data_Jul_catalogue.csv'
    start = timeit.default_timer()
    print("Starting up now: ", datetime.datetime.now())
    create_catalogue(base_directory, output_csv)
    print(f">> Catalogue created at {output_csv}")
    stop = timeit.default_timer()
    print('Done, Time: ', datetime.datetime.now())
    print('Run Time: ', (stop - start) / 60, 'Mins')
