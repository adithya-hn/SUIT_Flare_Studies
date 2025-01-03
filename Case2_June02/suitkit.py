"""
Created on 20 July 2024

This project contains a comprehensive module for suit image analysis, primarily focused on the processing and analysis of solar images in FITS format. 
The scripts and functions provided are designed to filter, process, and align these images and subsequently generate movies from the processed data. 

Author: Adithya H.N. 
Purpose: A comprehensive module for SUIT image analysis

Major function:
    - Image coalign
    - Save images
    - View maps
    - List the files in the time range

Modification History:

"""

from sunpy.time import parse_time
import os
from datetime import datetime
#import datetime
import csv
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
from sunkit_image.coalignment import calculate_match_template_shift, apply_shifts
from datetime import timedelta
import timeit
import pathlib
from colormap import filterColor
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import numpy as np
from PIL import Image
import ImagesToMovie_pkg

def suit_co_align_fd_imgs(search_fold,filter_name,add_logos,batch_size=10,rate=30,ref_idx=0,Save_fits=False,test_mode=False,start_idx=0,end_idx=21):
    """
    Function to co-align full-disk images, add logos, create movies, and save results.

    Parameters:
    - search_fold: Directory to search for images.
    - filter_name: Name of the filter.
    - add_logos: Whether to add logos to images.
    - batch_size: Number of images to process in each batch.
    - rate: Frame rate for the output movie.
    - ref_idx: Reference index for alignment.
    - Save_fits: Whether to save the aligned FITS files.
    - test_mode: Whether to run in test mode with a limited number of files.
    - start_idx: Start index for test mode.
    - end_idx: End index for test mode.
    """
    global base_fold
    start = timeit.default_timer()
    now = datetime.now() - timedelta(days=1)
    print("Starting up now: ", datetime.now())
    print(f'Searching for {filter_name} images in {search_fold} folder')
    #print('Path: ',os.getcwd())
    base_fold=os.getcwd()
    #Fetch and sort files
    files = get_sorted_files(search_fold, filter_name)
    if test_mode:
        files=files[start_idx:end_idx]
    print('Total files:', len(files))

    ref_img = sunpy.map.Map(files[0])
    ref_submap = get_submap(ref_img)

    # Create directories
    fl_date =datetime.fromisoformat(str(ref_img.date))
    fol_nm, jpg_fold, algn_dir = create_directories(fl_date,Save_fits)

    o_x, o_y, o_d, od, x_arry, y_arry, aln_imgs = [], [], [], [], [], [], []

    for i in range(0, len(files), batch_size):
        batch_files = [files[0]] + files[i:i + batch_size]
        batch_results = process_batch(batch_files, ref_submap, add_logos, jpg_fold,Save_fits,algn_dir)
        o_x.extend(batch_results[0])
        o_y.extend(batch_results[1])
        o_d.extend(batch_results[2])
        od.extend(batch_results[3])
        x_arry.extend(batch_results[4])
        y_arry.extend(batch_results[5])
        aln_imgs.extend(batch_results[6])

    # Save results and create movie
    save_results(fl_date, o_x, o_y, o_d, x_arry, y_arry, aln_imgs)
    create_movie(jpg_fold, fl_date,rate)

    stop = timeit.default_timer()
    print('Done, Time: ',datetime.now())
    print('Run Time: ', (stop - start) / 60, 'Mins')

def get_sorted_files(search_fold, filter_name):
    """
    Fetch and sort FITS files based on timestamp.

    Parameters:
    - search_fold: Directory to search for images.
    - filter_name: Name of the filter.

    Returns:
    - Sorted list of file paths.
    """
    fdir = search_fold 
    files = glob.glob(fdir + f'*2{filter_name}.fits')
    files = sorted(files, key=lambda file_name: datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    return files

def get_sorted_roi_files(search_fold, filter_name):
    """
    Fetch and sort FITS files based on timestamp.

    Parameters:
    - search_fold: Directory to search for images.
    - filter_name: Name of the filter.

    Returns:
    - Sorted list of file paths.
    """
    fdir = search_fold 
    files = glob.glob(fdir + f'*3{filter_name}.fits')
    files = sorted(files, key=lambda file_name: datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    return files

def get_sorted_all_files(search_fold, filter_name):
    """
    Fetch and sort FITS files based on timestamp.

    Parameters:
    - search_fold: Directory to search for images.
    - filter_name: Name of the filter.

    Returns:
    - Sorted list of file paths.
    """
    fdir = search_fold 
    files = glob.glob(fdir + f'*{filter_name}.fits')
    files = sorted(files, key=lambda file_name: datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    return files

def get_submap(ref_img):
    """
    Get a submap from the reference image for template matching.

    Parameters:
    - ref_img: Reference SunPy map image.

    Returns:
    - Submap of the reference image.
    """
    center_coord = SkyCoord(0 * u.arcsec, 950* u.arcsec, frame=ref_img.coordinate_frame) #54,157
    width = 1100 * u.arcsec
    height =300 * u.arcsec   
    
    #ref_img.meta.update({'CROTA2':0})
    #print('rot',ref_img.meta.get('CROTA2')) 
    offset_frame = SkyOffsetFrame(origin=center_coord, rotation=0*u.deg)
    rectangle = SkyCoord(lon=[-1/2, 1/2] * width, lat=[-1/2, 1/2] * height, frame=offset_frame)
    ref_submap = ref_img.submap(rectangle) #bottom_left, top_right=top_right)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=ref_submap)
    ref_submap.plot(axes=ax)
    plt.savefig('Template.jpg')
    #plt.show()
    plt.close()
    return ref_submap

def create_directories(fl_date,Save_fits):
    """
    Create necessary directories for saving results.

    Parameters:
    - fl_date: Date of the reference image.
    - Save_fits: Whether to save the aligned FITS files.

    Returns:
    - Directory paths for saving results.
    """
    fol_nm = f'{base_fold}/{str(fl_date.day).zfill(2)}_{str(fl_date.month).zfill(2)}_{str(fl_date.year).zfill(2)}'
    jpg_fold = f'{fol_nm}/FD_imgs'
    video_fold=f'{base_fold}/Website_Movies'
    algn_dir = f'{fol_nm}/Aligned_Fits'
    log_fold=f'{base_fold}/Drift_log'
    
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(video_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(log_fold).mkdir(parents=True, exist_ok=True)
    if Save_fits:
        pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)

    
    return fol_nm, jpg_fold, algn_dir

def process_batch(files, ref_submap, add_logos, jpg_fold,Save_fits,algn_dir):
    """
    Process a batch of files for alignment and enhancement.

    Parameters:
    - files: List of file paths.
    - ref_submap: Submap for template matching.
    - add_logos: Whether to add logos to images.
    - jpg_fold: Directory to save JPEG images.
    - Save_fits: Whether to save the aligned FITS files.
    - algn_dir: Directory to save aligned FITS files.

    Returns:
    - Results of the batch processing.
    """    
    o_x, o_y, o_d, od, x_arry, y_arry, aln_imgs = [], [], [], [], [], [], []
    ref_head=ref_submap.fits_header
    ref_cdel=ref_head['CDELT1']

    Sequence = sunpy.map.Map(files, sequence=True)
    for l in range(len(Sequence) - 1):
        o_x.append(Sequence[l + 1].meta.get('CRPIX1'))
        o_y.append(Sequence[l + 1].meta.get('CRPIX2'))
        o_d.append(datetime.fromisoformat(str(Sequence[l + 1].date)))
        od.append(str(Sequence[l + 1].date))

    align_shift = calculate_match_template_shift(Sequence, template=ref_submap)
    x_arry.extend(align_shift['x'].value[1:] / ref_cdel * -1)  # avoiding the inserted image data
    y_arry.extend(align_shift['y'].value[1:] / ref_cdel * -1)
    shift_xPix = align_shift['x'].value / ref_cdel * -1
    shift_yPix = align_shift['y'].value / ref_cdel * -1

    aligned_maps = apply_shifts(Sequence, yshift=shift_yPix * u.pixel, xshift=shift_xPix * u.pixel, clip=False)
    for k in range(len(aligned_maps) - 1):  # account for inserted image, j=0 is reference image
        j = k + 1
        img_head = aligned_maps[j].fits_header
        img_head['CRPIX1'] = ref_head['CRPIX1'] + shift_xPix[j]
        img_head['CRPIX2'] = ref_head['CRPIX2'] + shift_yPix[j]

        # Image enhancement
        limb_enh_data = enhance_image(aligned_maps[j], img_head)
        aligned_img = sunpy.map.Map(aligned_maps[j].data, img_head)
        if Save_fits:
            fits_fnm=f'{algn_dir}/{os.path.basename(files[j])}'
            aligned_img.save(fits_fnm,overwrite=True)

        aln_imgs.append(os.path.basename(files[j]))  # element is added to sequence not original file list so k is correct than j
        fl_nm = f'{jpg_fold}/{os.path.basename(files[j])[:-4]}jpg'
        save_image(aligned_img, add_logos, fl_nm)

    return o_x, o_y, o_d, od, x_arry, y_arry, aln_imgs

def enhance_image(aligned_map, img_head):
    """
    Enhance the image by applying amplification to off-limb features.

    Parameters:
    - aligned_map: Aligned SunPy map image.
    - img_head: Header of the image.

    Returns:
    - Enhanced image data.
    """
    off_limb_len = 100  # pixels beyond off limb
    amp = 10  # off disk amplification factor
    h, w = aligned_map.data.shape
    col, row, radius = img_head['CRPIX1'] - 3, img_head['CRPIX2'] - 4, img_head['R_SUN']
    mask = np.ones((h, w)) * create_circular_mask(h, w, col, row, radius)
    mask = np.where(mask < 1, amp, mask)  # contains amplification factor (amp) for off disk feature
    offdisk = aligned_map.data * mask
    outermask = np.ones((h, w)) * create_circular_mask(h, w, col, row, radius + off_limb_len)  # controls the size of outer mask
    limb_enh_data = offdisk * outermask
    return limb_enh_data

def create_circular_mask(h, w, col, row, radius):
    '''
    Creates circular mask of desired size.
    - h, w: height and width of canvas
    - col, row: column and row of circle center
    - radius: radius of circle
    '''
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - col) ** 2 + (Y - row) ** 2)
    mask = dist_from_center <= radius
    return mask

def save_image(aligned_img, add_logos, fl_nm):
    """
    Save the aligned image to a file.

    Parameters:
    - aligned_img: Aligned SunPy map image.
    - add_logos: Whether to add logos to images.
    - fl_nm: File name to save the image.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=aligned_img)
    aligned_img.plot(title=False, cmap=filterColor['NB03'], vmin=np.percentile(aligned_img.data,1), vmax=np.percentile(aligned_img.data, 99.75))
    plot_str = str(aligned_img.date) + '          ' + 'NB03 Mg II k 279 nm'
    
    if add_logos == 1:
        logo_paths = {
            #Picking Up Relevant Logos.
            "logo1": 'assets/suit_white.png',
            "logo2": 'assets/sun_iucaa.png',
            "logo3": 'assets/iucaaisro.png'
        }
        logo_image3 = Image.open(logo_paths['logo3'])  # suit_logo, top left
        logo_image2 = Image.open(logo_paths['logo2'])  # sun at iucaa top right
        logo_image1 = Image.open(logo_paths['logo1'])  # suitlogo bottom right

        logo_image3 = logo_image3.resize((100, 76))  # Adjust the size as needed
        logo_image2 = logo_image2.resize((53, 53))
        logo_image1 = logo_image1.resize((73, 50))

        logo_image1 = np.array(logo_image1)
        logo_image2 = np.array(logo_image2)
        logo_image3 = np.array(logo_image3)

        fig.figimage(logo_image1, xo=880, yo=20, zorder=3, alpha=1)
        fig.figimage(logo_image2, xo=890, yo=890, zorder=3, alpha=1)
        fig.figimage(logo_image3, xo=12, yo=880, zorder=3, alpha=1)

    ax.text(50, 50, plot_str, color='white', weight='bold', fontsize=18)

    plt.draw()
    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    plt.tight_layout()
    plt.savefig(fl_nm, bbox_inches='tight', transparent="True", pad_inches=0)
    plt.close()

def save_results(fl_date, o_x, o_y, o_d, x_arry, y_arry, aln_imgs):
    """
    Save results to log files.

    Parameters:
    - fl_date: Date of the reference image.
    - o_x: List of x-coordinates of original CRPIX1 values.
    - o_y: List of y-coordinates of original CRPIX2 values.
    - o_d: List of datetime objects of the images.
    - x_arry: List of x-coordinate shifts.
    - y_arry: List of y-coordinate shifts.
    - aln_imgs: List of aligned image file names.
    """
    data = np.column_stack((aln_imgs, o_d, x_arry, y_arry, o_x, o_y))
    months = ['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    output_file = f'{base_fold}/Drift_log/{months[int(fl_date.month)]}{str(fl_date.day).zfill(2)}_l1fits_Drift_shift.csv'
    np.savetxt(output_file, data, delimiter=',', header='File names, Date, delX, delY, CR-X, CR-Y', fmt='%s')

def create_movie(jpg_fold, fl_date,rate=30):
    """
    Create a movie from aligned images.

    Parameters:
    - jpg_fold: Directory containing JPEG images.
    - fl_date: Date of the reference image.
    - rate: Frame rate for the movie.
    """
    months = ['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    dat_dir=jpg_fold+'/'
    shiftFile=f'{base_fold}/Drift_log/'+f'{months[int(fl_date.month)]}{str(fl_date.day).zfill(2)}_l1fits_Drift_shift.csv'
    ImagesToMovie_pkg.Filter_imgs(shiftFile,dat_dir)
    movie_name = f'{base_fold}/Website_Movies/{months[int(fl_date.month)]}{str(fl_date.day).zfill(2)}_2K_NB03_l1_movie.mp4'
    ImagesToMovie_pkg.Make_movie(jpg_fold, movie_name, rate)

def find_fits_files(base_dir, start_time, end_time):
    """
    List FITS files in the directory based on the given time range without a catalogue file.
    Relatively time-consuming compared to "find_files_in_time_range"

    Parameters:
    - base_dir: Base directory to search for FITS files.
    - start_time: Start time for the filtering range.
    - end_time: End time for the filtering range.

    Returns:
    - List of selected file paths within the specified time range.
    """
    selected_files = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.fits'):
                file_path = os.path.join(root, file)
                try:
                    smap = sunpy.map.Map(file_path)
                    obs_time = parse_time(smap.date).datetime
                    if start_time <= obs_time <= end_time:
                        selected_files.append(file_path)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
    return selected_files

def list_files_in_time_range(csv_file, start_time, end_time,img_type='all',Filter=''):
    """
    List the files from the CSV catalogue based on the given time range.

    Parameters:
    - csv_file: Path to the CSV file.
    - start_time: Start time for the filtering range.
    - end_time: End time for the filtering range.
    - img_type: Image type ('ROI','2K', '4K', or 'all').
    - Filter: Specific filter to apply (default is '').

    Returns:
    - List of selected image file names within the specified time range.
    """
    if img_type== '2K':
        img_dim=2048
    if img_type=='4K':
        img_dim=4096
    if Filter=='':
        if img_type== '2K':
            Filter== 'NB03'
        else:
            Filter== 'all'
   
    print(f'Listing {img_type} images of {Filter} Filter, between time range {start_time} <--> {end_time}')


    selected_files = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            obs_time = parse_time(row['time']).datetime
            if start_time <= obs_time <= end_time:

                if img_type != 'all':
                    if img_type == 'ROI':
                        if row['NAXIS1']!='2048' and  row['NAXIS1']!='4096':
                            if Filter == '':
                                #print('Listing all ROIs')
                                selected_files.append(row['file_path'])
                            else:
                                if row['filter_name']==Filter:
                                    #print(f'Listing {Filter} filter ROIs')
                                    selected_files.append(row['file_path'])

                    else:
                        if row['NAXIS1']==str(img_dim) :
                            #print(f'Listing {img_type} images')
                            selected_files.append(row['file_path'])

                
                else:
                    #print('Listing all type of data')
                    selected_files.append(row['file_path'])
    return selected_files

def process_fits_files(file_list):
    """Process the selected FITS files."""
    for file in file_list:
        try:
            smap = sunpy.map.Map(file)
            # Example processing: Print the file name and observation date
            print(f'Processing {file}')
            print(f"Observation Date: {smap.date}")
        except Exception as e:
            print(f"Error processing {file}: {e}")

def save_jpg(sunpy_map, output_file,col_map='gray'):
    """
    Convert a SunPy map to a JPEG image.

    Parameters:
    sunpy_map (sunpy.map.Map): The SunPy map object to convert.
    output_file (str): The path to the output JPEG file.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection=sunpy_map)
    sunpy_map.plot(axes=ax, cmap=col_map,clip_interval=(1, 99.99) * u.percent)
    plt.colorbar().remove()
    ax.set_axis_off()
    plt.savefig(output_file, format='jpeg', bbox_inches='tight', pad_inches=0,dpi=300)
    plt.close(fig)

def view_map(sunpy_map,close_map=False,col_map='gray'):
    fig = plt.figure()
    ax = fig.add_subplot(projection=sunpy_map)
    sunpy_map.plot(axes=ax,cmap=col_map, clip_interval=(1, 99.99) * u.percent)
    output_file=sunpy_map.meta.get('F_NAME')[:-4]+'jpg'
    plt.savefig(output_file, format='jpeg',dpi=300)
    if close_map:
        plt.close(fig)
    if close_map is False:
        plt.show()

def map_movie(files,Filter='',out='SUIT_movie.mp4',start_time='',end_time=''):
    print('Total images: ',len(files))
    Files=[]
    if Filter != '':
        for fls in files:
            if fls.endswith(f'{Filter}.fits'):
                Files.append(fls)
                print('files:',fls)
                
    Sequence = sunpy.map.Map(Files, sequence=True)  
    fig = plt.figure()
    ax = fig.add_subplot(projection=Sequence.maps[0])
    ani = Sequence.plot(axes=ax,cmap=filterColor['NB04'])
    #plt.axis('off')
    if start_time != '':
        if Filter !='':
            ani.save(str(start_time)+'_'+str(end_time)+'_'+Filter+'_'+out)
        else:
            ani.save(str(start_time)+'_'+str(end_time)+'_'+out)
    if start_time =='':
        if Filter != '':
            ani.save(Filter+'_'+out)
        else:
            ani.save(out)

    plt.close()
    print('Done with making movies')
    

if __name__=="__main__":
    start = timeit.default_timer()
    print("Starting up now: ", datetime.now())
    '''
    search_fold = '/scratch/suit_data/level1fits/2024/07/21/normal_2k/'  # Custom Folder
    filter_name = 'NB03'
    add_logos=1 
    suit_co_align_fd_imgs(search_fold,filter_name,add_logos,test_mode=True)
    '''
    files=get_sorted_roi_files('/Analysis/Projects_Data/Flare_Data/June02_Flare_Data/Processed_Data/Aligned_images/NB03/','NB03')
    map_movie(files,Filter='NB03')

