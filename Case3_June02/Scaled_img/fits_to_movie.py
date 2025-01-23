import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from datetime import timedelta
import matplotlib.colors as colors
import timeit
import pathlib
import cv2
from colormap import filterColor

start = timeit.default_timer()


def fits_to_image(search_fold,pattern='',movie_name='Fits_Movie.mp4'):
    print(search_fold+f'*{pattern}.fits')
    files = sorted(glob.glob(search_fold+f'*{pattern}.fits'))
    print('Total files:',len(files))
    Sequence = sunpy.map.Map(files, sequence=True)
    fig = plt.figure()
    ax = fig.add_subplot(projection=Sequence.maps[0])
    #ani=Sequence.plot(axes=ax)
    ani = Sequence.plot(axes=ax,cmap=filterColor['NB08'],norm=colors.Normalize(vmin=0,vmax=10500,clip=True))
    #plt.axis('off')
    plt.colorbar()
    ani.save(movie_name)
    plt.close()

def Make_movie(image_folder, video_name, frame_rate):
    # Get list of images
    images = [img for img in os.listdir(image_folder) if img.endswith(".jpg") or img.endswith(".png")]
    images.sort()  # Ensure the images are in order

    # Get dimensions of the images
    first_image_path = os.path.join(image_folder, images[0])
    frame = cv2.imread(first_image_path)
    height, width, layers = frame.shape

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    video = cv2.VideoWriter(video_name, fourcc, frame_rate, (width, height))

    for image in images:
        image_path = os.path.join(image_folder, image)
        frame = cv2.imread(image_path)
        video.write(frame)

    # Release the video writer object
    video.release()
    print(f"Video {video_name} created successfully.")

if __name__ == "__main__":
    dest_fld='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/Scaled_img'
    jpg_fold=dest_fld+'/'+'Imgs'
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    search_fold='/Analysis/Projects_Data/Flare_Data/June02_Flare_Data2/Processed/Aligned_images/NB08/'
    #mv_nm=
    pattern='NB08'
    fits_to_image(search_fold,pattern='NB08')