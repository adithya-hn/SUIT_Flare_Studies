#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mon Nov  3 05:34:31 PM CET 2025
@author: janmejoy
@hostname: machine

DESCRIPTION
- Remove contaminant spots from SUIT RoI images.
- Uses cross correlation to align first 10 images of sequence.
- Flat field is generated based on first frame of sequence.
- Correction will only work as long as RoI position does not change due to tracking.
"""

import matplotlib.pyplot as plt 
import glob 
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
from sunpy.map import Map, MapSequence
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunkit_image.coalignment import apply_shifts, calculate_match_template_shift
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
import os
from matplotlib.widgets import RectangleSelector
import numpy as np
from astropy.io import fits
from collections import defaultdict

def select_roi_with_mouse(sunpy_map, cmap=None, norm=None):
    """
    DESCRIPTION: To select RoI template.
    INPUT: Sunpy map.
    RETURNS: Sunpy submap. 
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sunpy_map)
    ax.set_title("Select ROI (click and drag) then close the window")
    sunpy_map.plot(axes=ax)
    coords = []
    def onselect(eclick, erelease):
        coords.append((eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata))

    toggle_selector = RectangleSelector(ax, onselect, useblit=True,
                      button=[1], minspanx=5, minspany=5, spancoords='pixels',
                      interactive=True)
    plt.show()

    if not coords:
        raise RuntimeError("ROI selection cancelled or failed.")

    x1, y1, x2, y2 = coords[0]
    
    bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
    top_right = (max(x1, x2), max(y1, y2)) * u.pix

    submap = sunpy_map.submap(bottom_left=bottom_left, top_right=top_right)
    return submap

def align_maps(map_seq):
    """
    DESCRIPTION:
    - Align a map sequence to the reference frame.
    - First frame of a sequence is the reference frame.
    INPUT:Unaligned map sequence.
    RETURN: Aligned sunpy map sequence.
   """
    # Ref submap is taken from the first frame of the sequence
    ref_submap = select_roi_with_mouse(map_seq[0]) 
    nt = len(map_seq)
    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)
    shifts = calculate_match_template_shift(map_seq, template=ref_submap)
    xshift_arcseconds = shifts["x"]
    yshift_arcseconds = shifts["y"]
    for i, m in enumerate(map_seq):
        xshift_keep[i] = xshift_arcseconds[i] / m.scale[0]
        yshift_keep[i] = yshift_arcseconds[i] / m.scale[1]
    map_seq = apply_shifts(map_seq, -yshift_keep, -xshift_keep, clip=False)
    final_seq_2 = []
    for i,j in enumerate(map_seq):
        date = j.date.strftime('%H:%M:%S')
        dhobt_dt = j.meta['dhobt_dt']
        grt_dt = j.meta['grt_dt']
        mfgdate = j.meta['mfgdate']
        t_obs = j.meta['t_obs']
        date_obs = j.meta['date-obs']
        obs_strt = j.meta['obs_strt']
        obs_stp = j.meta['obs_stp']
        crtime = j.meta['crtime']
        exptime = j.meta['cmd_expt']
        meas_exptime = j.meta['meas_exp']
        p = Map(j.data, map_seq[0].meta)
        p.meta['dhobt_dt'] = dhobt_dt
        p.meta['grt_dt'] = grt_dt
        p.meta['mfgdate'] = mfgdate
        p.meta['t_obs'] = t_obs
        p.meta['date-obs'] = date_obs
        p.meta['obs_strt'] = obs_strt
        p.meta['obs_stp'] = obs_stp
        p.meta['crtime'] = crtime
        p.meta['cmd_expt'] = exptime
        p.meta['meas_exp'] = meas_exptime
        final_seq_2.append(p)
    final_seq_2 = Map(final_seq_2, sequence=True)
    return (final_seq_2)

def generate_flat(al_seq, savepath,SAVE=False,MODE='median'):
    """
    DESCRIPTION:
    - Generate a flat field of contaminant spots.
    - Returns the flat field with the option to save or not.
    INPUT: Aligned map sequence, option to SAVE (bool).
    RETURNS: Flat field image array.
    """
    ref_map=al_seq[0]
    aligned_maps= [m.data for m in al_seq]

    if (MODE=='median'):
        aligned_map_arr= np.stack(aligned_maps, axis=0)
        combined_image= np.median(aligned_map_arr, axis=0)
    elif (MODE=='max'):
        aligned_map_arr= np.stack(aligned_maps, axis=0)
        combined_image= np.max(aligned_map_arr, axis=0)
    else:
        print("Specify image combination mode")
    flat= ref_map.data/combined_image
    if SAVE:
        fits.writeto(savepath, flat, overwrite= True)
        print('Flat file saved at', savepath)
    return flat

def visualize(map1, flatframe, map3,ref_map):
    """
    DESCRIPTION:
    See a preview of the individual images.
    INPUT:
    - map1, map3: Sunpy Map
    - flatframe: numpy 2D array
    RETURNS: Visualization.
    """
    VMN= 0
    VMX= np.max(ref_map.data)
    fig, ax= plt.subplots(1,3, sharex=True, sharey=True)
    im0= ax[0].imshow(map1.data, origin='lower', vmin=VMN)
    ax[0].set_title('Raw')
    plt.colorbar(im0, ax=ax[0])
    
    im1= ax[1].imshow(flatframe, origin='lower', vmin=VMN, vmax=1.2)
    ax[1].set_title('Calibration frame')
    plt.colorbar(im1, ax=ax[1])
    
    im2= ax[2].imshow(map3.data, origin='lower', vmin=VMN, vmax=VMX)
    ax[2].set_title('Corrected map')
    plt.colorbar(im2, ax=ax[2])
    plt.show()

def visualize_seq (seq1,seq2,flat_img):

    n_frames = len(seq1)

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(15, 6), sharex=True, sharey=True
    )

    im1 = ax1.imshow(seq1[0].data, origin='lower', cmap='gray')
    im2 = ax2.imshow(seq2[0].data, origin='lower', cmap='gray')
    im3 = ax3.imshow(flat_img, origin='lower', cmap='gray')

    ax1.set_title("Original")
    ax2.set_title("Corrected")
    ax3.set_title("Flat Frame")

    plt.tight_layout()
    plt.pause(1)
    plt.show(block=True)

    for i in range(n_frames):
        im1.set_data(seq1[i].data)
        im2.set_data(seq2[i].data)
        # ref_img stays fixed

        fig.suptitle(f"Frame {i}", fontsize=14)
        plt.pause(0.05)    # playback speed

def play_sequences(seq1, seq2, ref_img):
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(15, 6), sharex=True, sharey=True
    )

    im1 = ax1.imshow(seq1[0].data, origin='lower', cmap='gray')
    im2 = ax2.imshow(seq2[0].data, origin='lower', cmap='gray')
    im3 = ax3.imshow(ref_img, origin='lower', cmap='gray')

    n = len(seq1)

    for i in range(n):
        im1.set_data(seq1[i].data)
        im2.set_data(seq2[i].data)
        fig.suptitle(f"Frame {i}")
        plt.pause(0.05)

    return fig, (ax1, ax2, ax3), (im1, im2, im3)

import time

def play_sequences_loop(seq1, seq2, ref_img, delay=0.05):
    n = len(seq1)

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(15, 6), sharex=True, sharey=True
    )

    im1 = ax1.imshow(seq1[0].data, origin='lower', cmap='gray')
    im2 = ax2.imshow(seq2[0].data, origin='lower', cmap='gray')
    im3 = ax3.imshow(ref_img, origin='lower', cmap='gray')

    ax1.set_title("Seq1")
    ax2.set_title("Seq2")
    ax3.set_title("Reference")
    plt.tight_layout()

    # Show window immediately (non-blocking)
    plt.show(block=False)

    while plt.fignum_exists(fig.number):      # run until user closes window
        for i in range(n):
            if not plt.fignum_exists(fig.number):
                break

            im1.set_data(seq1[i].data)
            im2.set_data(seq2[i].data)

            fig.suptitle(f"Frame {i}")
            fig.canvas.draw_idle()

            # Control playback speed
            time.sleep(delay)

    return fig, (ax1, ax2, ax3), (im1, im2, im3)


if __name__=='__main__':
    MODE='median' # options: 'median' or 'max'
    PLOT= True # Plot preview?
    SAVE=False # Save o/p fits?

    #project_path = os.path.abspath("..")
    #savepath= os.path.join(project_path, "data/interim/roi_flat_frame.fits")
    raw_files_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/raw/'
    flt='NB02'
    
    raw_files=glob.glob(raw_files_pth+f'*3{flt}.fits')
    groups = defaultdict(list)
    maps=Map(raw_files,sequence=True)
    # Read headers, see their  and Y position chanage as function of time
    for i in range(len(maps)):
        mp= maps[i]
        x1 = mp.meta['X1']
        x2 = mp.meta['X2']
        groups[(x1, x2)].append(raw_files)

    for (x1, x2), f_seq in groups.items():
        i+=1
        if len(f_seq) < 10:
            print('Not enough files to stcak')
            continue
        print(x1,x2,'--',len(f_seq))
        print(f_seq)
        savepath=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/contam_correct/roi_flat_frame_{i}.fits'
        #f_seq = sorted(glob.glob(os.path.join(project_path,'data/raw/*.fits')))
        seq = Map(f_seq[0], sequence=True)
        seq.peek()
        aligned_sequence= align_maps(seq[:10]) # Generate flat using first 10 images of sequence
        aligned_sequence.peek()
        flat_frame= generate_flat(aligned_sequence, savepath,SAVE)
        ref_map= aligned_sequence[0]
        corrected_map_ls=[]
        for m in seq: #Multiprocess to be implemented here
            corrected_map= Map(m.data/flat_frame, m.meta)
            img_savepath='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/Processed/contam_corr_data/'+m.meta['F_NAME'] #os.path.join(project_path, 'products/roi', m.meta['F_NAME'])
            corrected_map.save(img_savepath, overwrite=True)
            corrected_map_ls.append(corrected_map)
        #if PLOT: visualize(seq, flat_frame, corrected_map_ls,ref_map) 
        if PLOT: visualize_seq(seq, corrected_map_ls,flat_frame, corrected_map_ls)   
