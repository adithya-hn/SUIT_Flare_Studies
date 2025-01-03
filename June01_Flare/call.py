import ImagesToMovie_pkg

#shiftFile='/data/sreejith/MCNS_POC/Jitter_corrections/Jun01_1.1fits_Jitter_shift.csv'
dat_dir='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data/Coloured_ROI_imgs/NB08_video/'
#ImagesToMovie_pkg.Filter_imgs(shiftFile,dat_dir)
day='01'
out_movie_nm=f'June{day}ROI_NB08_lv1.1_movie.mp4'
ImagesToMovie_pkg.Make_movie(dat_dir,out_movie_nm,24)