import ffmpeg
import pandas as pd
import numpy as np
import glob
import os
from datetime import datetime
import pathlib
import shutil
#
def move_files(mv_list_files,folder_files,f_source,f_dest):
    rm_fl=[]
    for t in mv_list_files:
        idx=np.where(folder_files==t)
        #print(idx,t)
        try:
            shutil.move((os.path.join(f_source,(folder_files[idx[0]][0])+'.jpg')),f_dest)
            rm_fl.append(folder_files[idx[0]][0])
            #print(rm_fl)

        except:
            pass
    return rm_fl

def Filter_imgs(jitFile,DataFold):
    d=[]
    data= (pd.read_csv(jitFile, sep=',')).transpose()
    Data=data.values
    mvDir=os.path.join(DataFold,'Bad_images')
    pathlib.Path(mvDir).mkdir(parents=True, exist_ok=True)
    pos1=np.where(abs(Data[2])>50)
    files = sorted(glob.glob(DataFold+'*.jpg'))
    fold_fnm=np.array([ (os.path.basename(f))[:-4] for f in files]) #jpg
    list_fnm=np.array([ f[:-5] for f in Data[0]]) #fits
    pos1=np.array(pos1)

    mv_list_files=list_fnm[pos1[0]]
    print('====================')
    print("Found few off images: ",mv_list_files)
    rm_file=move_files(mv_list_files,fold_fnm,DataFold,mvDir)
    print('Moved Bad images (x):',rm_file)
    pos4=np.where(abs(Data[3])>50)
    mv_list_files2=list_fnm[pos4[0]]
    rm_file=move_files(mv_list_files2,fold_fnm,DataFold,mvDir)
    print("Found few off images: ",mv_list_files2)
    print('Moved Bad images (y):',rm_file)
    print('====================')


def Make_movie(input_pattern,output_file,rate): 
    input_pattern=input_pattern+'*.jpg'
    input_stream = ffmpeg.input(input_pattern, pattern_type='glob', framerate=rate)

    # Define the output, including scaling to ensure even dimensions
    output_stream = ffmpeg.output(
        input_stream,
        output_file,
        vf='scale=ceil(iw/2)*2:ceil(ih/2)*2',
        vcodec='libx264',
        crf=18,
        preset='slow'
    )

    # Run the ffmpeg command
    ffmpeg.run(output_stream)
