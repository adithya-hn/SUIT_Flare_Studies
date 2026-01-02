import cv2
import glob
import re


png_files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/localize_brightenings/har_wavlets/m_scal_conts/*.png"))
#print(png_files)
img = cv2.imread(png_files[0])
h, w, _ = img.shape

fps = 3
out = cv2.VideoWriter(
    "m_scale_conts.mp4",
    cv2.VideoWriter_fourcc(*"mp4v"),
    fps,
    (w, h)
)

for f in png_files:
    out.write(cv2.imread(f))

out.release()
