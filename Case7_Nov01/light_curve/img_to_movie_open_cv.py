import cv2
import glob

png_files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aligned_crop_pngs/nb4_png/*.png"))

img = cv2.imread(png_files[0])
h, w, _ = img.shape

fps = 10
out = cv2.VideoWriter(
    "output.mp4",
    cv2.VideoWriter_fourcc(*"mp4v"),
    fps,
    (w, h)
)

for f in png_files:
    out.write(cv2.imread(f))

out.release()
