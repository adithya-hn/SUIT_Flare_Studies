import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os
from matplotlib.patches import Rectangle
import math
valid_D=fits.getdata('NB04_full_first_diff_cube.fits')
cube=fits.getdata('NB04_full_aligned_cube.fits')

sigmaMap = 1.4826 * np.median(np.abs(cube - np.median(cube, axis=0)), axis=0)

plt.title('Aligned image MAD sigma along time axis')
plt.imshow(sigmaMap,origin='lower',vmin=100,vmax=4000)
plt.colorbar()
plt.savefig('aligned_image_sigma_temp_vmin_max.png')
plt.show()

# median_I = np.median(cube, axis=0)
# p10, p50, p90 = np.percentile(median_I, [10, 50, 90])

# def pick_pixel(target):
#     idx = np.argmin(np.abs(median_I - target))
#     return np.unravel_index(idx, median_I.shape)

# pixels = [
#     pick_pixel(p10),  # quiet
#     pick_pixel(p50),  # intermediate
#     pick_pixel(p90),  # bright
# ]
# print(pixels)
os._exit(0)

# example pixels
def plot_cutout(img, y, x, halfsize=20):
    cut = img[y-halfsize:y+halfsize+1,
              x-halfsize:x+halfsize+1]

    plt.figure(figsize=(3,3))
    plt.imshow(cut, origin='lower', cmap='gray')
    plt.plot(halfsize, halfsize, 'ro')
    plt.title("Local cutout")
    plt.savefig('local_cutout.png',dpi=300)
    plt.show()

def plot_slices(img, y, x, label):
    fig, axs = plt.subplots(1, 3, figsize=(12,3))

    axs[0].imshow(img, origin='lower', cmap='gray')
    axs[0].plot(x, y, 'ro')
    axs[0].set_title(label)

    axs[1].plot(img[y, :])
    axs[1].axvline(x, color='r')
    axs[1].set_title("Row slice")

    axs[2].plot(img[:, x])
    axs[2].axvline(y, color='r')
    axs[2].set_title("Column slice")

    plt.tight_layout()
    plt.show()




pixels = [(56,356), (101,311), (242,293),(260,230),(300,100),(460,280),(270,180),(270,181),(271,180)]

img = cube[0]

plt.figure(figsize=(6,6))
plt.imshow(img, origin='lower', cmap='gray')
for i, (y, x) in enumerate (pixels):
    plt.plot(x, y, 'ro', markersize=5)
    plt.text(x+5, y+5,f'{i} ' ,color='red', fontsize=8)
plt.title("Selected reference pixels")
plt.savefig('selected pixels.png',dpi=300)
plt.show()


bin_width = 10


# for i, (y, x) in enumerate(pixels):
#     d = valid_D[:, y, x]     # Haar differences at this pixel
#     d = d[np.isfinite(d)]

#     bins = np.arange(
#         np.floor(d.min()/bin_width)*bin_width,
#         np.ceil(d.max()/bin_width)*bin_width + bin_width,
#         bin_width
#     )
    
#     plt.subplot(1, len(pixels), i+1)
#     plt.hist(d, bins=bins, histtype='step', linewidth=1.5)
#     plt.axvline(0, color='k', lw=0.5)
#     plt.title(f"Pixel ({y},{x})")
#     plt.xlabel("ΔI (counts)")
#     if i == 0:
#         plt.ylabel("N")
n = len(pixels)                 # number of selected pixels
ncols = 3                       # choose (e.g. 3)
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(4*ncols, 3*nrows),
    squeeze=False
)
for i, ((y, x), ax) in enumerate(zip(pixels, axes.flat)):
    d = valid_D[:, y, x]
    d = d[np.isfinite(d)]

    bins = np.arange(
        np.floor(d.min()/bin_width)*bin_width,
        np.ceil(d.max()/bin_width)*bin_width + bin_width,
        bin_width
    )

    ax.hist(d, bins=bins, histtype='step', linewidth=1.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_title(f"{i} Pixel ({y},{x})")
    ax.set_xlabel("ΔI (counts)")
    ax.set_ylabel("N")
plt.tight_layout()
plt.savefig('asorted_pixels_first_difference_distrb.png',dpi=300)
plt.show()

# robust sigma per pixel
sigma_map = 1.4826 * np.median(np.abs(valid_D - np.median(valid_D, axis=0)), axis=0)

# zoom into a small region
ys, xs = slice(250,290), slice(160,200)
plt.figure(figsize=(6,6))
plt.imshow(img, origin='lower', cmap='gray')

# rectangle parameters
y0, y1 = ys.start, ys.stop
x0, x1 = xs.start, xs.stop

rect = Rectangle(
    (x0, y0),              # bottom-left corner
    x1 - x0,               # width
    y1 - y0,               # height
    edgecolor='red',
    facecolor='none',
    linewidth=2
)

plt.gca().add_patch(rect)
plt.title("Selected slice region")
plt.savefig('selected region.png',dpi=300)
plt.show()


plt.figure(figsize=(6,6))
plt.imshow(sigma_map, origin='lower', cmap='gray')

# rectangle parameters
y0, y1 = ys.start, ys.stop
x0, x1 = xs.start, xs.stop

rect = Rectangle(
    (x0, y0),              # bottom-left corner
    x1 - x0,               # width
    y1 - y0,               # height
    edgecolor='red',
    facecolor='none',
    linewidth=2
)

plt.gca().add_patch(rect)
plt.title("Selected slice region")
plt.savefig('selected region on sigma map.png',dpi=300)
plt.show()


plt.figure(figsize=(5,4))
plt.imshow(sigma_map[ys, xs], origin='lower')
plt.plot(20,20,'ro',markersize=2)
plt.text(20+2, 20+2,f'6 ' ,color='red', fontsize=8)

plt.plot(20,21,'bo',markersize=2)
plt.text(20, 20+3,f'7 ' ,color='blue', fontsize=8)

plt.plot(21,20,'ko',markersize=2)
plt.text(21+4, 20,f'8 ' ,color='k', fontsize=8)

plt.colorbar(label="σ (counts)")
plt.title("Local noise σ")
plt.savefig('Patch_sigm_map.png',dpi=300)
plt.show()

median_I = np.median(cube, axis=0)
mean_I = np.mean(cube, axis=0)
plt.figure(figsize=(10,6))
plt.scatter(median_I.flatten(), sigma_map.flatten(), s=1, alpha=0.3)
plt.xlabel("Median intensity (counts)")
plt.ylabel("Noise σ (counts)")
plt.title("Noise vs brightness")
plt.savefig('median_vs_sig.png',dpi=300)
plt.close()

plt.figure(figsize=(10,6))
plt.scatter(mean_I.flatten(), sigma_map.flatten(), s=1, alpha=0.3)
plt.xlabel("Mean intensity (counts)")
plt.ylabel("Noise σ (counts)")
plt.title("Noise vs brightness")
plt.savefig('mean_vs_sig.png',dpi=300)
plt.close()

# coefficient of variation
cv = np.std(sigma_map) / np.mean(sigma_map)
print("σ spatial CV:", cv)
