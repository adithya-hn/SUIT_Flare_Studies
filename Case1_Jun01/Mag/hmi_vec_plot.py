import numpy as np
import matplotlib.pyplot as plt
from sunpy.map import Map
from astropy.visualization import simple_norm

# Load the HMI SHARP maps
# br_map = Map("/Analysis/Research_Projects/HMI/Sharp/SHARP_CEA_Data/hmi.sharp_cea_720s.10799.20240219_001200_TAI/hmi.sharp_cea_720s.10799.20240219_001200_TAI.Br.fits")
# bt_map = Map("/Analysis/Research_Projects/HMI/Sharp/SHARP_CEA_Data/hmi.sharp_cea_720s.10799.20240219_001200_TAI/hmi.sharp_cea_720s.10799.20240219_001200_TAI.Bt.fits")
# bp_map = Map("/Analysis/Research_Projects/HMI/Sharp/SHARP_CEA_Data/hmi.sharp_cea_720s.10799.20240219_001200_TAI/hmi.sharp_cea_720s.10799.20240219_001200_TAI.Bp.fits")

br_map = Map("/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/hmi/SHARAP_data/13697/20240601_074800/hmi.sharp_cea_720s.11297.20240601_074800_TAI.Br.fits")
bt_map = Map("/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/hmi/SHARAP_data/13697/20240601_074800/hmi.sharp_cea_720s.11297.20240601_074800_TAI.Bt.fits")
bp_map = Map("/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/hmi/SHARAP_data/13697/20240601_074800/hmi.sharp_cea_720s.11297.20240601_074800_TAI.Bp.fits")

'''
# --- Extract data arrays ---
br = br_map.data
bt = bt_map.data
bp = bp_map.data

# --- Compute transverse field strength ---
b_trans = np.sqrt(bt**2 + bp**2)

# --- Downsample to avoid arrow crowding ---
step = 8
y, x = np.mgrid[0:br.shape[0]:step, 0:br.shape[1]:step]

# --- Vector components at sampled points ---
# In the CEA coordinate system:
#  - Bphi increases toward solar west → x direction
#  - Btheta increases toward south → need to flip sign for north-up plot
u = bp[::step, ::step]      # Bphi
v = -bt[::step, ::step]     # -Btheta (north-up)
mag = np.sqrt(u**2 + v**2)

# Normalize vectors to show direction but scale arrow length by strength
u_dir = u / mag
v_dir = v / mag

# --- Plot ---
fig = plt.figure(figsize=(9, 9))
ax = plt.subplot(projection=br_map)

# Background (Br)
norm = simple_norm(br, 'linear', percent=99)
im = ax.imshow(br, origin='lower', cmap='gray', norm=norm)

# Arrows: length ∝ |B_t|, direction = (Bphi, -Btheta)
q = ax.quiver(x, y, u_dir, v_dir, mag, scale=80, cmap='plasma', pivot='middle', width=0.001, alpha=0.8)

#ax.quiver(x, y, u, v, color='black', scale_units='xy', scale=1,angles='xy', headwidth=3, headlength=5)

# --- Colorbars & titles ---
cb1 = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb1.set_label(r'$B_r$ [G]', fontsize=12)
cb2 = plt.colorbar(q, ax=ax, fraction=0.046, pad=0.04)
#cb2.set_label(r'$B_t = \sqrt{B_\theta^2 + B_\phi^2}$ [G]', fontsize=12)

ax.set_title('HMI SHARP Vector Magnetic Field (Direction + Strength)', fontsize=13)
plt.tight_layout()
plt.show()'''


# --- Extract data ---
br = br_map.data
bt = bt_map.data
bp = bp_map.data

# --- Compute transverse field ---
b_trans = np.sqrt(bt**2 + bp**2)

# --- Downsample for clarity ---
step = 12
y, x = np.mgrid[0:br.shape[0]:step, 0:br.shape[1]:step]
u = bp[::step, ::step]      # Bphi → x direction
v = -bt[::step, ::step]     # -Btheta → y direction (north-up)
bmag = np.sqrt(u**2 + v**2)

# --- Normalize for arrow direction ---
u_dir = u / (bmag + 1e-6)
v_dir = v / (bmag + 1e-6)

# Scale factor controls arrow size
scale_factor = 0.008  # adjust this for longer/shorter arrows
u_scaled = u_dir * bmag * scale_factor
v_scaled = v_dir * bmag * scale_factor

# --- Plot ---
fig = plt.figure(figsize=(9, 9))
ax = plt.subplot(projection=br_map)

# Background: Br (radial field)
norm = simple_norm(br, 'linear', percent=99)
im = ax.imshow(br, origin='lower', cmap='RdBu_r', norm=norm)

# Arrows: length ∝ field strength, direction = (Bphi, -Btheta)
ax.quiver(x, y, u_scaled, v_scaled, color='k', angles='xy',
          scale_units='xy', scale=1, width=0.003, headwidth=4, headlength=6)

# Add colorbar and labels
cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label(r'$B_r$ [G]')

ax.set_title('HMI SHARP Vector Magnetic Field (Arrow length ∝ Strength)', fontsize=13)
plt.tight_layout()
plt.show()