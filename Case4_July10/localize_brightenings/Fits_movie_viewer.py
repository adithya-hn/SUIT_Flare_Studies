import sunpy.map as smap
import glob
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# -------------------------------------------------------
# USER INPUT
# -------------------------------------------------------
fits_files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/localize_brightenings/Diff_fits/NB04/*.fits"))     # list of FITS
maps = [smap.Map(f) for f in fits_files]
n_maps = len(maps)

# global index for current frame
index = 0
# -------------------------------------------------------

# compute common data range for sliders
data_min = min(m.data.min() for m in maps)
data_max = max(m.data.max() for m in maps)

initial_vmin = data_min
initial_vmax = 0.3 * data_max  # tweak as you like

#fig, ax = plt.subplots()
#plt.subplots_adjust(bottom=0.25)  # space for sliders

current_map = maps[index]
fig=plt.figure()
ax = fig.add_subplot(projection=maps[0])
im = current_map.plot(axes=ax, vmin=initial_vmin, vmax=initial_vmax)
ax.set_title(f"Frame {index+1}/{n_maps} {current_map.date} ")

# ---------------- Sliders for vmin / vmax ----------------
ax_vmin = plt.axes([0.15, 0.10, 0.7, 0.03])
ax_vmax = plt.axes([0.15, 0.05, 0.7, 0.03])

s_vmin = Slider(ax_vmin, 'vmin', data_min, data_max, valinit=initial_vmin)
s_vmax = Slider(ax_vmax, 'vmax', data_min, data_max, valinit=initial_vmax)

def update_scale(val):
    vmin = s_vmin.val
    vmax = s_vmax.val
    if vmin >= vmax:  # simple safety check
        return
    im.set_clim(vmin, vmax)
    fig.canvas.draw_idle()

s_vmin.on_changed(update_scale)
s_vmax.on_changed(update_scale)

# --------------- Keyboard navigation to move sequence ---------------
def on_key(event):
    global index, im
    if event.key == 'right':
        index = (index + 1) % n_maps
    elif event.key == 'left':
        index = (index - 1) % n_maps
    else:
        return

    #ax.cla()  # clear axes
    fig.clear()                                       # clear figure
    ax = plt.subplot(projection=maps[index])          # projection here
    im = maps[index].plot(axes=ax, 
                         vmin=s_vmin.val, vmax=s_vmax.val)
    ax.set_title(f"[{index+1}/{n_maps}]  { maps[index].date} ")
    plt.colorbar()
    fig.canvas.draw_idle()


fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()