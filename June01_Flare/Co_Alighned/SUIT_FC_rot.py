
import numpy as np
from scipy.interpolate import interp2d
import math
from PIL import Image
import math
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import sunpy.map
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord


def create_subpixels(matrix, scale):
  new_matrix = np.zeros((matrix.shape[0] * scale, matrix.shape[1] * scale))
  for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
      new_matrix[i * scale:(i + 1) * scale, j * scale:(j + 1) * scale] = matrix[i, j]
  return (new_matrix/pow(scale,2))

def bin_back(matrix,scale):
   new_matrix = np.zeros((int(matrix.shape[0] / scale), int(matrix.shape[1] / scale)))
   for i in range(new_matrix.shape[0]):
    for j in range(new_matrix.shape[1]):
       new_matrix[i, j]=np.sum(matrix[i * scale:(i + 1) * scale, j * scale:(j + 1) * scale])
   return new_matrix



#img.meta.update({'CROTA2':0})
#img.save('SUT_T24_0897_000441_Lev1.0_2024-06-25T14.56.07.464_0973NB03.fits')

img_=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Co_Alighned/Sample/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.22.42.472_0973NB03.fits')
hmi_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/SHARAP_data/13697/20240601_072400/hmi.sharp_cea_720s.11297.20240601_072400_TAI.magnetogram.fits')
hmi_map.peek()
suit_pos = get_horizons_coord(-21, img_.date)
print(suit_pos.rsun)
img_.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
img_.save('original.fits',overwrite=True)

img=sunpy.map.Map('original.fits')
Data=img.data
suit_header=img.fits_header
w = WCS(suit_header)

#------------------
'''
proj_map=img_.reproject_to(hmi_map.wcs)
fig = plt.figure(figsize=(5, 5))
ax1= fig.add_subplot(projection=proj_map)
ax1.imshow(proj_map.data,origin='lower',cmap='gray')
ax1.set_title('original image')
ax1.imshow(hmi_map.data,origin='lower',cmap='gray',alpha=0.6)
plt.show()
'''
#------------------

pscl=suit_header['CDELT1'] #plate scale
suit_header['CROTA2']=0
X=suit_header['CRVAL1']/pscl #ROI centre in pixels
Y=suit_header['CRVAL2']/pscl

crx=suit_header['CRPIX1']
cry=suit_header['CRPIX2']

x1=suit_header['X1']
x2=suit_header['X2']
y1=suit_header['Y1']
y2=suit_header['Y2']

ccdX=suit_header['NAXIS1']
ccdY=suit_header['NAXIS2']

#img3=np.zeros((4096,4096))
#img3[y1:(y1+704),x1:(x1+704)]=Data
imgDim=Data.shape
if imgDim[0]==4096:
   img_cenX=crx
   img_cenY=cry
else:
    img_cenX=x1+ccdX/2-X 
    img_cenY=y1+ccdY/2-Y
#--------------------------------

#print('-->',img_cenX,img_cenY)


angle=int(-7.6) #-6 to -10 is the range 


#--------------------------------
angle_rad=math.radians(angle)      #converting degrees to radians
cosine=math.cos(angle_rad)
sine=math.sin(angle_rad)



i=0
j=0
x=j-img_cenX-1
y=i-img_cenY-1
y0=(-x*sine+y*cosine)
x0=(x*cosine+y*sine)
cy0=np.floor((y0)+img_cenY-1)
cx0=np.floor((x0)+img_cenX-1)
print('origin>>',cx0,cy0)

i=ccdY
j=0
x=j-img_cenX-1
y=i-img_cenY-1
y0=(-x*sine+y*cosine)
x0=(x*cosine+y*sine)
cy1=np.floor((y0)+img_cenY-1)
cx1=np.floor((x0)+img_cenX-1)
print('TL>>',cx1,cy1)

i=0
j=ccdY
x=j-img_cenX-1
y=i-img_cenY-1
y0=(-x*sine+y*cosine)
x0=(x*cosine+y*sine)
cy3=np.floor((y0)+img_cenY-1)
cx3=np.floor((x0)+img_cenX-1)
print('BR>>',cx3,cy3)

i=ccdY
j=ccdX
x=j-img_cenX-1
y=i-img_cenY-1
y0=(-x*sine+y*cosine)
x0=(x*cosine+y*sine)
cy2=np.floor((y0)+img_cenY-1)
cx2=np.floor((x0)+img_cenX-1)
print('TR>>',cx2,cy2)

#midpix
i=ccdY/2
j=ccdX/2

x=j-img_cenX-1
y=i-img_cenY-1
y0=(-x*sine+y*cosine)
x0=(x*cosine+y*sine)



#print(cx2-cx0,cy1-cy3)
#print(cx3-cx1,cy2-cy0)
scale=1

if angle>0:
    new_width=(cx2-cx0)*scale   #4096*scale#
    new_height=(cy1-cy3)*scale  #4096*scale#)
    cy=np.floor((y0)+img_cenY-1)-cy3
    cx=np.floor((x0)+img_cenX-1)-cx0
    print('Positive Angle:',cx0,cy3,cx,cy)

if angle<0:
    new_width=(cx3-cx1)*scale #np.floor(new_X2-new_X1)+1+(img_cenY-1)
    new_height=(cy2-cy0)*scale
    print('diff_width : ',640-cx3-cx1, 'Diff_height :',640-cy2-cy0  )
    cy=np.floor((y0)+img_cenY-1)-cy3
    cx=np.floor((x0)+img_cenX-1)-cx0
    print('Negative Angle:',cx1,cy0,cx,cy)

print('cenRot>>',cx,cy)
sky = w.pixel_to_world(cx,cy)
print(sky.Tx.arcsec,sky.Ty.arcsec)

if imgDim[0]==4096:
    suit_header['CRVAL1']=0
    suit_header['CRVAL2']=0
    if angle<0:
        suit_header['CRPIX2']=img_cenY-cy0
        suit_header['CRPIX1']=img_cenX-cx1
    if angle>0:
        suit_header['CRPIX2']=img_cenY-cy3
        suit_header['CRPIX1']=img_cenX-cx0
else:  
    suit_header['CRVAL1']=sky.Tx.arcsec
    suit_header['CRVAL2']=sky.Ty.arcsec
    suit_header['CRPIX1']=new_width/(2*scale)
    suit_header['CRPIX2']=new_height/(2*scale)
#print(new_X1+img_cenX,new_Y1+img_cenY,new_X2+img_cenX,new_Y2+img_cenY)
print(new_height,new_width)


Rotated_img=np.zeros((int(new_height),int(new_width)))
subPixMap=create_subpixels(Data,scale)

for i in range(ccdY*scale):
    for j in range(ccdX*scale):
        x=j-img_cenX*scale-1
        y=i-img_cenY*scale-1


        y0=(-x*sine+y*cosine)
        x0=(x*cosine+y*sine)
       
        if angle>0:
            new_y=int((y0)+img_cenY*scale-1)-int(cy3*scale)#+y1*scale
            new_x=int((x0)+img_cenX*scale-1)-int(cx0*scale)#+x1*scale

        if angle<0:
            
            new_y=int((y0)+img_cenY*scale-1)-int(cy0*scale)#+y1*scale
            new_x=int((x0)+img_cenX*scale-1)-int(cx1*scale)#+x1*scale

        if 0 <= new_x < new_width and 0 <= new_y < new_height :
            Rotated_img[new_y,new_x]+=subPixMap[i,j]


Rotated_img=bin_back(Rotated_img,scale)

hdl=fits.PrimaryHDU(Rotated_img)
hdl.header=suit_header
hdl.writeto('Rotated.fits',overwrite=True)
print('Error:',(np.sum(Rotated_img)-np.sum(Data))/np.sum(Data)*100)

fig = plt.figure(figsize=(10, 5))
ax1, ax2= fig.subplots(1, 2)
ax1.imshow(Data,origin='lower',cmap='gray')
ax1.set_title('original image')
#ax2.imshow(Data2,origin='lower',cmap='gray',alpha=0.6)
ax2.imshow(Rotated_img,origin='lower',cmap='gray')
#ax2.imshow(img3,origin='lower',cmap='gray',alpha=0.4)
ax2.set_title('Rotated by {} deg:'.format(angle))
plt.savefig('Rotated_ROI_Maps.png')
plt.show()