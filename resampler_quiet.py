#!/usr/bin/env python

#!/usr/local/Python-3.6.1/python

# This program subtracts two fits images of the same size and
# can put the output into a file. 
#
# It also pops up an image of the subtracted images.
#
# It will cycle through images in a multifits file.
#
# It will resample the image on pressing the 'l' key.
# It will resample the image to 5" pixel on pressing 5.

# Usage:
#  ./resampler.py Ha.fit r.fit diff.fit
#
# rebin subroutine from
# https://stackoverflow.com/questions/384759/how-to-convert-a-pil-image-into-a-numpy-array
#
# Need to following FITS keywords to calculate pixel size.
# From the main header
# HIERARCH ESO TEL FOCU SCALE        14.267 / Focal scale (arcsec                         mm)
#   
# From a subimage header:
# HIERARCH ESO DET CHIP PSZX   15.0                    Size of pixel in X
# HIERARCH ESO DET CHIP PSZY   15.0                     Size of pixel in Y
# Assume size is in micrometers, um so that gives 0.214005 arcseconds/pixel

# Load packages
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from astropy.io import fits

# Need this for accessing command line parameters
import sys

# Images are stored as numpy arrays
import numpy as np

# Need these for coordinate conversion
from astropy import units as u
from astropy.coordinates import SkyCoord

from mpl_toolkits.mplot3d import Axes3D

# For resampling
import PIL
from PIL import Image, ImageOps

from image_registration import chi2_shift
from image_registration import register_images
from image_registration.fft_tools import shift

from skimage import exposure

from astropy import wcs
import os

#################################################################
# input parameters are difference array, new shape.
# Want to return an array.

# rebin

# subroutine from
# https://stackoverflow.com/questions/384759/how-to-convert-a-pil-image-into-a-numpy-array


#def rebin(a, *args):
#    shape = a.shape

#    print('Rebinning array')

#    lenShape = len(shape)
#    factor = asarray(shape)/asarray(args)
#    evList = ['a.reshape('] + \
#             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
#             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
#    print( ''.join(evList))
#   return eval(''.join(evList))


#########################################################################
def rebin(numpy_array_image, boxsize):
    # convert nympy array image to PIL.Image
    # ratio is new_size/old_size eg 0.5
 
    print('Rebinning array using PIL image.resize')
    image = Image.fromarray(numpy_array_image)
    old_width = float(image.size[0])
    old_height = float(image.size[1])
    #    ratio = float( new_height / old_height)
    new_width = int(old_width / boxsize)
    new_height = int(old_height / boxsize)
#    image = image.resize((new_width, new_height), PIL.Image.ANTIALIAS)
    image = image.resize((new_width, new_height), PIL.Image.BILINEAR)
    # convert PIL.Image into nympy array back again
    return np.array(image)


###################################################################

def block_average(numpy_array_image, boxsize):

# never got this routine to work, don't use!!
    print('Block averaging array using np.reshape -- does not work!!')
    arrsize = numpy_array_image.shape
#    print("Size of array is", arrsize)
#    boxsize = ratio
#    print("boxsize is ",boxsize)
    new_dimensions = [x/boxsize for x in arrsize]
    new_xdim=int(arrsize[0]/boxsize)
    new_ydim=int(arrsize[1]/boxsize)
#    print("new dimensions are: ",new_xdim, new_ydim)
    # Clever stuff with numpy below                   
    temparr = numpy_array_image.reshape( (new_xdim, boxsize*boxsize, new_ydim) )
#    print(temparr)  
#    blk_av_image =  np.mean(temparr, axis=1)
#    print("new image\n",blk_av_image)
    return np.array(blk_av_image)

###############################################################
def block_average2(numpy_array_image, boxsize):

# Use loops to calculate block mean. Probably slow.

    print('Rebinning array using loops in python')
    arrsize = numpy_array_image.shape
    print("Size of array is", arrsize)
#    boxsize = ratio
    old_xdim = arrsize[0]
    old_ydim = arrsize[1]
    new_xdim = int(old_xdim/boxsize)
    new_ydim = int(old_ydim/boxsize)

    newarr = np.empty( (new_xdim, new_ydim) )
    print("New dimensions will be: ",new_xdim, new_ydim)
    for i in range(0,new_xdim):
        for j in range (0,new_ydim):
            block_sum = 0.0
            kstart = (i*boxsize)
            kend = (i*boxsize) + boxsize
            lstart = (j*boxsize)
            lend = (j*boxsize) + boxsize
#            print("kstart, kend-1:", kstart, kend-1)
#            print("lstart, lend-1:", lstart, lend-1)
            for k in range( kstart, kend ):
                for l in range( lstart, lend ):
#                    print("i,j,k,l:",i,j,k,l)
                    block_sum = block_sum + numpy_array_image[k,l]
#                    print("block_sum is ",block_sum)
#                    print("element "


            newarr[i,j] = block_sum/(boxsize * boxsize)
                    
    return newarr

#########################
#                       #
#   Start of Program.   #
#                       # 
#########################

# If we don't have 3 arguments, exit with warning message
if ( len(sys.argv) > 4 or len(sys.argv) < 4 ):
    print("There is a problem with your arguments.")
    print("Argument1 should be a multifits file.")
    sys.exit('Parameters are wrong')

#print("\nImage 1")

# Open image 1
image1 = sys.argv[1]
print("opening ",image1)
hdu_list1 = fits.open(image1)

filter_name = 'Unknown'
# Get filter name from FITS header
try:
    filter_name = hdu_list1[0].header['HIERARCH ESO INS FILT1 NAME']
except KeyError:
    print("Oops FITS header HIERARCH ESO INS FILT1 NAME does not exist")
#    filter_name = 'Unkown filter '

# Try looking for a FILTER keyword
try:
    filter_name = hdu_list1[0].header['FILTER']
except KeyError:
    print("Oops FITS header FILTER does not exist")
#    filter_name = 'Unknown filter'

print("Filter name is:",filter_name )

# How many subimages do we have?
try:
    nimages = hdu_list1[0].header['HIERARCH ESO DET CHIPS']
except KeyError:
    print("Oops FITS header HIERARCH ESO DET CHIPS does not exist, setting nimages to 1")
    nimages = 1

print("number of images is ",nimages,"\n")

# Get coordinates of centre of first image whole big image
ra_whole_image = '99 99 99'

try:
    ra_whole_image = hdu_list1[0].header['RA']
except KeyError:
    print("Oops FITS header RA does not exist")
#    ra_whole_image = '99 99 99'

try:
    ra_whole_image = hdu_list1[0].header['OBJCTRA']
except KeyError:
    print("Oops FITS header OBJCTRA does not exist")
#    ra_whole_image = '99 99 99'
print("RA of image centre is ",ra_whole_image)

# Now do the Dec
dec_whole_image = '99 99 99'
try:
    dec_whole_image = hdu_list1[0].header['DEC']
except KeyError:
    print("Oops FITS header DEC does not exist")

try:
    dec_whole_image = hdu_list1[0].header['OBJCTDEC']
except KeyError:
    print("Oops FITS header OBJCTDEC does not exist")

#print("Dec of image centre is ",dec_whole_image)


print("Coordiantes of image1 centre are ",ra_whole_image, dec_whole_image)

# Problem here!
# Convert ra, dec in degrees to hms, dms
#c = SkyCoord(ra=ra_whole_image*u.degree, dec=dec_whole_image*u.degree)
c = SkyCoord(ra_whole_image, dec_whole_image, unit=(u.hourangle, u.deg))
hmsdms_coords = c.to_string('hmsdms')

#print("Coordinates of image1 centre are ",hmsdms_coords)
print("\n")

#######################################################################

# Open image 2
image2 = sys.argv[2]
print("opening ",image2)
hdu_list2 = fits.open(image2)

# Get filter name from FITS header

filter_name2 = 'Unkown'
try:
    filter_name2 = hdu_list2[0].header['HIERARCH ESO INS FILT1 NAME']
except:
    print("Oops FITS header HIERARCH ESO INS FILT1 NAME does not exist")

# Try looking for a FILTER keyword
try:
    filter_name2 = hdu_list2[0].header['FILTER']
except KeyError:
    print("Oops FITS header FILTER does not exist")

print("Image 2 Filter is ",filter_name2)

# How many subimages do we have in the second image?
try:
    nimages2 = hdu_list2[0].header['HIERARCH ESO DET CHIPS']
except KeyError:
    print("Oops FITS header HIERARCH ESO DET CHIPS does not exist, setting nimages to 1")
    nimages2 = 1

print("number of images in file 2 is ",nimages2,"\n")

# Get coordinates of centre of second whole big image
ra_whole_image2 = '99 99 99'

try:
    ra_whole_image2 = hdu_list2[0].header['RA']
except KeyError:
    print("Oops FITS header2 RA does not exist")

try:
    ra_whole_image2 = hdu_list2[0].header['OBJCTRA']
except KeyError:
    print("Oops FITS header2 OBJCTRA does not exist")

#print("RA of image2 centre is ",ra_whole_image2)

# Now do the Dec
dec_whole_image = '99 99 99'

try:
    dec_whole_image2 = hdu_list2[0].header['DEC']
except KeyError:
        print("Oops FITS header DEC does not exist")

try:
    dec_whole_image2 = hdu_list2[0].header['OBJCTDEC']
except KeyError:
    print("Oops FITS header OBJCTDEC does not exist")

#print("Dec of image centre is ",dec_whole_image2)


print("Coordiantes of image2 centre are ",ra_whole_image2, dec_whole_image2)

# Convert ra, dec in degrees to hms, dms
c2 = SkyCoord(ra_whole_image2, dec_whole_image2,unit=(u.hourangle, u.deg))
hmsdms_coords2 = c2.to_string('hmsdms')

print("Second image:Coordinates of image2 centre are ",hmsdms_coords2)

# Check images are aligned
#print("Got to here 1")
centre_seperation = c.separation(c2)
#print("Got to here 2")

print ("\nCentre separation is ",centre_seperation.arcsecond,' arcseconds')

if (centre_seperation.arcsecond > 40.0 ):
    print("\nimages are not aligned\n")
    sys.exit("Quitting")

# How many subimages do we have?
#nimages2 = hdu_list1[0].header['HIERARCH ESO DET CHIPS']
#print("Second image: number of images is ",nimages2)

if (nimages != nimages2):
    sys.exit('Files have different number of subimages. Quitting...')

# Get the plate scale
f_platescale = 999
try:
    f_platescale = hdu_list1[0].header['HIERARCH ESO TEL FOCU SCALE']
except KeyError:
    print("keyword HIERARCH ESO TEL FOCU SCALE does not exist.")

try:
    f_platescale = hdu_list1[0].header['PLTSCALE']
except KeyError:
    print("keyword PLTSCALE does not exist.")

    print("Focal platescale is ",f_platescale)

# Now need the physical size of each pixel in microns
pixel_size_nm = 999
try:
    pixel_size_nm = hdu_list1[1].header['HIERARCH ESO DET CHIP PSZX']
except KeyError:
    print("keyword HIERARCH ESO DET CHIP PSZX does not exist.")

try:
    pixel_size_nm = hdu_list1[1].header['XPIXELSZ']
except KeyError:
    print("XPIXELSZ header does not exist.")

if (pixel_size_nm != 999):
    print("Size of pixel in nm",pixel_size_nm)
    pixel_size_arcsec = f_platescale*(pixel_size_nm/1000.0)
    print("Pixel size is ",pixel_size_arcsec, "arcsec.")

# This is crucial, interactive mode on
#plt.ion()

# Get WCS info
wcs_info = wcs.WCS(hdu_list1[1].header)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection=wcs_info)

# Instructions on program use.
#print("\n")
#print("Instructions:")
#print("q to quit")
#print("r to resample with PIL image.rebin.")
#print("b to resample with python loops")
#print("w to write fits image.")
#print("return to move on.")
#print("\n")

#for x in range(1,nimages):
x = 1
while x <= nimages:
    print("Extension is ",x)
    image_data1 = fits.getdata(sys.argv[1], ext=x)
    image_data2 = fits.getdata(sys.argv[2], ext=x)
    # subtract two numpy arrays
    #    diff_arr = np.subtract(image_data1, image_data2)
    xoff, yoff = register_images(image_data1, image_data2, zeromean='true', maxoff=30)
    print("Offsets are ",xoff, yoff)
    shifted_image_data2 = shift.shiftnd(image_data2, (-yoff, -xoff))
    print("Subtraction is Ha - 0.4*r \n")
    diff_arr = np.subtract(image_data1, 0.4*shifted_image_data2)

    # Remove median of image to get rid of sky background
    background = np.median(diff_arr)
    diff_arr = diff_arr - background
    ########################################
    # Write out subtracted image

    # Get WCS from first image
    wcs_info = wcs.WCS(hdu_list1[x].header)
    #    print("wcs is",w)
    wcs_head = wcs_info.to_header()
    
    # we create a PrimaryHDU object to encapsulate the data:
    hdu = fits.PrimaryHDU(diff_arr,header=wcs_head)
    
    # Create output file name. Add fits suffix if not in output filename.
    outfilename = sys.argv[3]
    containsfits = False
    # Tried to do this with loops and it got complicated fast especially
    # as matching .fits also means .fit matches
    if (".fits" in outfilename):
        basename=os.path.splitext(outfilename)[0]
        outfilename = basename+str(x)+".fits"
        containsfits = True

    if (".fit" in outfilename and containsfits == False):
        basename=os.path.splitext(outfilename)[0]
        outfilename = basename+str(x)+".fit"
        containsfits = True

    if (".fts" in outfilename and containsfits == False):
        basename=os.path.splitext(outfilename)[0]
        outfilename = basename+str(x)+".fts"
        containsfits = True

    # We can't write compressed images at the moment so change to fits
    if (".fits.fz" in outfilename and containsfits == False):
        basename=os.path.splitext(outfilename)[0]
        outfilename = basename+str(x)+".fts"
        containsfits = True

    if (containsfits == False):
        outfilename = outfilename+str(x)+".fits"
    

    print("output file name will be ",outfilename)

    hdu.writeto(outfilename,overwrite=True)

    #########################################


#    hdu = fits.PrimaryHDU(latest_binned_diff_arr)
#    hdu.writeto('resampled_r.fits', overwrite=True)
                #                continue

#                keyval = input("Press Enter to continue...")


#    elif (keyval == 'b'):
#        print("Block averaging image.")
#        current_diff_arr = diff_arr
#        current_pixel_size = pixel_size_arcsec
        
#        while(keyval == 'b'):        
#            # while keypress is b, keep resampling by factor of 2
#            boxsize = 2
#            #        rebin array
#
#            latest_binned_diff_arr = block_average2(current_diff_arr,boxsize)
#            print("After block av: file shape is",latest_binned_diff_arr.shape)
#            current_pixel_size = current_pixel_size* boxsize
#            print("After block av, pixel size is", current_pixel_size," arc seconds.")

            #        Plot image before rebinning on left
#            fig.clf()
#            fig.add_subplot(121,projection=wcs_info)

            #diff_arr_histeq = np.sort(np.ravel(diff_arr)).searchsorted(diff_arr)
#            diff_arr_histeq = exposure.equalize_hist(diff_arr)
#            lowval = np.percentile(diff_arr_histeq, 10)
#            highval = np.percentile(diff_arr_histeq, 90)
#            plt.imshow(diff_arr_histeq, cmap='hot',vmin=lowval, vmax=highval, aspect='equal',interpolation='bilinear',origin='lower')

            # Plot image after rebinning on right
#            fig.add_subplot(122,projection=wcs_info)
#            fig.add_subplot(122)
            #        plt.clf()
            #arr_histeq = np.sort(np.ravel(latest_binned_diff_arr)).searchsorted(latest_binned_diff_arr)
#            arr_histeq = exposure.equalize_hist(latest_binned_diff_arr)
#            lowval = np.percentile(arr_histeq, 10)
#            highval = np.percentile(arr_histeq, 90)

            # Display second image
#            plt.imshow(arr_histeq, cmap='hot',vmin=lowval, vmax=highval, aspect='equal',interpolation='bilinear',origin='lower')
#            current_diff_arr = latest_binned_diff_arr
#            keyval = input("Press Enter to continue...")
#            if keyval == 'w':
#                print('Writing resampled differene image to fits file resampled_b.fits')
#	        # write latest_binned_diff_arr to fits file
#                hdu = fits.PrimaryHDU(latest_binned_diff_arr)
#                hdu.writeto('resampled_b.fits', overwrite=True)
#                #                continue
#
#                keyval = input("Press Enter to continue...")


    # elif (keyval == '5'):
    #     print('Resampling array to 5 arcsec pixels...')
    #     current_diff_arr = diff_arr

    #     while(keyval == '5'):
    #         # while keypress i r, keep resampling by factor of 2
    #         boxsize = 25
    #         #        rebin array
    #         latest_binned_diff_arr = rebin(current_diff_arr,boxsize)
    #         #        Plot image before rebinning on left
    #         fig.add_subplot(121, projection=wcs_info)
    #         lowval = np.percentile(diff_arr, 10)
    #         highval = np.percentile(diff_arr, 90)
    #         plt.imshow(diff_arr, cmap='hot',vmin=lowval, vmax=highval, aspect='equal',interpolation='bilinear',origin='lower')
    #         # Plot image after rebinning on right
    #         fig.add_subplot(122, projection=wcs_info)
    #         #        plt.clf()
    #         lowval = np.percentile(latest_binned_diff_arr, 10)
    #         highval = np.percentile(latest_binned_diff_arr, 90)
    #         plt.imshow(latest_binned_diff_arr, cmap='hot',vmin=lowval, vmax=highval, aspect='equal',interpolation='bilinear',origin='lower')
    #         current_diff_arr = latest_binned_diff_arr
    #         keyval = input("Press Enter to continue...")

#    elif (keyval == 'w'):
#        print('Writing difference image to fits file difference.fits')
#        # write latest_binned_diff_arr to fits file
#        hdu = fits.PrimaryHDU(latest_binned_diff_arr)
#        hdu.writeto('difference.fits', overwrite=True)


    x = x + 1            
    continue
    hdu_list.close




