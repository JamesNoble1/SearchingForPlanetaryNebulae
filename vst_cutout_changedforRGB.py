#!/usr/bin/env python

import matplotlib.pyplot as plt
from astropy.io import fits
# Need this for accessing command line parameters
import sys
# Images are stored as numpy arrays
import numpy as np
# Need these for coordinate conversion
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.nddata.utils import Cutout2D
import os

####################################################################
# Usage
# vst_cutout.py file.fits rah ram ras decdeg decmin decsecs boxsize
#
# Test with
# ./vst_cutout.py ADP.2015-05-11T10_19_31.813.fits.fz 17 40 52.0 -26 30 50.0 50
# coordinates should be in first image, boxsize 50 x 50
#
# Find the PN
# ./vst_cutout.py ADP.2015-05-11T10_19_31.813.fits.fz 17 43 57.4 -26 11 52.0 100
# boxsize is 100 x 100 pixels
###################################################################

# If we don't have 8 arguments, exit with warning message
'''if ( len(sys.argv) != 9 ):
    print("There is a problem with your arguments.")
    print("Argument1 should be a multifits file followed by coords and box size.")
    sys.exit('Parameters are wrong')'''

sys.argv = [0,'G:/RGB_image/Fits',17,48,25,-22,11,53,300]
os.chdir(sys.argv[1])
for f in os.listdir(sys.argv[1]):
    
    image1 = os.path.abspath(f)
    print("opening ",image1)
    hdu_list1 = fits.open(image1)
    f_name, f_ext = os.path.splitext(f)
    # How many subimages do we have?
    try:
        nimages = hdu_list1[0].header['HIERARCH ESO DET CHIPS']
    except KeyError:
        print("Oops FITS header HIERARCH ESO DET CHIPS does not exist, setting nimages to 1")
        nimages = 1
    
    print("number of images is ",nimages,"\n")
    
    # Use contained_by(wcs[,image[) from
    # https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    
    # Assign variables to input parameters
    ra_hours = str(sys.argv[2])
    ra_mins = str(sys.argv[3])
    ra_secs = str(sys.argv[4])
    ra_str = ra_hours + " " + ra_mins + " " + ra_secs
    
    print("RA string is ",ra_str)
    
    dec_degrees = str(sys.argv[5])
    dec_mins = str(sys.argv[6])
    dec_secs = str(sys.argv[7])
    dec_str = dec_degrees + " " + dec_mins + " " + dec_secs
    
    print("Dec string is ",dec_str)
    
    # Now create skycoord class
    obj_coord = SkyCoord(ra_str, dec_str, frame='icrs',unit=(u.hourangle,u.deg) )
    print("RA in degrees is ",obj_coord.ra.degree)
    print("Dec in degrees is ",obj_coord.dec.degree)
    
    # Get header data for whole image and subimages
    header_list = fits.open(image1)
    
    
    isinimage = 'False'
    # Loop over subimages
    x = 1
    while x <= nimages:
    #    print("Extension is ",x)
    #    image_data = fits.getdata(sys.argv[1], ext=x)
        wcs_info = wcs.WCS(header_list[x].header)
    
        isinimage = obj_coord.contained_by(wcs_info)
        if isinimage == True:
            print("These coordinates are in subimage ",x)
    #        print("The type of side is ",type(sys.argv[8]) )
    #  Extract a cutout box
            side = int(sys.argv[8])
            boxsize = [side, side ]
    #        print("box size is ",boxsize)
    #  Get the data for the subimage containing our coordinates.
            image_data = fits.getdata(image1, ext=x)
    #  Extract the data        
            extracted_data = Cutout2D(image_data, obj_coord, boxsize, wcs=wcs_info )
    #  Create a popup figure window
            fig=plt.figure(figsize=(6,6) )
    #  set up axes using wcs (world coordinate system) of extracted data        
            ax = fig.add_subplot(111, projection=extracted_data.wcs)
            plt.axis('off') #Removes axis
    #  Longitute ie RA
            lon = ax.coords[0]
            lon.set_major_formatter('hh:mm:ss.s')
            lon.set_ticks(color='black')
            lon.set_ticks_position('b')
            lon.set_ticklabel_position('b')
            lon.set_axislabel_position('t')
            lon.set_axislabel('RA')
    
    #  Latitude ie Dec
            lat = ax.coords[1]
            lat.set_major_formatter('dd:mm:ss')
            lat.set_ticks(color='black')
            lat.set_ticks_position('l')
            lat.set_ticklabel_position('l')
            lat.set_axislabel_position('r')
            lat.set_axislabel('Dec')
    
    # Histogram equalise the data for display        
            extracted_data_equalised = np.sort(np.ravel(extracted_data.data)).searchsorted(\
    extracted_data.data)
            fig.add_axes(ax)
            ax.imshow(extracted_data_equalised, cmap='gray',origin='lower')
    # Create jpg file name
            jpgname = f_name + '.jpg'
            print("Creating jpg ",jpgname)
    # Create jpg        
            plt.savefig(jpgname, bbox_inches='tight', pad_inches = 0)
            
            fitsname = f +'.fits'
            print("Creating fits file ",fitsname)
            hdu = fits.PrimaryHDU(extracted_data.data)
            hdu.writeto(fitsname,overwrite=True)
    
            plt.show()
    
    
        x = x + 1


