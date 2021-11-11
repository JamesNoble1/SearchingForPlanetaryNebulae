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
import csv
import os



#Subdirectory of subbed fits = sys.argv[1]    
#PN_File = sys.arg[2]  ==== CSV file of knwown PN  
#box size = sys.argv[4]
#sys.argv = [0,'G:/Subtracted_Sample','C:/Users/james/OneDrive/Desktop/Project/TruePN.csv', 'G:/Images',50 ] #==Test parameters

#sys.argv[1] = 'G:/Subtracted_Sample' ====input folder of single fits files
#sys.argv[3] = output directory for images

All_PN = []

with open(sys.argv[2]) as csvfile:
    reader = csv.reader(csvfile) # Reads CSV and puts content in a list
    for row in reader: 
        All_PN.append(row)

All_PN = np.delete(All_PN, (0), axis=0) #Removes the title of columns


#All_PN.append([0,0,0,'17:50:40','-24:58:46.49']) #=== tester - these coords should be in first image in sample folder
os.chdir(sys.argv[1]) #Change to subtracted subdirectory and cycle through
for f in os.listdir(sys.argv[1]):
    #print('checking file', f)
    image1 = os.path.abspath(f)
    
    
    print("opening ",f)
    hdu_list1 = fits.open(image1)
    print(hdu_list1)
    print('Checking if known PN in', image1)
    for i in All_PN: #Iterate through list
        
        ra = i[3]
        dec = i[4]
        
        
        

        # Use contained_by(wcs[,image[) from
        # https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
        
        
        # Assign variables to input parameters
        ra_hours = ra.split(':',3)[0]
        ra_mins = ra.split(':',3)[1]
        ra_secs = ra.split(':',3)[2]
        ra_str = ra_hours + " " + ra_mins + " " + ra_secs
        
        #print("RA string is ",ra_str)
        
        dec_degrees = ra.split(':',3)[0]
        dec_mins = ra.split(':',3)[1]
        dec_secs = ra.split(':',3)[2]
        dec_str = dec_degrees + " " + dec_mins + " " + dec_secs
        
        #print("Dec string is ",dec_str)
        
        # Now create skycoord class
        obj_coord = SkyCoord(ra_str, dec_str, frame='icrs',unit=(u.hourangle,u.deg) )
        #print("RA in degrees is ",obj_coord.ra.degree)
        #print("Dec in degrees is ",obj_coord.dec.degree)
        
        # Get header data for subimages
        header_list = fits.open(image1)
        
        
        
        isinimage = 'False'
        #print(ra_str,dec_str)
        # Loop over subimages
        #x = 1
        #while x <= nimages:
        #    print("Extension is ",x)
        #    image_data = fits.getdata(sys.argv[1], ext=x)
        wcs_info = wcs.WCS(header_list[0].header)
        print (wcs_info)
        #print(wcs_info)
        #print()
        isinimage = obj_coord.contained_by(wcs_info)
        
        if isinimage == True:
            print("These coordinates",ra,dec," are in image ",f)
    #        print("The type of side is ",type(sys.argv[8]) )
    #  Extract a cutout box
            side = int(sys.argv[3])
            boxsize = [side, side ]
    #        print("box size is ",boxsize)
    #  Get the data for the subimage containing our coordinates.
            image_data = fits.getdata(image1)
    #  Extract the data        
            extracted_data = Cutout2D(image_data, obj_coord, boxsize, wcs=wcs_info )
    #  Create a popup figure window
            fig=plt.figure(figsize=(6,6) )
    #  set up axes using wcs (world coordinate system) of extracted data        
            ax = fig.add_subplot(111, projection=extracted_data.wcs)
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
            ax.imshow(extracted_data_equalised, cmap='hot',origin='lower')
    # Create jpg file name
            jpgname = 'cutout'+ra_str+dec_str+'.jpg'
            jpgname = jpgname.replace(" ","")
            jpgname = os.path.join(sys.argv[3], jpgname)
            print("Creating jpg ",jpgname)
    # Create jpg        
            plt.savefig(jpgname)
    # Create fits file name
            fitsname = 'cutout'+'_'+ra_str+'_'+dec_str+'.fits'
            fitsname = fitsname.replace(" ","")
            print("Creating fits file ",fitsname)
            hdu = fits.PrimaryHDU(extracted_data.data)
            hdu.writeto(fitsname,overwrite=True)
    
    
            plt.show()
    
        else:
            pass
        #x = x + 1
    

