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

####################################################################
#Outline
#Cycles through subtracted images and checks for PN coordinates within the image
#Creates a cutout image if any are found
#This is to build an image database for the machine learning model

#Usage
#arguments:
#vst_cutout_all.py, subtracted_fits_files_directory, csv_list_of_PN, Image_output_directory, Image_Size
###################################################################
#cdsys.argv = [0,'G:/Subtracted_Sample','C:/Users/james/OneDrive/Desktop/Project/TruePN.csv', 'G:/Images',500 ] #==Test parameters
All_PN = []

with open(sys.argv[2]) as csvfile:
    reader = csv.reader(csvfile) # Reads CSV and puts content in a list
    for row in reader: 
        All_PN.append(row)

All_PN = np.delete(All_PN, (0), axis=0) #Removes the title of columns
All_PN_formatted = []
for i in All_PN:
        #Columns in CSV where iinfo is stored    
        ra = i[3]  ###############################Note####################################
        dec = i[4]#################################These 2 parameters depend on your CSV file
        
        ra_hours = ra.split(':',3)[0]
        ra_mins = ra.split(':',3)[1]
        ra_secs = ra.split(':',3)[2]
        ra_str = str(ra_hours + " " + ra_mins + " " + ra_secs)
        #print(ra_str)
        dec_degrees = dec.split(':',3)[0]
        dec_mins = dec.split(':',3)[1]
        dec_secs = dec.split(':',3)[2]
        dec_str = str(dec_degrees + " " + dec_mins + " " + dec_secs)
        #print(dec_str)
        All_PN_formatted.append([ra_str,dec_str])



#All_PN.append([0,0,0,'17:50:40','-24:58:46.49']) #=== tester - these coords should be in first image in sample folder
os.chdir(sys.argv[1]) #Change to subtracted subdirectory and cycle through

path1 = os.path.join(sys.argv[3], 'JPEG') #Make jpeg and fits directories
if os.path.isdir(path1) == False: 
    os.mkdir(path1)
else:
    pass

path2 = os.path.join(sys.argv[3], 'FITS')
if os.path.isdir(path2) == False:
    os.mkdir(path2)
else:
    pass

N = 0

for f in os.listdir(sys.argv[1]):
    #print('checking file', f)
    image1 = os.path.abspath(f)
    print("opening ",f)
    hdu_list1 = fits.open(image1)
    print(hdu_list1)
    print('Checking if known PN in', image1)
            
    # How many subimages do we have?

    try:
        nimages = hdu_list1[0].header['HIERARCH ESO DET CHIPS']
    except KeyError:
        print("Oops FITS header HIERARCH ESO DET CHIPS does not exist, setting nimages to 1")
        nimages = 1
        
    print("number of images is ",nimages,"\n")
    for i in All_PN_formatted: #Iterate through list

        # Now create skycoord class
        obj_coord = SkyCoord(i[0], i[1], frame='icrs',unit=(u.hourangle,u.deg) )

        isinimage = 'False'

        x = 0

        wcs_info = wcs.WCS(hdu_list1[x].header)
        
        isinimage = obj_coord.contained_by(wcs_info)
        
        if isinimage == True:
            print("These coordinates are in subimage ",x)
        #        print("The type of side is ",type(sys.argv[8]) )
        #  Extract a cutout box
            side = int(sys.argv[4])
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
            extracted_data_equalised = np.sort(np.ravel(extracted_data.data)).searchsorted(extracted_data.data)
            fig.add_axes(ax)
            ax.imshow(extracted_data_equalised, cmap='gray',origin='lower')
            # Create jpg file name
            jpgname = 'cutout'+'_'+i[0]+'_'+i[1]+'.jpg'
            jpgname = jpgname.replace(" ","")
            jpgname = os.path.join(path1,jpgname)
            print("Creating jpg ",jpgname)
            # Create jpg        
            plt.savefig(jpgname)
            # Create fits file name
            fitsname = 'cutout'+'_'+i[0]+'_'+i[1]+'.fits'
            fitsname = fitsname.replace(" ","")
            fitsname = os.path.join(path2,fitsname)
            print("Creating fits file ",fitsname)
            hdu = fits.PrimaryHDU(extracted_data.data)
            hdu.writeto(fitsname,overwrite=True)
            
            
            #plt.show()
            
            N = N + 1
print('Total cutouts = ', N)
        


