# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 16:11:42 2021

@author: james
"""

from astropy.io import fits
import numpy as np
import sys
import os
from astropy import SkyCoord
from astropy import u





#Usage 
#Hips_from_Fits.py Library_of_fits_files 
#Find if filter is r or Ha
#Find coordinates of images
#If cordinates are same(to some degree) and Filter is different, rename files as **filter_Coordinate**
#Send each filter_Cordinate pair to Resampler
#Compile resampled files to directory
#Send directory to AladinBeta.jar
 

sys.argv[1] = os.path('C:‪\\Users\\james\\Documents\\Fits_TEST')
All_Fits = os.listdir(sys.argv[1])


#N = Total_Number_Files

#Initialising file count
N = 0

#Counts number of files
for files in os.walk(sys.argv[1]):
    for files in sys.argv[1]:
        N = N+1
print('Total number of files = ', N)

#Initialise iteration value
i=0


#Loop to find filter type
for i in range (N):
    




print("opening ",Fits_1)
hdu_list1 = fits.open(Fits_1)

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






# Get coordinates of centre of image whole big image
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


#Loop for all files in directory starting from 0
i = 0
for i in range (N):
#Renames Filename dependent on filter
    if filter_name == 'r_sdss':
        file_name[i] = 'r_'
    
    elif filter_name == 'NB_659':
        file_name[i] = 'H_'
    else:
        file_name[i] = 'unknown_'
#Adds coordinates to file name
    file_name[i] = file_name[i] + c
    #nested loop
    x = 0
    for x in range (i):
        if 
    
i=0
#for i in range (N):
#If coordinate 
        



#N = Total_Number_Files

#Initialising file count and iteration
N = 0 #Total number of files
i = 0 #File number/iteration number

'''
#Counts number of files
for files in os.walk(sys.argv[1]):
    for files in sys.argv[1]:
        
        #os.rename(r'file path\OLD file name.file type',r'file path\NEW file name.file type')
        os.rename(r'sys.argv[1]\OLD file name.file type',r'file path\NEW file name.file type')
        
        N = N+1
        
'''        
        
print('Total number of files = ', N)