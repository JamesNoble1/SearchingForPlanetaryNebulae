#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#"""
#Created on Thu Mar  4 17:12:37 2021
#
#@author: james
#"""


##########Imports########################
from astropy.io import fits
import numpy as np
import sys
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
import shutil
#########################################

#Command line usage 

#ReOrg4.py           Library_of_fits_files          Output_directory_of_sorted_fits     Output of subtracted_fits
#Use Full file paths

#########################################

#Purpose of script

#To match the telescope image file in to pairs.  These pairs will be of the same area of sky (same coordinates) but different wavebands
#Send these paired files to a subtraction program written by Dr Rhys Morris
#Output all files to a single directory.  These image files can then be stitched together (using AladinBeta.jar) and viewable on Aladin

#Outline

#Find information about the fits file - VPHAS object, filter type (red or H alpha), coordinates
#Rename file (for ease of use and sorting) as 'VPHASobject___filtertype___RAcoord___DECdoord'
#Match the file pairs by sorting by image separation
#Send each filter_Cordinate pair to subtraction program
#Compile resampled files to directory

#########################################

#the subtraction prgram - Resampler_quiet must be in same directory
 
#Will need to change SYS.ARGV arguments for unix

#####Sorts list by thrid column
def Sort(sub_li): 
    l = len(sub_li) 
    for i in range(0, l): 
        for j in range(0, l-i-1): 
            if (sub_li[j][2] > sub_li[j + 1][2]): 
                tempo = sub_li[j] 
                sub_li[j]= sub_li[j + 1] 
                sub_li[j + 1]= tempo 
    return sub_li 


# If we don't have 3 arguments, exit with warning message
if ( len(sys.argv) > 4 or len(sys.argv) < 4 ):
    print("There is a problem with your arguments.")
    sys.exit('Parameters are wrong')


os.chdir(sys.argv[1]) #Change to fits library


for f in os.listdir(): 
    filename = os.path.basename(f) #Gets filename
    f_name, f_ext = os.path.splitext(f)  #Splits filename and extension
    if f_ext == '.fz': #If there is a compressed fits file (fz) then exit
        print('This prgram does not accept fits.fz files. Please unpack them first!')
        sys.exit()
    if f_ext == '.fits': #If fits then continue 
        pass
    else:
        print('unknown file type')
        sys.exit()

    print('opening ', filename)
    hdu_list1 = fits.open(filename)
    
    Object_name = 'unknown'
    #Get object name from header
    try:
        Object_name = hdu_list1[0].header['OBJECT']
    except KeyError:
        print("Oops FITS header OBJECT NAME does not exist")
    #    filter_name = 'Unkown filter '
    
    # Try looking for a FILTER keyword
    try:
        Object_name = hdu_list1[0].header['Object']
    except KeyError:
        print("Oops FITS header Object does not exist")
    #    filter_name = 'Unknown filter'
    print('Object = ', Object_name)

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
    print('filter = ', filter_name)
    
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
        
    print("DEC of image centre is ",dec_whole_image)
    
    hdu_list1.close() 

    #Creates file name with vphas object, ra and dec with triple underscore seperation
    f_name = str((Object_name)) + '___' + (filter_name) + '___' + str(format(ra_whole_image, '.5f')) + '___' + str(format(dec_whole_image, '.5f'))
    new_name = '{}{}'.format(f_name, f_ext) #Creates new name with file extension
    print(new_name)
    if os.path.isfile(new_name) == True: ##If file already exists, delete the duplicate
        print(new_name, 'already exists, deleting the duplicate')
        os.remove(os.path.abspath(f))
    else:        
        os.rename(f, new_name) #If not, rename the file with the new name
    
N = len(os.listdir()) #Total number of files
i = 0 #Initiate variable
print('Total number of files = ',N)
#print(New_File_List)

#finds filter, ra and dec from the new file name
#Does this for one file
Possible_Files = [] #initialise possible array, this contans all possible file pairs
for f in os.listdir():  #iterate through directory
    filename1 = f
    f_name1, f_ext = os.path.splitext(f) #split filename and extension
    Object_name1 = (f_name1.split('___',3)[0]) #find vphas object
    filter1 = (f_name1.split('___',3)[1]) #find filter
    ra1 = float(f_name1.split('___',3)[2]) #find ra
    format(ra1, '.5f') #format ra to 5 sig fig
    dec1 = float(f_name1.split('___',3)[3]) #find dec
    format(dec1, '.5f') #Format dec to 5 sig fig
    c1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg)) #creaste skycoord
     
#Then does this for a different file
    for i in os.listdir():
        filename2 = i
        f_name2, f_ext = os.path.splitext(i) 
        Object_name2 = (f_name2.split('___',3)[0])
        filter2 = (f_name2.split('___',3)[1])
        ra2 = float(f_name2.split('___',3)[2])
        format(ra2, '.5f')
        dec2 = float(f_name2.split('___',3)[3])
        format(dec2, '.5f')
        c2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg))
        separation = c1.separation(c2) #find angular separation of the two image centers
#if filenames are different, separation less than 4000 arcsec, filter different, vphas object the same, then add to list of all possible matches
        if (f_name1 != f_name2) and (separation.arcsecond <= 4000) and (filter1 != filter2) and (Object_name1 == Object_name2): 
            print('matched', f_name1, 'and', f_name2)
            #append details to a list, separation is 3rd column
            NewLine = ([f_name1, f_name2, separation.arcsecond, ra1, dec1, ra2, dec2, Object_name1]) 
            Possible_Files.append(NewLine)
            
print(Possible_Files)
print('***************************Sorting in to ascending separation*****************************')
sub_li = Possible_Files
Sort(sub_li) #Sort list of possible matches by ascending separation, so closest files are matched first
print(Possible_Files)
print('********')

sub_li = np.array(sub_li) #make array
i = 0

#Send matched pairs to their own subdirectory

Parent_directory = sys.argv[2] #go to output of sorted fits
for i in range (len(sub_li)):     #iterate through list (closest first) sending matched pairs to a new subdirectory. 
        ra1 = format(float(sub_li[i,3]),'.5f') ######Formatting all to 5dp
        dec1 = format(float(sub_li[i,4]),'.5f')
        ra2 = format(float(sub_li[i,5]),'.5f')
        dec2 = format(float(sub_li[i,6]),'.5f')
        Object_name = sub_li[i,7]
                                    #Output directory specified by user input
        SubName = str(Object_name) + '___' + str(ra1) + '___' + str(dec1)+'___' +  str(ra2) + '___' + str(dec2) #Sub directory name 
        NewSub1 = os.path.join(Parent_directory, SubName)
        file1 ='{}{}'.format(sub_li[i,0], f_ext)
        file2 ='{}{}'.format(sub_li[i,1], f_ext)
        Path1 = os.path.abspath(file1) #old path
        Path2 = os.path.abspath(file2) #old path
        NewPath1 = os.path.join(NewSub1,file1)         #New Path for first file
        NewPath2 = os.path.join(NewSub1,file2)      #New Path for second file
       
        if os.path.isdir(NewSub1) == True: #If it is already a directiry, then skip file
            print(NewSub1, 'is already a directory')
            pass
        elif os.path.isfile(Path1) == False: ##If file already moved, skip
            #print(os.path.isfile(os.path.abspath(sub_li[i,0])))
            print(Path1, 'The file has already been moved')
            pass
        elif os.path.isfile(Path2) == False: ##If file already moved, skip
            print(Path2, 'The file has already been moved')
            pass
        elif os.path.isfile(NewPath1) == True: #If New File already exists, skip
            print(NewPath1, 'The new file already exists')
            pass
        elif os.path.isfile(NewPath2) == True: #If new file already exists, skip
            print(NewPath2, 'The new file already exists')
            pass
        elif os.path.isdir(NewSub1) == False: #If not, then make the sub directory and move the files to it
            os.mkdir(NewSub1)
            print('Moving',Path1 ,'and', Path2)
		
            shutil.move (Path1,NewPath1)  #First file moved to new sub directory
            shutil.move (Path2,NewPath2)    #Second File moved to new sub directory

            

os.chdir(sys.argv[2])

#Iterate through the sub directories sending the pairs to Rhys' subtraction program

Subtracted = sys.argv[3]
y = list(sys.argv)
i = 1
for d in os.listdir():
    print('**********************************************************')
    print('Completed: ', i,'/',len(os.listdir()))
    i = i + 1 #Progress Tally
    for f in os.listdir(d):
        #print('i', os.getcwd())
      
        Path1 = os.path.join(os.getcwd(),d)
        #print('Path1 = ', Path1)
 
        if (f.split('___',3)[1]) == 'NB_659':
            Path2 = os.path.join(Path1,f)
            
            y[1] = Path2
            sys.argv = tuple(y)
            #print('NB', y[1])
            continue
            
        
        if (f.split('___',3)[1]) == 'r_SDSS':
            Path2 = os.path.join(Path1,f)
            y[2] = Path2
            #print('r', y[2])
            sys.argv = tuple(y)
            continue
	
    d_ext = '___.fits'
    d_name = d + d_ext   
    Path3 = os.path.join(Subtracted,d_name) #this ensures correct arguments for resampler
    y[3] = Path3 #Output for subtracted files

    if os.path.isfile(Path3+'1'):
        continue
    else:
   		 exec(open('/data/cluster2/PN_project/JamesMax2021/Dataset_1/resampler_quiet.py').read())
        







