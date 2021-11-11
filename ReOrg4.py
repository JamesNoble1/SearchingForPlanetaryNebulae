#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#"""
#Created on Thu Mar  4 17:12:37 2021
#
#@author: james
#"""

from astropy.io import fits
import numpy as np
import sys
import os
from os import path
from astropy import units as u
from astropy.coordinates import SkyCoord
import shutil





#Usage 
#Reorg.py           Library_of_fits_files          Output_directory_of_sorted_fits     Output of subtracted_fits
#Use Full file paths

#Find if filter is r or Ha
#Find coordinates of images
#If cordinates are same(to some degree) and Filter is different, rename files as **filter_Coordinate**
#Send each filter_Cordinate pair to Resampler
#Compile resampled files to directory
#Send directory to AladinBeta.jar - probably unecessary

#########Resampler_quiet must be in same directory
 
######Will need to change SYS.ARGV stuff for unix



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
    print("Argument1 should be a library of fits files, argument 2 should be a directory for the sorted fits files and argument 4 should be the output folder")

    sys.exit('Parameters are wrong')


# Python program to rename all file 
# names in your directory  

  
os.chdir(sys.argv[1]) 


for f in os.listdir(): 
    filename = os.path.basename(f)
    f_name, f_ext = os.path.splitext(f) 
    #print(filename)
    #print(f_name, f_ext)
    if f_ext == '.fz':
        print('This prgram does not accept fits.fz files. Please unpack them first!')
        sys.exit()
    if f_ext == '.fits':
        pass
    else:
        print('unknown file type')
        sys.exit()
    

#########Rhys' Code############## Altered version

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
    
    #c = SkyCoord(ra_whole_image, dec_whole_image, unit=(u.hourangle, u.deg))
    #hmsdms_coords = c.to_string('hmsdms')
    #print ('c is ', c, 'hmsdms is ', hmsdms_coords)
    
##############################################################################
#James' Code    
    #Renames file with ra and dec after triple underscore
    f_name = str((Object_name)) + '___' + (filter_name) + '___' + str(format(ra_whole_image, '.5f')) + '___' + str(format(dec_whole_image, '.5f'))
    #increment() 
  
    new_name = '{}{}'.format(f_name, f_ext) 
    print(new_name)
    if os.path.isfile(new_name) == True: ##If file already moved, skip
            #print(os.path.isfile(os.path.abspath(sub_li[i,0])))
        print(new_name, 'already exists, deleting the duplicate')
        os.remove(os.path.abspath(f))
    else:        
        os.rename(f, new_name)
    
    ###NOTE filter name is characters 1-6

#New_File_List = []
#New_File_List = os.listdir()
N = len(os.listdir())
i = 0
print('Total number of files = ',N)
#print(New_File_List)

#If separation is less than 0.0111111 degrees = 40 arcseconds = Rhys' Resampler_quiet.py parameter   
#If this part of the file name is different, and this other part of the file name is the same as another file then move them into their
#own list and then send that list to resampler_quiet.py



#finds filter, ra and dec from the new file name

#Does this for one file
Possible_Files = [] #initialise possible array
for f in os.listdir():
    filename1 = f
    f_name1, f_ext = os.path.splitext(f) 
    Object_name1 = (f_name1.split('___',3)[0])
    filter1 = (f_name1.split('___',3)[1])
    ra1 = float(f_name1.split('___',3)[2])
    format(ra1, '.5f')
    dec1 = float(f_name1.split('___',3)[3])
    format(dec1, '.5f')
    c1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
    
    
    
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
        #separation = ((ra1-ra2)**2 + (dec1-dec2)**2)**(1/2)
        c2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg))
        separation = c1.separation(c2)
#Compares files names #speration was 40
        if (f_name1 != f_name2) and (separation.arcsecond <= 4000) and (filter1 != filter2) and (Object_name1 == Object_name2): #If file 2 is nearer than 40 arcseconds and different filter, add to list
            print('matched', f_name1, 'and', f_name2)
            NewLine = ([f_name1, f_name2, separation.arcsecond, ra1, dec1, ra2, dec2, Object_name1])
           
            Possible_Files.append(NewLine)
            
print(Possible_Files)
print('***************************Sorting in to decreasing separation*****************************')

#print(Possible_Files)
sub_li = Possible_Files
Sort(sub_li) #Sorts files by ascending separation *********SORT BY DATE?
print(Possible_Files)
print('********')

sub_li = np.array(sub_li)


i = 0

Parent_directory = sys.argv[2]
for i in range (len(sub_li)):     
        ra1 = format(float(sub_li[i,3]),'.5f') ######Formatting all to 5dp
        dec1 = format(float(sub_li[i,4]),'.5f')
        ra2 = format(float(sub_li[i,5]),'.5f')
        dec2 = format(float(sub_li[i,6]),'.5f')
        Object_name = sub_li[i,7]
                                    #Output directory specified by user input
        SubName = str(Object_name) + '___' + str(ra1) + '___' + str(dec1)+'___' +  str(ra2) + '___' + str(dec2) #Sub directory name 
        #print(SubName)
       
        NewSub1 = os.path.join(Parent_directory, SubName)
        #print(NewSub1)
        file1 ='{}{}'.format(sub_li[i,0], f_ext)
        file2 ='{}{}'.format(sub_li[i,1], f_ext)
        Path1 = os.path.abspath(file1)
        Path2 = os.path.abspath(file2)
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
            #print(Path1)
            #print(NewPath1)
            #print(Path2)
            #print(NewPath2)
            #print('new1 = ', NewPath1 , 'old = ', Path1)
            print('Moving',Path1 ,'and', Path2)
		
            shutil.move (Path1,NewPath1)  #First file moved to new sub directory
            shutil.move (Path2,NewPath2)    #Second File moved to new sub directory
  
        

            



            
#Change to the new directory
os.chdir(sys.argv[2])

#sys.arg[3] = output of subtracted fits
#NB file will always be at the top, and r will always be after it so this will work, but there is probably a more elegant solution
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
    #print('d2 = ', d)
    #print(d)
    d_ext = '___.fits'
    d_name = d + d_ext   
    Path3 = os.path.join(Subtracted,d_name)
    y[3] = Path3

    if os.path.isfile(Path3+'1'):
        continue
    else:
    #sys.argv = tuple(y)
    #print(sys.argv[1])
   		 exec(open('/data/cluster2/PN_project/JamesMax2021/Dataset_1/resampler_quiet.py').read())
        







