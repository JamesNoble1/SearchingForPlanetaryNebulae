#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:06:10 2021

@author: james
"""
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os


sys.argv = [0,'G:/Subtracted_Sample/FITS','G:/Subtracted_Sample/JPEG']

os.chdir(sys.argv[1])
for f in os.listdir(sys.argv[1]):
    #print('checking file', f)
    image = os.path.abspath(f)
    print("opening ",f)
    hdu_list1 = fits.open(image)
    #print(hdu_list1)
    print(hdu_list1.info())
    image_data = fits.getdata(image, ext=0)
    data_equalised = np.sort(np.ravel(image_data.data)).searchsorted(image_data.data)
    #fig=plt.figure(figsize=(7,7) )
    plt.axis('off') #Removes axis
    #plt.imshow(data_equalised, cmap='gray',origin='lower') #Shows image
    #
    #Create JPG Path
    jpgname = os.path.join(sys.argv[2], f) + '.jpg' 
    
    #plt.savefig(jpgname, bbox_inches='tight', pad_inches = 0) #Saves fig and removes white border
    #plt.imsave(fname=jpgname, arr=data_equalised, cmap='gray', format='jpg')
