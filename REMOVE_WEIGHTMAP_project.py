#!/usr/bin/python
"""
Created on Tue Mar  9 15:22:35 2021

@author: james
"""
import os
import sys

#####################################

#Arguments:

file1 = sys.argv [1] #Download Script
file2 = sys.argv [2] #Readme file from ESO website for requested download
Directory = sys.argv [3] #optional arg that contains directory of already downloaded files
###########Use      sed -i -e 's/\r$//' *******.sh      to convert .sh from windows to unix

################################################################################################################
#Original definitions from https://thispointer.com/python-how-to-delete-specific-lines-in-a-file-in-a-memory-efficient-way/

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

def delete_line_by_full_match(original_file, line_to_delete):
    """ In a file, delete the lines at line number in given list"""
    is_skipped = False
    dummy_file = original_file + '.bak'
    # Open original file in read only mode and dummy file in write mode
    with open(original_file, 'r') as read_obj:
        with open(dummy_file, 'w') as write_obj:
        # Line by line copy data from original file to dummy file
            for line in read_obj:
                line_to_match = line
                if line[-1] == '\n':
                    line_to_match = line[:-1]
                # if current line matches with the given line then skip that line
                if line_to_match != line_to_delete:
                    write_obj.write(line)
                else:
                    is_skipped = True
                    print('delete = ',line)
        # If any line is skipped then rename dummy file as original file
        if is_skipped:
            os.remove(original_file)
            os.rename(dummy_file, original_file)
        else:
            os.remove(dummy_file)
    
def check_string():
    with open(file1) as temp_f:
        datafile = temp_f.readlines()
    for line in datafile:
        if string in line:
            print('c', string)
            line = line.strip()
            ALL_LINES.append(line)
###################################################################################
        
f_name = [] #Initialise list

if len(sys.argv) == 4: #If optional argument included do this
    Parent = os.getcwd() #Get current directory
    file3 = open('FileList.txt', 'w+')    #Create file as writable
    os.chdir(Directory) #Change to Directory
    for f in os.listdir(Directory):
        print(f)
        file3.write(f+'\n') #Write files to txt
    os.chdir(Parent) #Change back to Parent directory

string = 'ANCILLARY.WEIGHTMAP' #Unwanted file type

with open(file2, "r") as fp: #Make a list of filenames to be deleted -from readme
    for line in lines_that_contain(string, fp):
        print (line[3:38])
        f_name.append(line[3:38])

if len(sys.argv) == 4: #Add all lines in file3 to a list
    file3 = open('FileList.txt')
    LINES = file3.readlines()
    for i in LINES:
        f_name.append(i)

ALL_LINES = [] #Initialise list   
print ('f_name =', f_name)

for i in f_name:
    string = i
    string = string.strip() #Removes \n  - Pytho new line character
    check_string()
  
with open(file1,'r') as fp:#Add Readme File to delete list
        for line in lines_that_contain('READ', fp):
            print(line)
            line = line.strip()
            ALL_LINES.append(line)

for i in range (len(ALL_LINES)):
    delete_line_by_full_match(file1, ALL_LINES[i]) #Delete lines

