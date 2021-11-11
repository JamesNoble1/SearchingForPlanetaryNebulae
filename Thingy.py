#!/usr/bin/python
"""
Created on Tue Mar  9 15:22:35 2021

@author: james
"""
import os
import sys
#file1 = 'C:/Users/james/Documents/Python_TEST/Text_TEST/downloadRequest600807script.sh'
#file2 = 'C:/Users/james/Documents/Python_TEST/Text_TEST/README_600807.txt'
#DELETE = 0
#################Input 1 is .sh, input 2 is README
file1 = sys.argv [1]
file2 = sys.argv [2]

string = 'ANCILLARY.WEIGHTMAP'

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

def delete_line_by_full_match(original_file, line_to_delete):
    """ In a file, delete the lines at line number in given list"""
    is_skipped = False
    dummy_file = original_file + '.bak'
    # Open original file in read only mode and dummy file in write mode
    with open(original_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
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
    # If any line is skipped then rename dummy file as original file
    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)

f_name = []
with open(file2, "r") as fp:
    for line in lines_that_contain(string, fp):
        print (line[3:38])
        f_name.append(line[3:38])
    
ALL_LINES = []    
print (f_name)

for i in range (len(f_name)):
    with open(file1) as fp:
        for line in lines_that_contain(f_name[i], fp):
            line = line.strip()
            ALL_LINES.append(line)
         

        
with open(file1) as fp:
        for line in lines_that_contain('READ', fp):
            print(line)
            line = line.strip()
            ALL_LINES.append(line)

print(ALL_LINES)

for i in range (len(ALL_LINES)):
    delete_line_by_full_match(file1, ALL_LINES[i])
