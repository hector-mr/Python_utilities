# I created this script in order to append several FITS tables. Right down here there are a couple of examples of whatever can be done
# with FITS tables. In the end, I also print the table to a text file (which could even have a LaTex format!!!)

# coding: utf-8


import numpy as np
import scipy as scipy
import scipy.stats as stats
import astropy.stats as astats
import numpy.random as random
from astropy.table import Table, vstack
import os
from os import walk
import glob
import operator




def sort_table(table, cols):
    # sort a table by multiple columns
    # table: a list of lists (or tuple of tuples) where each inner list represents a row
    # cols:  a list (or tuple) specifying the column numbers to sort by
    # e.g. (1,0) would sort by column 1, then by column 0
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table




# COMBINE ALL THE .txt FILES PRESENT IN THE CURRENT DIRECTORY


path=os.getcwd()
title="newtable"


# Find all the .txt files in the current directory


filelist=[]
for file in os.listdir(path):
	if file.endswith(".txt"):
		if file!=title+".txt":  # Allows to overwrite the final table and not to append it to the other ones
			filelist.append(file)
			filelist=sorted(filelist) # Sorts the files by name, given that Python reads them in a random way


# Open the tables


tables=[]
for i in range(len(filelist)):
	tables.append(i)
	tables[i] = Table.read(filelist[i],format='ascii.fixed_width')

t=Table(vstack([x for x in tables]))
t.write(title+".txt",format='ascii.fixed_width') 