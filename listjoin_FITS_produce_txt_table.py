# I created this script in order to append several FITS tables. Right down here there are a couple of examples of whatever can be done+
# with FITS tables. In the end, I also print the table to a text file (which could even have a LaTex format!!!)

# coding: utf-8


import numpy as np
import scipy as scipy
import scipy.stats as stats
import astropy.stats as astats
import numpy.random as random
from astropy.table import Table
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table
import pyfits
import os
from os import walk
import glob
from tabulate import tabulate
import atpy



import operator

def sort_table(table, cols):
    # sort a table by multiple columns
    # table: a list of lists (or tuple of tuples) where each inner list represents a row
    # cols:  a list (or tuple) specifying the column numbers to sort by
    # e.g. (1,0) would sort by column 1, then by column 0
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table



# SELECT RECORDS FROM FITS TABLE
# In the next example, assuming the table’s second field having the name ‘magnitude’, 
# an output table containing all the records of magnitude > 5 from the input table is generated:

# t = pyfits.open('table.fits')
# tbdata = t[1].data
# mask = tbdata.['magnitude'] > 5
# newtbdata = tbdata[mask]
# hdu = pyfits.BinTableHDU(data=newtbdata)
# hdu.writeto('newtable.fits',clobber=True) # clobber=True allows to overwrite the existing FITS table



#t1 = pyfits.open('table1.fits')
#t2 = pyfits.open('table2.fits')



# MERGE FITS TABLES
# columns = t1[1].columns+t2[1].columns   
# hdu = pyfits.BinTableHDU.from_columns(columns)
# hdu.writeto('newtable.fits',clobber=True) # clobber=True allows to overwrite the existing FITS table



# APPEND  FITS TABLES   
# This cannot be done with cat *.fits > newtable.fits in the terminal!!!!!!
# ls *.fits > list.lst would create a file with all the names

#nrows1 = t1[1].data.shape[0]
#nrows2 = t2[1].data.shape[0]
#nrows = nrows1 + nrows2
#hdu = pyfits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows) 
#for colname in t1[1].columns.names:
#	hdu.data[colname][nrows1:] = t2[1].data[colname]
#hdu.writeto('newtable.fits',clobber=True) # clobber=True allows to overwrite the existing FITS table








# COMBINE ALL THE FITS FILES PRESENT IN THE CURRENT DIRECTORY


path=os.getcwd()
title="newtable"


# Find all the .fits files in the current directory


filelist=[]
for file in os.listdir(path):
	if file.endswith(".fits"):
		if file!=title+".fits":  # Allows to overwrite the final table and not to append it to the other ones
			filelist.append(file)
			filelist=sorted(filelist) # Sorts the files by name, given that Python reads them in a random way


# Open the tables


tables=[]
rows=[]
#columns=[]
for i in range(len(filelist)):
	tables.append(i)
	rows.append(i)
	#columns.append(i)
	tables[i] = pyfits.open(filelist[i])
	rows[i] = tables[i][1].data.shape[0]
	#columns[i]=tables[i][1].columns


nrows=sum(rows)
hdu = pyfits.BinTableHDU.from_columns(tables[0][1].columns, nrows=nrows) 
for colname in tables[0][1].columns.names:
	for i in range(len(tables)):
		hdu.data[colname][i] = tables[i][1].data[0][colname]
hdu.writeto(title+".fits",clobber=True) # clobber=True allows to overwrite the existing FITS table




# Print the final table to a .txt file and open it to write the header belonging to the previous tables and then overwrite it as a FITS again


tablefinal=pyfits.open(title+".fits")
tablefortext=Table.read(title+".fits",format='fits')


tableh=fits.open(filelist[0]) # In order to get the header for the final FITS table
prihdrr = tableh[0].header # I wrote the units as a comment when I created the individual tables. Header from one of them!




tbhdu = pyfits.BinTableHDU.from_columns(tablefinal[1].columns, nrows=nrows)
prihdr = pyfits.Header()
prihdr['COMMENT']=str(prihdrr['COMMENT'][0])+str(prihdrr['COMMENT'][1])   # It was divided into two. Weird. I am joining both elements
prihdu = pyfits.PrimaryHDU(header=prihdr)
tbhdulist = pyfits.HDUList([prihdu, tbhdu])
tbhdulist.writeto(title+".fits",clobber=True) #clobber=True allows to overwrite the existing FITS table

# prihdr['COMMENT'] was ivided into two after writing the table...




f=open(title+".txt","w")


print >>f, tabulate(tablefortext,headers=[str(i) for i in tables[0][1].columns.names])
# See the options of Pypi tabulate, the package which is being used. It allows to write a LaTex table, e.g. (tablefmt="latex")


print >>f,'\n'+'#'+str(prihdrr['COMMENT'][0])+str(prihdrr['COMMENT'][1]) # I include a comment with the units in the end of the table

f.close()



tabletext=Table.read(title+".txt",format='ascii')
tabletext_definitive=sort_table(tabletext,(0,1,2)) # Sort the table by column 0, then by column 1 and finally by column 2


g=open(title+".txt","w")

print >>g, tabulate(tabletext_definitive,headers=[str(i) for i in tables[0][1].columns.names,numalign="right"]) # Plain text

print >>g,'\n'+'#'+str(prihdrr['COMMENT'][0])+str(prihdrr['COMMENT'][1]) # I include a comment with the units in the end of the table

g.close()


h=open(title+"_latex.txt","w")

print >>h, tabulate(tabletext_definitive,headers=[str(i) for i in tables[0][1].columns.names],tablefmt="latex") # LaTex format

print >>h,'\n'+'#'+str(prihdrr['COMMENT'][0])+str(prihdrr['COMMENT'][1]) # I include a comment with the units in the end of the table

h.close()