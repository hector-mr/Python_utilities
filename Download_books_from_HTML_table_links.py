#!/usr/bin/env python
# coding: utf-8


import lxml.html as lh
import pandas as pd
import os
import requests
import shutil
import string
import subprocess
import urllib2

'''
Scrape HTML table with download links to free books. Download them using wget and rename using the table.
'''



# FOREWORD

# Define table url

url = 'https://docs.google.com/spreadsheets/d/1HzdumNltTj2SHmCv3SRdoub8SvpIEn75fa4Q23x0keU/htmlview#gid=793911758'

# Handle the contents of the website

page = requests.get(url)


# Store the contents of the website under doc

doc = lh.fromstring(page.content)


# Parse data that are stored between <tr>..</tr> of HTML.
# This was achieved by right-clicking on a row and hitting "inspect element"

tr_elements = doc.xpath('//tr')


# Check the length of the first 12 rows. In this case, it should be 5:
# Index BookTitle Author EnglishPackageName OpenURL

#print([len(T) for T in tr_elements[:12]])


# PARSE TABLE HEADER

# Create empty list

col = []

# For each row, store each first element (header) and an empty list

i = 0 


# Careful: on top of column headers, there are more: ["", "B", "L", "S"]
# Otherwise, select tr_elements[0]

for t in tr_elements[1]:
    i += 1
    name = t.text_content()
    #print('%d:"%s"'%(i, name))
    col.append((name,[]))


#print(len(col))


# CREATE PANDAS DATA FRAME

# Since the first two rows are the header, data is stored on the third row and onwards

for j in range(2, len(tr_elements)):
    
    # T is our j'th row

    T = tr_elements[j]
    
    # If row is not of len(col), the //tr data is not from our table

    if len(T) != len(col):

        break
    
    # i is the index of our column

    i = 0
    
    # Iterate through each element of the row

    for t in T.iterchildren():

        data = t.text_content()
        
        # Check if row is empty

        if i > 0:

        # Convert any numerical value to integers

            try:

                data = int(data)

            except:

                pass
            
        # Append the data to the empty list of the i'th column

        col[i][1].append(data)
        
        # Increment i for the next column

        i += 1


#print([len(C) for (title,C) in col])



Dict = {title:column for (title,column) in col}
df = pd.DataFrame(Dict)

# Drop original index column

df = df.drop("1", axis=1)

# Drop empty blank line at the beginning

df = df.drop(0, axis=0).reset_index(drop=True)


#print(df.head(5))
#print(df['English Package Name'].unique())


# MAIN LOOP TO DOWNLOAD BOOKS BY CATEGORIES


for c, ca in enumerate(df['English Package Name'].unique()):
    
    # Replace commas and blank spaces by underscores

    ctg = ca.replace(",", " ").replace(" ", "_").replace("__", "_")
    #print(ctg)
    
    
    # Create subdirectories for each category. In general, shell=True should NOT be used, but it must here

    book_dir = '/home/hector/Downloads/Springer_ebooks/'
    subprocess.Popen('mkdir {}'.format(book_dir), shell=True, 
                     stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    
    
    book_subd = '/home/hector/Downloads/Springer_ebooks/{}/'.format(ctg)
    p = subprocess.Popen('mkdir {}'.format(book_subd), shell=True, 
                     stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    #p.communicate()
    #p.wait()
    
    
    
    # Donwload command: wget
    # https://www.computerhope.com/unix/wget.htm
    # http://www.gnu.org/software/wget/manual/wget.html#Recursive-Accept_002fReject-Options
    
    # Explanation:
    # -nd: do not create directory structure, and rather download everything to the present directory
    # -r -l 0: download recursively, do not go beyond level 0
    # -A .pdf,epub: download pdf and epub files only
    # -R *_*.pdf,bbm*pdf,bfm*pdf: exclude chapter files (*_*.pdf) and bibliography files (bbm*pdf,bfm*pdf)
    # -e robots=off: execution without ancillary files
    
    
    # Download books and rename
        
    for link in df[df['English Package Name'] == ca]['OpenURL']:

        # URL gets redirected, so it is necessary to process it: https://stackoverflow.com/questions/4902523/how-to-get-the-url-of-a-redirect-with-python

        req = urllib2.Request(link)
        res = urllib2.urlopen(req)
        finalurl = res.geturl()

        # Download using wget
        # Careful: do NOT include blank spaces in any of these elements. Otherwise, subprocess.Popen will add double quotation marks

        wget_cmd = ["wget",
                    "-nd",
                    "-r",
                    "-l",
                    "0",
                    "-A",
                    ".pdf,.epub",
                    "-R",
                    "*_*.pdf,bbm*pdf,bfm*pdf",
                    "-e",
                    "robots=off",
                    "-P",
                    "'{}'".format(book_subd),
                    "'{}'".format(finalurl)]
        #print(wget_cmd)


        # Download files

        wg = subprocess.Popen(" ".join(wget_cmd), shell=True, 
            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        #print(subprocess.list2cmdline(wget_cmd))
        print(wg.communicate()[1])  # wg.args in Python3
        wg.wait()
        #print(link)


        # Original book names

        book_name = df[(df['English Package Name'] == ca) & (df['OpenURL'] == link)]['Book Title'].values[0]


        # Remove copyright, registered and trademark symbols
        # https://stackoverflow.com/questions/92438/stripping-non-printable-characters-from-a-string-in-python
        # https://www.utf8-chartable.de/unicode-utf8-table.pl?start=128&number=1024&utf8=string-literal&unicodeinhtml=hex

        book_name = filter(lambda x: x in string.printable, book_name)


        # Replace commas and blank spaces by underscores. Replace hyphens and dashes by underscores

        mapping = [ (",", " "), (" ", "_"), ("__", "_"), ("-", "_"), ("/", "_")]

        for k, v in mapping:

            book_name = book_name.replace(k, v)


        # Downloaded book names (part of the URL)

        downloaded_book_name = finalurl.split('/')[-1]

            
        # Rename files

        try:

            shutil.move('{}{}.pdf'.format(book_subd, downloaded_book_name), '{}{}.pdf'.format(book_subd, book_name))

        except IOError: # File does not exist

            pass

        except UnicodeEncodeError: # In case any characters beyond utf-8 remained in the book names. Should have been removed by string.printable, though

            print(book_name)
            pass

        try:

            shutil.move('{}{}.epub'.format(book_subd, downloaded_book_name), '{}{}.epub'.format(book_subd, book_name))

        except IOError:
            
            pass

        except UnicodeEncodeError:

            print(book_name)
            pass