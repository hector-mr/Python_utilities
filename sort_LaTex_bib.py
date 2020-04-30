#!/usr/bin/env python
# coding: utf-8

'''
Sort LaTex bibliography elements alphabetically based on first author, not on type (book, article, inproceedings,...)
'''

import argparse

def arg_parser():
    """
    Arguments to parse for model if run with console.
    :return args: dict of param names (str) and values (str)
    """
    parser = argparse.ArgumentParser(description='Argument parser')
    parser.add_argument('-f', '--file_name', metavar='str', type=str, help='Bibliography file name, including .bib', required=True)
    args = vars(parser.parse_args())
    return args


file_name = arg_parser()['file_name']


# Ancillary text parsing functions

def extract_data_lines(filename, start_text, end_text, include_start=False,
                       include_end=False):
    """
    open `filename`, and yield the lines between
    the line that contains `start_text` and the line that contains `end_text`.
    """
    started = False
    with open(filename) as fh:
        for line in fh:
            if started:
                if end_text in line:
                    if include_end:
                        yield line
                    break
                yield line
            elif start_text in line:
                started = True
                if include_start:
                    yield line


def find_data_lines(filename, text):
    """
    open `filename`, and yield the line that starts with `text`
    """
    with open(filename) as fh:

        li = []

        for line in fh:

            if line.startswith("@"):

                li.append(line.split(",")[0]) # Remove ",\n" from lines, only retrieve text
                
    return li


# FIND TEXT BETWEEN EACH ENTRY

start_entries = find_data_lines(file_name, '@')

# For the last entry, use the first one as end_text

end_entries = start_entries[1:] + [start_entries[0]]

# Empty list of empty lists

ads_groups = [ [] for i in range(len(start_entries)) ]

for s, st in enumerate(start_entries):
    
    for line in extract_data_lines(file_name, st, end_entries[s],
                                   include_start=True, include_end=False):
        
        ads_groups[s].append(line)


# SORT GROUPS

sorted_ads_groups = sorted(ads_groups, key=lambda x: x[0].split('{')[1])


# WRITE TO FINAL .BIB FILE

with open('{}_sorted.bib'.format(file_name.split('.')[0]), 'w') as text_file:
    
    for group in sorted_ads_groups:
        
        for line in group:
    
            text_file.write(line)
