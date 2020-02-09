# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python

import math
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]
f3 = sys.argv[3]

out = open(f3,'w')

# length of slice to average over
L = 100000

# step
S = L

coverages = []
coordinates = []
cov = []

with open(f1) as handle:
    for line in handle:
        line = line.strip().split("\t")
 #       print line
        x = int(line[1])
        c = int(line[3])

        coordinates.append(x)
        coverages.append(c)
#print coverages
# normalize coverages?
m = sum(coverages)
#print m
for i in coverages:
    temp = int(i)/int(m)
    cov.append(temp)

#print cov

# divide by total ??


def slices(array, l=L, step=S):
    '''
    rolling slices
    '''

    for i in range(0, len(array), step):
    #    print i
        slice = array[i:i+l]
     #   print slice
        yield slice


for x, slice in zip(slices(coordinates),
                    slices(coverages)):

    # get the average coverage by slice
    y = sum(slice)/float(len(slice))

    # clip

#    print('chr8 {} {} {}'.format(
#        min(x),
#        max(x),
#        y
#    ))	
    temp =  str(f2) + "\t" + str(min(x)) + "\t" + str(max(x)) + "\t" + str(y)
    out.write(temp + "\n")
