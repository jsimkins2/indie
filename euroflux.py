# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:34:12 2016

@author: James
"""
import numpy as np
import matplotlib.pyplot as plt
def getColumns(inFile, delim="\t", header=True):
    """
    Get columns of data from inFile. The order of the rows is respected
    
    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that 
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in 
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the 
    heading names (i.e. the keys in the cols dict)
    """
    cols = {}
    name = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    name[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    name[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[name[i]] += [cell]
                i += 1
                
    return cols, name
year = (2006,2008,2010,2011,2012)
year = map(int, year)
chco2  = []


# have to convert to float before you can run np.nan
for i in year:
    f = file('EFDC_L2_Flx_CHCha_'+str(i)+'_v02_30m.txt', 'r')
    cols, name = getColumns(f, delim=",")
    c = cols['CO2']
    c = np.array(c)
    c = c.astype('float')
    c[c == -9999] = np.nan
    chco2.append(c)
    
'''
for n,i in enumerate(chco2[0]):
    if i==-9999:
      chco2[n]=np.nan

for i in chco2[]:
[4 if x==1 else x for x in a]
    chco2= np.array(chco2, dtype = float)
    #chco2 = chco2.astype('float')
    chco2[chco2 == -9999] = np.nan
-
co2 = cols['CO2']
co2 = np.array(chco2[0])
co2 = co2.astype('float')
co2[co2 == -9999] = np.nan
'''