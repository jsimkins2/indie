# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:34:12 2016

@author: James
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools
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
year = (2006,2007,2008,2010,2011,2012)
year = map(int, year)

neco2  = []
netime = []
nefc = []
newd = []
news = []
neta = []


# have to convert to float before you can run np.nan
for i in year:
    f = file('EFDC_L2_Flx_ATNeu_'+str(i)+'_v02_30m.txt', 'r')
    cols, name = getColumns(f, delim=",")
    #define each variable column
    c = cols['CO2']
    c = np.array(c)
    c = c.astype('float')
    c[c == -9999] = np.nan
    neco2.append(c)
    
    time = cols['DTime']
    time = np.array(time)
    time = time.astype('float')
    time[time == -9999] = np.nan
    netime.append(time)
    
    fc = cols['FC']
    fc = np.array(fc)
    fc = fc.astype('float')
    fc[fc == -9999] = np.nan
    nefc.append(fc)
    
    wd = cols['WD']
    wd = np.array(wd)
    wd = wd.astype('float')
    wd[wd == -9999] = np.nan
    newd.append(wd)
    
    ws = cols['WS']
    ws = np.array(ws)
    ws = ws.astype('float')
    ws[ws == -9999] = np.nan
    news.append(ws)
    
    ta = cols['TA']
    ta = np.array(ta)
    ta = ta.astype('float')
    ta[ta == -9999] = np.nan
    neta.append(ta)   
f.close()

chco2  = []
chtime = []
chfc = []
chwd = []
chws = []
chta = []


# have to convert to float before you can run np.nan
for i in year:
    f = file('EFDC_L2_Flx_CHCha_'+str(i)+'_v02_30m.txt', 'r')
    cols, name = getColumns(f, delim=",")
    #define each variable column
    c = cols['CO2']
    c = np.array(c)
    c = c.astype('float')
    c[c == -9999] = np.nan
    chco2.append(c)
    
    time = cols['DTime']
    time = np.array(time)
    time = time.astype('float')
    time[time == -9999] = np.nan
    chtime.append(time)
    
    fc = cols['FC']
    fc = np.array(fc)
    fc = fc.astype('float')
    fc[fc == -9999] = np.nan
    chfc.append(fc)
    
    wd = cols['WD']
    wd = np.array(wd)
    wd = wd.astype('float')
    wd[wd == -9999] = np.nan
    chwd.append(wd)
    
    ws = cols['WS']
    ws = np.array(ws)
    ws = ws.astype('float')
    ws[ws == -9999] = np.nan
    chws.append(ws)
    
    ta = cols['TA']
    ta = np.array(ta)
    ta = ta.astype('float')
    ta[ta == -9999] = np.nan
    chta.append(ta)   
f.close()
### Time for plotting ###

chwd = list(itertools.chain(*chwd))
chco2 = list(itertools.chain(*chco2))
neco2 = list(itertools.chain(*neco2))
    
# x1 - x2
#np.subtract(x1, x2)



    
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