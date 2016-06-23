# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:59:34 2016

@author: James
"""

import pandas as pd
import numpy as np

site = ['US-WCr', 'US-Whs', 'US-Wkg']
year_range = ['1999-2014', '2007-2014', '2004-2014']

#filename = '/air/projects/fluxnet2015/tier1/fullset/hh/FLX_US-'+site+'_FLUXNET2015_FULLSET_HH_'+year_range+'_1-1.csv'
filename = list()
for i in xrange(len(site)):
    z = 'FLX_'+site[i]+'_FLUXNET2015_FULLSET_HH_'+year_range[i]+'_1-1.csv'
    filename.append(z)

data = pd.DataFrame({'' : []})

for i in xrange(len(filename)):
    d = pd.read_csv(filename[i],na_values=['-9999'])
    data.append(d)

list(data.columns.values)

description = data.describe()
