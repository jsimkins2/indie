# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 12:01:53 2016

@author: James
"""

import numpy as np
import numpy.fft as FFT
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series, DataFrame, Panel

site = ['US-WCr', 'US-Whs', 'US-Wkg']
year_range = ['1999-2014', '2007-2014', '2004-2014']

filename = list()
for i in xrange(len(site)):
    z = 'FLX_'+site[i]+'_FLUXNET2015_SUBSET_HH_'+year_range[i]+'_1-1.csv'
    filename.append(z)


data = pd.read_csv(filename[0],na_values=['-9999'])

var = data['NEE_VUT_REF']
dates = pd.date_range('1999-01-01 00:00', periods=280512, freq='30min').time

s = pd.Series(data['NEE_VUT_REF'].values, index = dates)
