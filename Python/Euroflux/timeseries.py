# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 21:11:14 2016

@author: James
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 19:15:22 2016

@author: James
"""

import pandas as pd
import matplotlib
from matplotlib import dates as d
import datetime as dt
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

year = (2006,2007,2008,2010,2011,2012)
year = map(int, year)

data = pd.read_excel('CH4_Neustiftfull.xlsx', 'CH4_Neustift10', skiprows = 0, index_col = 0)

#data['Time'] = data.index.map(lambda x: x.strftime("%H:%M"))
Wnspeed = []
data = data.loc[data['CH4_con'] > -90]
data = data.loc[data['WS'] >= 1]
Wndir = []
for row in data['WD']:
    if row >= 22.5 and row < 67.5:
        Wndir.append(1)
    elif row >= 67.5 and row < 122.5:
        Wndir.append(2)
    elif row >= 122.5 and row < 157.5:
        Wndir.append(3)
    elif row >= 157.5 and row < 202.5:
        Wndir.append(4)
    elif row >= 202.5 and row < 247.5:
        Wndir.append(5)
    elif row >= 247.5 and row < 292.5:
        Wndir.append(6)
    elif row >= 292.5 and row < 337.5:
        Wndir.append(7)
    elif row >= 337.5:
        Wndir.append(8)
    elif row < 22.5:
        Wndir.append(8)
    else:
        Wndir.append(9)

data['Wndir'] = Wndir


        
winds = ['NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']

corr1 = data.corr(method='pearson', min_periods=1)
#print corr1
#data.CH4_con[data.CH4_con == '#NV'] = np.nan


data = data.groupby('Wndir').describe().unstack()
#data = data.drop([9])

corr2 = data.corr(method='pearson', min_periods=1)


#data.index = pd.to_datetime(data.index.astype(str))

#corr2 = data.corr(method='pearson', min_periods=1)
#print corr2

#data.index = data.index.astype(list)


#corr = np.corrcoef(data['Wndir'], data['CH4_con'])


fig, ax = plt.subplots(1, figsize=(12,6))
ax.set_title('Chamau CO2 Flux vs. Wind Direction', fontsize=14, weight='bold')
ax.set_ylabel('CO2 (umol m-2 s-1)', fontsize=14, weight='bold')
ax.set_xlabel('Direction', fontsize=14, weight = 'bold')
#plt.xticks(data.index, winds)
#plt.ylim(0, 9000)
ax.set_xticklabels(winds)

#plt.scatter(data['WD'], data['CH4_con'])
ax.plot(data.index, data['CH4_con']['mean'], 'g', linewidth=2.0)
ax.plot(data.index, data['CH4_con']['75%'], color='g', linestyle = 'dashed')
ax.plot(data.index, data['CH4_con']['25%'], color='g', linestyle = 'dashed')
#ax.fill_between(data.index, data['CH4_con']['mean'], data['CH4_con']['75%'], alpha=.5, facecolor='g')
#ax.fill_between(data.index, data['CH4_con']['mean'], data['CH4_con']['25%'], alpha=.5, facecolor='g')

plt.tight_layout()
plt.show()

#dt.strptime(data.index, '%d.%m.%Y %H:%M').strftime('%H:%M')
