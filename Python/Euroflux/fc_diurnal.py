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

data = pd.read_excel('Neustift.xlsx', 'Neu', skiprows = 0, index_col = 0)

data['Time'] = data.index.map(lambda x: x.strftime("%H:%M"))
#data.CH4_con[data.CH4_con == '#NV'] = np.nan

data = data.groupby('WD').describe().unstack()


data.index = pd.to_datetime(data.index.astype(str))

#data['means']=data[['CH4_con']].mean(axis=1)
fig, ax = plt.subplots(1, figsize=(12,6))
ax.set_title('CO2 Flux Diurnal Profile', fontsize=14, weight='bold')
ax.set_ylabel('CO2 Flux (umol m-2 s-1)', fontsize=14, weight='bold')
ax.set_xlabel('Time of Day', fontsize=14, weight = 'bold')

ax.plot(data.index, data['FC']['mean'], 'g', linewidth=2.0)
ax.plot(data.index, data['FC']['75%'], color='g')
ax.plot(data.index, data['FC']['25%'], color='g')
ax.fill_between(data.index, data['FC']['mean'], data['FC']['75%'], alpha=.5, facecolor='g')
ax.fill_between(data.index, data['FC']['mean'], data['FC']['25%'], alpha=.5, facecolor='g')

plt.tight_layout()
plt.show()

#dt.strptime(data.index, '%d.%m.%Y %H:%M').strftime('%H:%M')
