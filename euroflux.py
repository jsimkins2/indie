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



f = file('CH4_Neustift10.csv', 'r')
cols, name = getColumns(f, delim=",")

methc = []

me = cols['MeanCH4_[ppb]']
me = np.array(me)
me[me == '#NV'] = np.nan
me = me.astype('float')
methc.append(me)
f.close()

year1 = (2010,2011,2012)
year1 = map(int, year1)

methco2  = []
methtime = []
methfc = []
methwd = []
methws = []
methta = []


# have to convert to float before you can run np.nan
for i in year1:
    f = file('EFDC_L2_Flx_ATNeu_'+str(i)+'_v02_30m.txt', 'r')
    cols, name = getColumns(f, delim=",")
    #define each variable column
    c = cols['CO2']
    c = np.array(c)
    c = c.astype('float')
    c[c == -9999] = np.nan
    methco2.append(c)
    
    time = cols['DTime']
    time = np.array(time)
    time = time.astype('float')
    time[time == -9999] = np.nan
    methtime.append(time)
    
    fc = cols['FC']
    fc = np.array(fc)
    fc = fc.astype('float')
    fc[fc == -9999] = np.nan
    methfc.append(fc)
    
    wd = cols['WD']
    wd = np.array(wd)
    wd = wd.astype('float')
    wd[wd == -9999] = np.nan
    methwd.append(wd)
    
    ws = cols['WS']
    ws = np.array(ws)
    ws = ws.astype('float')
    ws[ws == -9999] = np.nan
    methws.append(ws)
    
    ta = cols['TA']
    ta = np.array(ta)
    ta = ta.astype('float')
    ta[ta == -9999] = np.nan
    methta.append(ta)   
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
### appending the individual year data to one long continuous set
# Chamau data
chwd_com = list(itertools.chain(*chwd))
chco2_com = list(itertools.chain(*chco2))
chfc_com = list(itertools.chain(*chfc))
chws_com = list(itertools.chain(*chws))
newd_com = list(itertools.chain(*newd))
neco2_com = list(itertools.chain(*neco2))
nefc_com = list(itertools.chain(*nefc))
news_com = list(itertools.chain(*news))
methwd_com = list(itertools.chain(*methwd))
methco2_com = list(itertools.chain(*methco2))
methfc_com = list(itertools.chain(*methfc))
methws_com = list(itertools.chain(*methws))
methtime_com = list(itertools.chain(*methtime))


chwd_com = np.array(chwd_com)
chco2_com = np.array(chco2_com)
chfc_com = np.array(chfc_com)
chws_com = np.array(chws_com)
newd_com = np.array(newd_com)
neco2_com = np.array(newd_com)
nefc_com = np.array(nefc_com)
news_com = np.array(news_com)
methwd_com = np.array(methwd_com)
methco2_com = np.array(methwd_com)
methfc_com = np.array(methfc_com)
methws_com = np.array(methws_com)
methc_com = np.array(methc)
methtime_com = np.array(methtime_com)
print np.nanmean(nefc_com)
print np.nanstd(nefc_com)

dtime = range(0, 105217)
dtime = np.array(dtime)

cham = np.row_stack((chwd_com, chws_com,chfc_com, dtime))
neu = np.row_stack((newd_com,news_com,nefc_com, dtime))
wolf = np.row_stack((methwd_com,methws_com,methfc_com,methc_com))

#northerly
ch_northerly = np.where((cham[0] > 337.5) & (cham[0] < 22.5) & (cham[1] > 4))
ne_northerly = np.where((neu[0] > 337.5) & (neu[0] < 22.5) & (neu[1] > 4))
wolf_northerly = np.where((wolf[0] > 337.5) & (wolf[0] < 22.5) & (wolf[1] > 4))
#northeasterly
ch_northeasterly = np.where((cham[0] > 22.5) & (cham[0] < 67.5) & (cham[1] > 4))
ne_northeasterly = np.where((neu[0] > 22.5) & (neu[0] < 67.5) & (neu[1] > 4))
wolf_northeasterly = np.where((wolf[0] > 22.5) & (wolf[0] < 67.5) & (wolf[1] > 4))
#easterly
ch_easterly = np.where((cham[0] > 67.5) & (cham[0] < 112.5) & (cham[1] > 4))
ne_easterly = np.where((neu[0] > 67.5) & (neu[0] < 112.5) & (neu[1] > 4))
wolf_easterly = np.where((wolf[0] > 67.5) & (wolf[0] < 112.5) & (wolf[1] > 4))
#Southeasterly
ch_southeasterly = np.where((cham[0] > 112.5) & (cham[0] < 157.5) & (cham[1] > 4))
ne_southeasterly = np.where((neu[0] > 112.5) & (neu[0] < 157.5) & (neu[1] > 4))
wolf_southeasterly = np.where((wolf[0] > 112.5) & (wolf[0] < 157.5) & (wolf[1] > 4))
#southerly
ch_southerly = np.where((cham[0] > 157.5) & (cham[0] < 202.5) & (cham[1] > 4))
ne_southerly = np.where((neu[0] > 157.5) & (neu[0] < 202.5) & (neu[1] > 4))
wolf_southerly = np.where((wolf[0] > 157.5) & (wolf[0] < 202.5) & (wolf[1] > 4))
#southwesterly
ch_southwesterly = np.where((cham[0] > 202.5) & (cham[0] < 247.5) & (cham[1] > 4))
ne_southwesterly = np.where((neu[0] > 202.5) & (neu[0] < 247.5) & (neu[1] > 4))
wolf_southwesterly = np.where((wolf[0] > 202.5) & (wolf[0] < 247.5) & (wolf[1] > 4))
#westerly
ch_westerly = np.where((cham[0] > 247.5) & (cham[0] < 292.5) & (cham[1] > 4))
ne_westerly = np.where((neu[0] > 247.5) & (neu[0] < 292.5) & (neu[1] > 4))
wolf_westerly = np.where((wolf[0] > 247.5) & (wolf[0] < 292.5) & (wolf[1] > 4))
#northwesterly
ch_northwesterly = np.where((cham[0] > 292.5) & (cham[0] < 337.5) & (cham[1] > 4))
ne_northwesterly = np.where((neu[0] > 292.5) & (neu[0] < 337.5) & (neu[1] > 4))
wolf_northwesterly = np.where((wolf[0] > 292.5) & (wolf[0] < 337.5) & (wolf[1] > 4))
###############################################################################
chfc_n = np.take(cham[2], ch_northerly)
chwd_n = np.take(cham[0], ch_northerly)
nefc_n = np.take(neu[2], ne_northerly)
newd_n = np.take(neu[0], ne_northerly)
methfc_n = np.take(wolf[2], wolf_northerly)
methwd_n = np.take(wolf[0], wolf_northerly)
methc_n = np.take(wolf[3], wolf_northerly)

chfc_ne = np.take(cham[2], ch_northeasterly)
chwd_ne = np.take(cham[0], ch_northeasterly)
nefc_ne = np.take(neu[2], ne_northeasterly)
newd_ne = np.take(neu[0], ne_northeasterly)
methfc_ne = np.take(wolf[2], wolf_northeasterly)
methwd_ne = np.take(wolf[0], wolf_northeasterly)
methc_ne = np.take(wolf[3], wolf_northeasterly)

chfc_e = np.take(cham[2], ch_easterly)
chwd_e = np.take(cham[0], ch_easterly)
nefc_e = np.take(neu[2], ne_easterly)
newd_e = np.take(neu[0], ne_easterly)
methfc_e = np.take(wolf[2], wolf_easterly)
methwd_e = np.take(wolf[0], wolf_easterly)
methc_e = np.take(wolf[3], wolf_easterly)

chfc_se = np.take(cham[2], ch_southeasterly)
chwd_se = np.take(cham[0], ch_southeasterly)
nefc_se = np.take(neu[2], ne_southeasterly)
newd_se = np.take(neu[0], ne_southeasterly)
methfc_se = np.take(wolf[2], wolf_southeasterly)
methwd_se = np.take(wolf[0], wolf_southeasterly)
methc_se = np.take(wolf[3], wolf_southeasterly)

chfc_s = np.take(cham[2], ch_southerly)
chwd_s = np.take(cham[0], ch_southerly)
nefc_s = np.take(neu[2], ne_southerly)
newd_s = np.take(neu[0], ne_southerly)
methfc_s = np.take(wolf[2], wolf_southerly)
methwd_s = np.take(wolf[0], wolf_southerly)
methc_s = np.take(wolf[3], wolf_southerly)

chfc_sw = np.take(cham[2], ch_southwesterly)
chwd_sw = np.take(cham[0], ch_southwesterly)
nefc_sw = np.take(neu[2], ne_southwesterly)
newd_sw = np.take(neu[0], ne_southwesterly)
methfc_sw = np.take(wolf[2], wolf_southwesterly)
methwd_sw = np.take(wolf[0], wolf_southwesterly)
methc_sw = np.take(wolf[3], wolf_southwesterly)

chfc_w = np.take(cham[2], ch_westerly)
chwd_w = np.take(cham[0], ch_westerly)
nefc_w = np.take(neu[2], ne_westerly)
newd_w = np.take(neu[0], ne_westerly)
methfc_w = np.take(wolf[2], wolf_westerly)
methwd_w = np.take(wolf[0], wolf_westerly)
methc_w = np.take(wolf[3], wolf_westerly)

chfc_nw = np.take(cham[2], ch_northwesterly)
chwd_nw = np.take(cham[0], ch_northwesterly)
nefc_nw = np.take(neu[2], ne_northwesterly)
newd_nw = np.take(neu[0], ne_northwesterly)
methfc_nw = np.take(wolf[2], wolf_northwesterly)
methwd_nw = np.take(wolf[0], wolf_northwesterly)
methc_nw = np.take(wolf[3], wolf_northwesterly)

################################################################
ch_n = np.row_stack((chwd_n, chfc_n))
ne_n = np.row_stack((newd_n, nefc_n))
meth_n = np.row_stack((methwd_n,methfc_n, methc_n))

ch_ne = np.row_stack((chwd_ne, chfc_ne))
ne_ne = np.row_stack((newd_ne, nefc_ne))
meth_ne = np.row_stack((methwd_ne,methfc_ne, methc_ne))

ch_e = np.row_stack((chwd_e, chfc_e))
ne_e = np.row_stack((newd_e, nefc_e))
meth_e = np.row_stack((methwd_e,methfc_e, methc_e))

ch_se = np.row_stack((chwd_se, chfc_se))
ne_se = np.row_stack((newd_se, nefc_se))
meth_se = np.row_stack((methwd_se,methfc_se, methc_se))

ch_s = np.row_stack((chwd_s, chfc_s))
ne_s = np.row_stack((newd_s, nefc_s))
meth_s = np.row_stack((methwd_s,methfc_s, methc_s))

ch_sw = np.row_stack((chwd_sw, chfc_sw))
ne_sw = np.row_stack((newd_sw, nefc_sw))
meth_sw = np.row_stack((methwd_sw,methfc_sw, methc_sw))

ch_w = np.row_stack((chwd_w, chfc_w))
ne_w = np.row_stack((newd_w, nefc_w))
meth_w = np.row_stack((methwd_w,methfc_w, methc_w))

ch_nw = np.row_stack((chwd_nw, chfc_nw))
ne_nw = np.row_stack((newd_nw, nefc_nw))
meth_nw = np.row_stack((methwd_nw,methfc_nw, methc_nw))

#print np.nanmean(ne_w[1])
# x1 - x2
#np.subtract(x1, x2)

#### PLotting

nm = len(meth_n[2])
nem = len(meth_ne[2])
em = len(meth_e[2])
sem = len(meth_se[2])
sm = len(meth_s[2])
swm = len(meth_sw[2])
wm = len(meth_w[2])
nwm = len(meth_nw[2])
y = [nm, nem, em, sem, sm, swm, wm, nwm]
N = len(y)
x = range(N)
width = 1/1.5
fig1 = plt.bar(x, y, width, color="blue")
plt.show(fig1)
### Plot co2 flux in wind direction bins
nn = len(ne_n[1])
nen = len(ne_ne[1])
en = len(ne_e[1])
sen = len(ne_se[1])
sn = len(ne_s[1])
swn = len(ne_sw[1])
wn = len(ne_w[1])
nwn = len(ne_nw[1])
y = [nn, nen, en, sen, sn, swn, wn, nwn]
N = len(y)
x = range(N)
width = 1/1.5
fig2 = plt.bar(x, y, width, color="red")
plt.show(fig2)


### Plot time series of Methane Concentration
fig3 = plt.plot(methc_com[0], color = "blue")
#fig3 = plt.plot(methfc_com, color = 'red')


### Plot differences between co2 directionality dependencies between chamau and nestift
### Plot co2 flux diurnal cycle and wd flux diurnal cycle

sety = range(1,13)
setb = range(40,61)
days = range(0,17537)
for i in days:
    ch_ind = np.where((dtime == sety) & (dtime == (28*i + setb)))
    ch_night = np.take(chfc_com, ch_ind)

ch_night = np.where((cham[0] > 337.5) & (cham[0] < 22.5) & (cham[1] > 4))
ne_northerly = np.where((neu[0] > 337.5) & (neu[0] < 22.5) & (neu[1] > 4))
wolf_northerly = np.where((wolf[0] > 337.5) & (wolf[0] < 22.5) & (wolf[1] > 4))

    