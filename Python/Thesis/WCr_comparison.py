# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from netCDF4 import Dataset
import pandas as pd
import numpy as np

#Define function to get variables
def get_var(sds, name):
    return np.array(sds.variables[name])


gfdl = Dataset('Data/GFDL.CM3.rcp45.r1i1p1.2006.nc')
print gfdl.variables.keys()

dap = ["tas","rlds","ps","rsds","uas","vas","huss","pr"]
var = ["air_temperature","surface_downwelling_longwave_flux_in_air","air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"]

emp = list()
for i in xrange(len(var)):
    z = get_var(gfdl,var[i])
    #reshape removes the extra two dimensions (i.e. z.shape = (2920L, 1L, 1L,) so it removes the 1L's)
    z = z.reshape(z.shape[0])
    emp.append(z)

dates = pd.date_range('2006-01-01 00:00', periods=2920, freq='3H')

# writing the variable arrays to pd.series with the index dates so we can make the dataframe and then resample
ta = pd.Series(emp[0], index = dates)
rlds = pd.Series(emp[1], index = dates)
ps = pd.Series(emp[2], index = dates)
ps = ps/100
rsds = pd.Series(emp[3], index = dates)
pr = pd.Series(emp[7], index = dates)
df = pd.DataFrame({'ta' : ta, 'rlds' : rlds, 'ps' : ps, 'rsds' : rsds, 'pr' : pr})

dates2 = pd.date_range('2006-12-31 21:30:00', periods=5, freq='3Min')
df2 = pd.DataFrame([[np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]], columns = ('ta', 'rlds', 'ps', 'rsds', 'pr'), index = dates2)
#resample the dataframe to 30 minutes so we can compare this with willow creek flux data. as.freq() fills new values with nan
gfdl = df.resample('30Min').asfreq()
gfdl = gfdl.append(df2)
###################################Now read in the flux data ##################################

flux = pd.read_csv('Data/FLX_US-WCr_FLUXNET2015_SUBSET_HH_1999-2014_1-1.csv')
flux_dates = pd.date_range('1999-01-01 00:00', periods=280512, freq='30Min')

ta_ = pd.Series(flux['TA_F'].values, index = flux_dates)
ta_ = ta_ + 273.15
rlds_ = pd.Series(flux['LW_IN_F'].values, index = flux_dates)
ps_ = pd.Series(flux['PA_F'].values, index = flux_dates)
ps_ = ps_*10
rsds_ = pd.Series(flux['SW_IN_F'].values, index = flux_dates)
pr_ = pd.Series(flux['P_F'].values, index = flux_dates)

tower = pd.DataFrame({'ta' : ta_, 'rlds' : rlds_, 'ps' : ps_, 'rsds' : rsds_, 'pr' : pr_})
tower = tower['2006']

tower['ta'].plot.hist(stacked=True, bins=1000)