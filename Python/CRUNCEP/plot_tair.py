from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from copy import copy
import matplotlib

def get_var(sds, name):
    return np.array(sds.variables[name])

pathname= 'http://thredds.daac.ornl.gov/thredds/dodsC/ornldaac/1220/mstmip_driver_global_hd_climate_tair_2008_v1.nc4?tair[0:1:20][0:1:359][0:1:719],lat[0:1:359],lon[0:1:719],time[0:1:20]'
filehandle=Dataset(pathname,'r',format="NETCDF4")

## [scan lines along track direction, pixel elements along scan direction]
temp11 = get_var(filehandle, 'tair')[20]
lats = get_var(filehandle, 'lat')
lons = get_var(filehandle, 'lon')

print 'getting data'

# super basic image display - no coords, maps, etc.
#plt.imshow(temp11[0:500,:])
#plt.colorbar()
#plt.show()

# attempt a nice mapped image
cmap = copy(matplotlib.cm.get_cmap('gray'))
plt.figure(figsize=(20,9),dpi=240)
lons, lats = np.meshgrid(lons,lats)
m = Basemap(projection='cyl',lon_0=0)
m.drawcoastlines(color='aqua', linewidth=0.5)
m.drawcountries(color='aqua', linewidth=0.5)
print 'done drawing lines'
im = m.pcolormesh(lons,lats,temp11,latlon=True,vmin=210,vmax=310)
print 'done with pcolor'
plt.title('6 Hourly Instantaneous Air Temperature at 2m For 2008')
plt.colorbar()
#plt.savefig('globalmet.png', bbox_inches='tight')
plt.show()
print 'done saving'

