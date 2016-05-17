import netCDF4
from netCDF4 import Dataset

filroot1 = 'C:\Users\James Simkins\Documents\Wisconsin\\2016\Bioclimatology\euroflux\EFDC_L2_Flx_ATNeu_'
filroot2 = '_v02_30m'

for i in range(2005,2006,1):
    i = str(i)
    filin = [filroot1+i+filroot2]
    Adata = open('filin', 'r')


    

#CHCha_path = 'C:\Users\James Simkins\Documents\Wisconsin\\2016\Bioclimatology\euroflux\EFDC_L2_Flx_CHCha_',i,'_v06_30m'

#EFDC_L2_Flx_ATNeu_2005_v02_30m