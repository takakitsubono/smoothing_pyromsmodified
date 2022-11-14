# cording: utf-8
import os
import numpy as np
import netCDF4

#if pyroms has been installed
import pyroms
import pyroms_toolbox
from bathy_smoother import *
#if pyroms has been installed

from smlp  import smthng_laplacian_rx0
from smlp  import RoughMtrx
import time

grd_name = 'grid.nc'
nc = netCDF4.Dataset(grd_name,"r")
q=nc.variables['hraw'][:,:]
msk_rho=nc.variables['mask_rho'][:,:]
area=nc.variables['area'][:,:]
nc.close()
nc.close

print('### Max Roughness cal. ####')
print('=== original pyroms progmram ====')
st_time=time.perf_counter()
RoughMat_p = bathy_tools.RoughnessMatrix(q, msk_rho )
print ('Max Roughness value is: ', RoughMat_p.max())
ed_time=time.perf_counter()
pyrm_time=ed_time-st_time
print('Time_pyrms: ',pyrm_time)

print(' ')
print('=== converted program  ====')
st_time=time.perf_counter()
RoughMat_s = RoughMtrx(q, msk_rho )
print ('Max Roughness value is: ', RoughMat_s.max())
ed_time=time.perf_counter()
smlp_time=ed_time-st_time
print('Time_smlp : ',smlp_time)
print(' ')
print('Ratio(pyrms/smlp): ',pyrm_time/smlp_time)

print(' ')
rx0_max=0.3
print('### Smoothing cal.(rx=0.3) ####')
print('=== original pyroms progmram ====')
st_time=time.perf_counter()
bbb = bathy_smoothing.smoothing_Laplacian_rx0(msk_rho, q, rx0_max)
ed_time=time.perf_counter()
pyrm_time=ed_time-st_time
print('Time_pyrms: ',pyrm_time)

print(' ')
print('=== converted program  ====')
st_time=time.perf_counter()
aaa=smthng_laplacian_rx0(msk_rho,q,rx0_max)
ed_time=time.perf_counter()
smlp_time=ed_time-st_time
print('Time_smlp : ',smlp_time)
print(' ')
print('Ratio(pyrms/smlp): ',pyrm_time/smlp_time)

print('Num of diff. point :',np.size(np.where(aaa!=bbb),1))

print(' ')
rx0_max=0.15
print('### Smoothing  cal.(rx=0.15) ####')
print('=== original pyroms progmram ====')
st_time=time.perf_counter()
bbb = bathy_smoothing.smoothing_Laplacian_rx0(msk_rho, q, rx0_max)
ed_time=time.perf_counter()
pyrm_time=ed_time-st_time
print('Time_pyrms: ',pyrm_time)

print(' ')
print('=== converted program  ====')
st_time=time.perf_counter()
aaa=smthng_laplacian_rx0(msk_rho,q,rx0_max)
ed_time=time.perf_counter()
smlp_time=ed_time-st_time
print('Time_smlp : ',smlp_time)
print(' ')
print('Ratio(pyrms/smlp): ',pyrm_time/smlp_time)
print('Num of diff. point :',np.size(np.where(aaa!=bbb),1))






