import numpy as np
from netCDF4 import Dataset
import sys

#------------------read NO2 from netcdf files-------------------------
filename = '/scratch/d6/groupYim/b137225/adjoint_output2/base/CCTM_e2a.CONC.e2a.19990702'
readfile = Dataset(filename, mode='r', open=True)
readVar  = np.squeeze(readfile.variables['NO2'][:][:][:][:]) 
TFLAG    = np.squeeze(readfile.variables['TFLAG'][:][:][:])
#print(readVar.shape)
getVar   = readVar[24,0,:,:] # layer 1 and last tstep
getVar   = getVar * 1000     # ppmv to ppbv
col = getVar.shape[1]
row = getVar.shape[0]
dim = col * row
Max = np.amax(getVar)
Min = np.amin(getVar)
shapeVar = np.reshape(getVar,dim)

#------------------------------Apply-KF-filter---------------------------
from gh_internal import plot_g_h_results
import book_plots
import numpy as np
import matplotlib.pyplot as plt

def g_h_filter(data,x0,dx,g,h,dt=1.,pred=None):
    x=x0
    results=[]
    for z in data:
        x_est=x+(dx*dt)
        dx=dx
        if pred is not None:
            pred.append(x_est)
        residual=z-x_est
        dx=dx+h*(residual)/dt
        x=x_est+g*residual
        results.append(x)
    return np.array(results)

weights = shapeVar
book_plots.plot_track([0, dim], [Min, Max], label='Actual weight')
data = g_h_filter(data=weights, x0=1, dx=1, g=3./10, h=1./3, dt=1.)
#plot_g_h_results(weights, data);
#plt.show()
#print(np.amin(data))

#-----------------------further-adjustment---------------------------
# for negative conc values
CMIN  = 0
i = 0
for index in range(0,dim):
    if ( data[i] < CMIN and i < dim):
        data[i] = CMIN
#    print(data[i])
    i = i + 1

#----------------------Dump-var-to-new-netcdf-file-------------------
dumpVar = np.reshape(data,(row,col))
varOut  = dumpVar/1000  # ppbv to ppmv

ltime,llay = (1,1)
lrow, lcol = (row,col)
nvars  = 1
ldattim = 2

attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]

datapath = 'varoutput.nc'
varDump = Dataset(datapath, 'w', format='NETCDF3_64BIT')

for attr in attrs:
    # Read the attribute
    if hasattr(readfile, attr):
        attrVal = getattr(readfile, attr);
        # Write to the output file
        setattr(varDump, attr, attrVal)

setattr(varDump,"NVARS", nvars)
setattr(varDump,"NLAYS", llay)
setattr(varDump,"UPNAM","CALFILE")
setattr(varDump,"VAR-LIST","NO2             ")
setattr(varDump,"FILEDESC","output variables after applying kalman filter")

# create dimensions
varDump.description = 'read variables from checkpoint 5'
varDump.createDimension('TSTEP', ltime)
varDump.createDimension('DATE-TIME', ldattim)
varDump.createDimension('LAY', llay)
varDump.createDimension('VAR', 1)
varDump.createDimension('ROW', lrow)
varDump.createDimension('COL', lcol)

varDump.sync()

dattim = np.squeeze(readfile.variables['TFLAG'][:][:])

#create variables
file_tflag = varDump.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
file_tflag[:] = np.zeros([2,ltime])
varattrs=["long_name","units","var_desc"]
for varattr in varattrs:
    # Read the attribute
    if hasattr(readfile.variables['TFLAG'], varattr):
        varattrVal = getattr(readfile.variables['TFLAG'], varattr);
        # Write it to the new file
        setattr(file_tflag, varattr, varattrVal)

for tt in range(0,ltime):
    for vv in range(0,nvars):
        for dd in range(0,ldattim):
            file_tflag[tt,vv,dd]  =  TFLAG[tt,vv,dd]

Var     = varDump.createVariable('NO2','f4',('TSTEP','LAY','ROW','COL'))
varattrs=["long_name","units","var_desc"]
Var[:]  = np.zeros([ltime,llay,lrow,lcol])
Var[:]  = np.reshape(varOut,[ltime,llay,lrow,lcol])

for varattr in varattrs:
    # Read the attribute
    if hasattr(readfile.variables['NO2'], varattr):
        varattrVal = getattr(readfile.variables['NO2'], varattr);
        # Write to the new file
        setattr(Var, varattr, varattrVal)

varDump.close()


