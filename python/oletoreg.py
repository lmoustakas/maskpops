#!/usr/bin/env python
import numpy as np
import sys

if __name__ == "__main__":
    datafile=sys.argv[1]
else:
    datafile='arcs_phz_fromole.dat'
data = np.loadtxt(datafile)
id,rah,ram,ras,dd,dm,ds=data.T

regfile=datafile.split('.dat')[0]+'.reg'
r=open(regfile,'w')
r.write('global color=purple\n')

for ii in np.arange(rah.size):
    hms='%.2i'%rah[ii]+':'+\
        '%.2i'%ram[ii]+':'+\
        '%.3f'%ras[ii]
    dms='%.2i'%dd[ii]+':'+\
        '%.2i'%dm[ii]+':'+\
        '%.2i'%np.trunc(ds[ii])+'.'+'%.3i'%(1000*(ds[ii]-np.trunc(ds[ii])))
    idstr = '%s' % id[ii]
    r.write('FK5;circle('+hms+','+dms+',1") # text={'+idstr+'}\n')
r.close()

#
#radeg=rah*15.0+ram/60.0+ras/3600.0
#decdeg=np.sign(dd)*(np.abs(dd)+dm/60.0+ds/3600.0)
#onearcsec = 1./3600.
#radstr = '%.6f' % onearcsec
#for ii in np.arange(radeg.size):
#    radegstr = '%.8f' % radeg[ii]
#    decdegstr = '%.6f' % decdeg[ii]
#    idstr = '%s' % id[ii]
#    r.write('circle('+radegstr+','+decdegstr+','+radstr+') # text={'+idstr+'}\n')
#r.close()
