#!/usr/bin/env python
import numpy as np
import pyregion
import sys

if __name__ == "__main__":
    datafile=sys.argv[1]
else:
    datafile='macs1206_phzrois_green.reg'

f=open(datafile,'r')
for line in f:
    try:
        tagname=(line.split('text={')[1]).split('}')[0]
        regfile=datafile.split('.reg')[0]+'_'+tagname+'.reg'
        o=open(regfile,'w')
        o.write('FK5\n')
        o.write(line)
        o.close()
    except:
        continue
f.close()

# startfile=datafile.split('.reg')[0]+'_startup.dat'
# s=open(startfile,'w')
# s.write('# 1 rootname\n')
# s.write('# 2 redshift\n')
# s.write('# 3 cluster\n')
# s.write('# 4 datestamp\n')
# s.write('# 5 scaledir\n')
# s.close()
