#!/usr/bin/env python
"""
Given coordinates and a redshift for a source, an image, and a lens
model in the form of the deflections calculated at a specific
reference source redshift, construct the source plane model. 

L. Moustakas, JPL/Caltech, 11
May 2012 CLASH discrete, do not circulate outside collaboration.
"""
import matplotlib.pyplot as plt
import numpy as np
# import scipy.interpolate as si
import pyfits
import os
import asciitable
import cosmolopy as cp
import cosmolopy.distance as cd

cdad = cd.angular_diameter_distance
cosmo = {'omega_M_0' : 0.27,
         'omega_lambda_0' : 0.73,
         'omega_k_0': 0.0,
         'h' : 0.70}
extras = {'w' : -1.0}
cosmo.update(extras)

# read in setup information from a specific file
lmdir='/Users/leonidas/research/projects/CLASH/lensedvolumes/data/lensmodels/'
clist=asciitable.read(lmdir+'lensmodels.lis',
                      names=['alias','clusterdir','deflxfile','deflyfile','zcluster','zmodel'])
# <loop over all bands>
# * read in image plane cutout in an area that encompasses the full mask
# * read in deflection data over the same area. 
# * renormalize it to the appropriate source redshift
# * calculate the corresponding magnification map
# * calculate the full range of area needed in the source plane and create a zero array
"""
  * Now the tricky bit. We want to use the deflection and
  magnification to properly populate the source plane.  We can get
  fancy, but let's keep it simple.
   o Begin with identifying each image plane pixel.
   o Record the pixel's bounding coordinates [which coordinate system?
   alpha-dec?]
   o Transform each of those points to the source plane. This is the
   source-projected pixel. [What coordinate system? If still in
   alpha-dec, we will have managed to avoid any gridding choices
   still.]
   o

   http://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
   http://matplotlib.sourceforge.net/faq/howto_faq.html#test-whether-a-point-is-inside-a-polygon
  
"""

#   use deflection and magnification to populate the source plane image
#   -- can we drizzle properly... 

# Observed F160W AB magnitude is 25.7. Its restframe B-band magnitude
# is approximately -22.4, uncorrected for the magnification of ~15.5.
# So the corrected absolute magnitude is -22.4+2.5*log(15.5) = -19.4.

# target source redshift for rescaling model deflection maps
ztarget = np.linspace(2.5,11.5,10)

# a reasonable grid of magnifications
mu=np.linspace(0.2,120,100)
# fiducial threshold magnification
mufiducial = 1.0 
mustr = '%.1f' % mufiducial

# set the anchor point of Mref <-> muref
Mref, muref = -19.4, 15.5
Mlim = Mref+2.5*np.log10(mu/muref)
def getmu(Mag):
    return muref*10**(0.4*(Mag-Mref))

TotalVolarray = np.zeros(mu.shape)
TotalRedarray = np.zeros(ztarget.shape)

# Plot stuff
fig=plt.figure()
MFig=fig.add_subplot(111)
MFig.set_xlabel('Redshift [z]')
MFig.set_ylabel('Effective Volume ($\mu>$'+mustr+') [Mpc$^3$]')
MFig.set_xlim(2.0,12)
MFig.set_ylim(1.0,2e5)
MFig.set_yscale('log')

# # some annotations
# MFig.text(30,2e4,'z=[9,10]',size=13)
# ytext = 1e5
# MFig.text(0.22,ytext,'Unlensed M$_{B}$ [AB] limit:',size=10)
# plotmus=[-18.0,Mref,-21.0,-22.0]
# for pp in plotmus: MFig.text(getmu(pp)/1.2,ytext,'%.1f' % pp,size=10)

# outdata = open('redshiftdata.dat','w')

for ii in np.arange(clist.size):

    alias      = clist['alias'][ii]
    clusterdir = clist['clusterdir'][ii]
    deflxfile  = clist['deflxfile'][ii]
    deflyfile  = clist['deflyfile'][ii]
    zcluster   = clist['zcluster'][ii]
    zmodel     = clist['zmodel'][ii]

# read in the deflections from Adi's lens models 
    axlist = pyfits.open(clusterdir+'/'+deflxfile)
    aylist = pyfits.open(clusterdir+'/'+deflyfile)
#    ax, ay = rescale*axlist[0].data, rescale*aylist[0].data
    ax, ay = axlist[0].data, aylist[0].data

# do some initializations etc
    Unity  = np.zeros(ax.shape)+1.0
    header = axlist[0].header
    try:
        cdelt1, cdelt2 = header['CDELT1'], header['CDELT2']
        header['CDELT1'], header['CDELT2'] = cdelt1, cdelt2
        pxsc = np.abs(cdelt1)*3600.0 # in arcseconds per linear pixel dimension
    except:
        cd1_1, cd2_2 = header['CD1_1'], header['CD2_2']
        pxsc = np.abs(cd1_1)*3600.0 # in arcseconds per linear pixel dimension

    Redarray=np.zeros(ztarget.shape)
    for zi in np.arange(ztarget.size):
        print alias, ztarget[zi]
        rescale = \
            (cdad(ztarget[zi], zcluster, **cosmo) / cdad(ztarget[zi], **cosmo)) * \
            (cdad(zmodel, **cosmo) / cdad(zmodel, zcluster, **cosmo))
        
# Calculate the magnification from the Jacobian.  Note that it first
# calculates the gradient with respect to the 'y' axis, and then wrt
# the 'x' axis. 
        axy, axx = np.gradient(rescale*ax)
        ayy, ayx = np.gradient(rescale*ay)
        Jxx = Unity -axx
        Jxy =       -axy
        Jyx =       -ayx
        Jyy = Unity -ayy
    
        Mu    = Unity / (Jxx*Jyy - Jxy*Jyx)
        AbsMu = np.abs(Mu)

# The diff_comoving_volume is the differential comoving volume per
# unit redshift per unit solid angle, in units of Mpc**3 ster**-1.  So
# we convert that to the Volume per (unlensed) pixel, for a source
# plane at z=9-10.
        VolPix = (pxsc/206265.0)**2 * cd.diff_comoving_volume(ztarget[zi], **cosmo)

# PixMap will now be the Mpc**3 volume that each pixel corresponds to
# in the z=9-10 source plane.
        
        PixMap = VolPix / AbsMu 

        Redarray[zi]=np.sum(PixMap[(AbsMu>mufiducial)])

    TotalRedarray += Redarray
    MFig.plot(ztarget,Redarray,label=alias)
    MFig.legend(loc=3,prop={'size':8})
        
# # Now let us calculate the source-plane volume for each of the
# # magnification lower limits we established.
#         
#         Volarray = np.zeros(mu.shape)
#         for jj in np.arange(mu.size):
#             Volarray[jj] = np.sum( PixMap[(AbsMu>mufiducial)] )
#             TotalVolarray += Volarray
# 
#             s=si.UnivariateSpline(mu,Volarray,s=5)
#             print alias,clusterdir,s(muref)
#             outstring = '%s %s %.1f ' % (alias,clusterdir,s(muref))
#             outdata.write(outstring+'\n')
#             MFig.plot(mu,Volarray,label=alias)
#             MFig.legend(loc=3,prop={'size':8})

# s=si.UnivariateSpline(mu,TotalVolarray,s=5)
# print 'Total Volume at WeiOne Magnification / Luminosity: ',s(muref)
# outstring = 'Total Total %.1f ' % (s(muref))
# outdata.write(outstring+'\n')
# outdata.close()

MFig.plot(ztarget,TotalRedarray,linewidth=4,color='b',label='CLASH 12')
MFig.legend(loc=3,prop={'size':10})
#for pp in plotmus: MFig.text(getmu(pp)/1.2,ytext,'%.1f' % pp,size=10)

# MFig.plot(mu,TotalVolarray,linewidth=4,color='b',label='CLASH 12')
# MFig.legend(loc=3,prop={'size':10})
# for pp in plotmus: MFig.text(getmu(pp)/1.2,ytext,'%.1f' % pp,size=10)

fig.savefig('masterredshift.png')

