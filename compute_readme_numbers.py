import astropy.table as atpy
import numpy as np
import glob
import astropy.io.fits as pyfits

T = atpy.Table().read('data/mwsall-pix-iron.fits')
TG = atpy.Table().read('data/mwsall-pix-iron.fits', 'GAIA')
print(len(T))
#    The total number of fitted coadded spectra included in the release: 6,372,607
qq = T[T['RR_SPECTYPE'] == 'STAR']
print(len(np.unique(qq['TARGETID'])))
#    The number of unique DESI sources, classified as stars by Redrock: 4,586,247

print(len(np.unique(TG['SOURCE_ID'])))
#    The number of unique Gaia sources: 4,920,596
fs = glob.glob('../rv_variability/rvtabs_iron/rvpix_exp-*')
print(sum([len(pyfits.getdata(_)) for _ in fs]))
#    The number of single-epoch spectra with stellar atmospheric parameters and radial velocities: 10,012,925
R = {}
for f in fs:
    TX = atpy.Table().read(f, 'GAIA', format='fits')
    for s in TX['SOURCE_ID']:
        if s > 0:
            if s not in R:
                R[s] = 0
            R[s] += 1
cnt = 0
for k, v in R.items():
    if v > 1:
        cnt += 1
print(cnt)
#    The number of Gaia sources with more than one radial velocity/stellar parameter measurement: 1,718,305
print('backup', ((T['SURVEY'] == 'main') & (T['PROGRAM'] == 'backup')).sum())
print('bright', ((T['SURVEY'] == 'main') & (T['PROGRAM'] == 'bright')).sum())
print('dark', ((T['SURVEY'] == 'main') & (T['PROGRAM'] == 'dark')).sum())
#    The largest number of sources in the catalog has been observed in the bright program of the main survey (3,070,120), followed by the dark program of the main survey (1,357,709), and the backup program of the main survey (1,218,152)
