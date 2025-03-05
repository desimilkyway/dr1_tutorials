import astropy.table as atpy
import numpy as np

RV_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
# Colour-magnitude distribution by survey/program
surveys = np.unique(RV_T['SURVEY'][:])
mysurveys = ['cmx', 'special', 'sv1', 'sv2', 'sv3', 'main']
assert (tuple(sorted(mysurveys)) == tuple(surveys))
for s in mysurveys:
    xind = RV_T['SURVEY'] == s
    for p in np.unique(RV_T['PROGRAM'][xind]):
        xind1 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM'] == p)
        xind2 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM'] == p) & (
            RV_T['RR_SPECTYPE'] == 'STAR')
        xind3 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM']
                                         == p) & (RV_T['RVS_WARN'] == 0)
        print('%s & %s &  %s & %s & %s \\\\' %
              (s, p, format(xind1.sum(), ','), format(
                  xind2.sum(), ','), format(xind3.sum(), ',')))
