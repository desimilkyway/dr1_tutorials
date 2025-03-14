import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import plot_preamb as pp
import scipy.optimize
import scipy.stats

pp.run()

RV_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
SP_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'SPTAB',
                         mask_invalid=False)
FM_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'FIBERMAP',
                         mask_invalid=False)
G_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                        'GAIA',
                        mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM']
                                                    == 'backup')


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']
HOST = open('WSDB', 'r').read()

plt.clf()
fig = plt.figure(figsize=(3.37 * 1., 3.37 * .7))
delt = RV_T['VRAD'] - G_T['RADIAL_VELOCITY']
cur_sel = np.isfinite(delt) & cur_sel0

xr = 0, 360
yr = -32, 90
bins = [360 // 2, 120 // 2]
S = scipy.stats.binned_statistic_2d(RV_T['TARGET_RA'][cur_sel],
                                    RV_T['TARGET_DEC'][cur_sel],
                                    delt[cur_sel],
                                    'median',
                                    range=[xr, yr],
                                    bins=bins)
SC = scipy.stats.binned_statistic_2d(RV_T['TARGET_RA'][cur_sel],
                                     RV_T['TARGET_DEC'][cur_sel],
                                     delt[cur_sel],
                                     'count',
                                     range=[xr, yr],
                                     bins=bins)
minval = 5
stat = S.statistic
stat[SC.statistic < minval] = np.nan
plt.imshow(stat.T,
           extent=list(xr) + list(yr),
           origin='lower',
           cmap='turbo',
           aspect='auto',
           vmax=25,
           vmin=-25)
plt.gci().set_rasterized(True)
plt.colorbar(label=r'$\delta V_{rad}$ [km/s]')
plt.xlim(360, 0)
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$\delta$ [deg]')
plt.tight_layout()
plt.savefig('plots/backup_delt.pdf')
