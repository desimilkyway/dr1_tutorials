import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp
import crossmatcher
import match_lists
from matplotlib.colors import TABLEAU_COLORS

pp.run()

RV_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
SP_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'SPTAB',
                         mask_invalid=False)
FM_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'FIBERMAP',
                         mask_invalid=False)
G_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                        'GAIA',
                        mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
cnt = 0
# cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
#    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)
cur_sel0 = main_sel & (RV_T['SN_R'] > 10)


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']
SAGAT = atpy.Table().read(
    'saga_cleaned_catalog.tsv',
    format='ascii',
)
DD, xind = match_lists.match_lists(ra, dec, SAGAT['RAdeg'], SAGAT['DECdeg'],
                                   1. / 3600)
SAGA_R = {'fe_h': np.zeros(len(ra)) + np.nan}
SAGA_R['fe_h'][np.isfinite(DD)] = SAGAT['[M/H]'].filled(
    np.nan)[xind[np.isfinite(DD)]]
HOST = open('WSDB', 'r').read()
D_GA = crossmatcher.doit(
    'galah_dr4.allstar',
    ra,
    dec,
    'fe_h,teff,logg,fe_h,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp',
    host=HOST,
    db='wsdb',
    asDict=True)

D_AP = crossmatcher.doit(
    'apogee_dr17.allstar',
    ra,
    dec,
    '''alpha_m,fe_h,c_fe,n_fe,o_fe,na_fe,mg_fe,si_fe,ca_fe,ti_fe,mn_fe,ni_fe,ce_fe,vhelio_avg,logg,teff,teff_spec,logg_spec,
    aspcapflag,starflag, fe_h_flag''',
    host=HOST,
    db='wsdb',
    asDict=True)
D_AP['fe_h'][(D_AP['fe_h_flag'] != 0) | (D_AP['starflag'] != 0) |
             (D_AP['aspcapflag'] != 0)] = np.nan
D_GA['fe_h'][(D_GA['flag_fe_h'] != 0) | (D_GA['flag_sp'] != 0)] = np.nan

minteff, maxteff = 4500, 7000
plt.clf()
fig = plt.figure(figsize=(3.37 * 2, 3.37 * 1.5))
for cnt in range(6):
    T = [RV_T, SP_T][cnt // 3]
    titl = ['RVS', 'SP'][cnt // 3]
    if cnt == 0:
        cur_sel = cur_sel0 & (T['FEH_ERR'] < 0.1) & (T['VSINI'] < 30)
    else:
        cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1') & (
            SP_T['COVAR'][:, 0, 0]**.5 < .1)
    cur_sel = cur_sel & betw(T['TEFF'], minteff, maxteff)
    plt.subplot(2, 3, cnt + 1)
    COMP = [D_GA, D_AP, SAGA_R][(cnt % 3)]

    if cnt % 3 < 2:
        plt.hist2d(COMP['fe_h'][cur_sel],
                   T['FEH'][cur_sel],
                   range=[[-3, 1], [-3, 1]],
                   bins=[50, 50],
                   norm=maco.LogNorm())
        plt.plot([-3, 1], [-3, 1], color='black')
        plt.gca().set_rasterized(True)
        plt.text(-2, 0.5, ['GALAH', 'APOGEE'][cnt % 3])
    else:
        plt.plot(COMP['fe_h'][cur_sel], T['FEH'][cur_sel], '.')
        plt.xlim(-4, 0)
        plt.ylim(-4, 0)
        plt.plot([-4, 1], [-4, 1], color='black')
        plt.text(-3, -.5, 'SAGA')

    plt.xlabel('[Fe/H]$_{survey}$ [dex]')

    #plt.title(titl)
    if cnt % 3 == 0:
        plt.annotate(titl, (0.1, .9), xycoords='axes fraction')
        plt.ylabel(r'[Fe/H]$_{DESI}$ [dex]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0.04, hspace=0.25)
plt.savefig('plots/feh_compar.pdf')

fig = plt.figure(figsize=(3.37 * 2, 3.37 * .7))
for cnt in range(3):
    plt.subplot(1, 3, cnt + 1)
    COMP = [D_GA, D_AP, SAGA_R][(cnt % 3)]
    tit = ['GALAH', 'APOGEE', 'SAGA'][cnt % 3]

    for i, (T, label) in enumerate(zip([RV_T, SP_T], ['RVS', 'SP'])):
        if i == 0:
            cur_sel = cur_sel0 & (T['FEH_ERR'] < .1) & (T['VSINI'] < 30)
        else:
            cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1') & (
                SP_T['COVAR'][:, 0, 0]**.5 < .1)
        cur_sel = cur_sel & betw(T['TEFF'], minteff,
                                 maxteff)  # & (T['FEH_ERR'] < .1)
        delt = T['FEH'][cur_sel] - COMP['fe_h'][cur_sel]
        delt = delt[np.isfinite(delt)]
        percs = [np.percentile(delt, _) for _ in [16, 50, 84]]
        med = percs[1]

        plt.hist(
            delt,
            range=[-1, 1],
            linestyle=[':', '--'][i],
            histtype='step',
            bins=[100, 100, 10][cnt],
            label=label,
        )
        plt.annotate('$%.2f_{%.2f}^{+%.2f}$' %
                     (percs[1], percs[0] - percs[1], percs[2] - percs[1]),
                     (.7, .8 - .1 * i),
                     color=list(TABLEAU_COLORS.values())[i],
                     xycoords='axes fraction')
    plt.annotate(tit, (.5, .9), xycoords='axes fraction')
    plt.xlabel(r'$\delta$ [Fe/H] [dex]')
    if cnt == 0:
        plt.legend()
plt.tight_layout()
#plt.subplots_adjust(wspace=0.04, hspace=0.25)
plt.savefig('plots/feh_compar2.pdf')
