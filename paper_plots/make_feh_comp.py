import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp
import crossmatcher
import match_lists
from matplotlib.colors import TABLEAU_COLORS
import scipy.optimize


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


def combiner(*args):
    # fill nansby us
    ret = args[0] * 1
    for x in args[1:]:
        ind = ~np.isfinite(ret)
        ret[ind] = x[ind]
    return ret


teff_ref = 5000
logteff_scale = 0.1


def func(p, X, Y):
    return np.mean(
        np.abs(Y - np.poly1d(p)((X - np.log10(teff_ref)) / logteff_scale)))


def fitter(teff, feh_ref, feh_obs):
    X, Y = np.log10(teff), feh_obs - feh_ref
    aind = main_sel & np.isfinite(X + Y) & betw(teff, minteff, maxteff)
    X, Y = X[aind], Y[aind]
    R1 = scipy.optimize.minimize(func, [0, 0, 0], args=(X, Y))
    return R1.x


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
cnt = 0

# cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
#    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)

cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['SN_R'] > 10)

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

D_GA = crossmatcher.doit('galah_dr4.allstar',
                         ra,
                         dec,
                         'fe_h,teff,logg,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp',
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

coeff_sp = fitter(SP_T['TEFF'], combiner(D_GA["fe_h"], D_AP['fe_h']),
                  SP_T['FEH'])
coeff_rv = fitter(RV_T['TEFF'], combiner(D_GA["fe_h"], D_AP['fe_h']),
                  RV_T['FEH'])
coeff_sp = np.round(coeff_sp, 3)
coeff_rv = np.round(coeff_rv, 3)

SP_T['FEH_CALIB'] = SP_T['FEH'] - np.poly1d(coeff_sp)(
    np.log10(SP_T['TEFF'] / teff_ref) / logteff_scale)
RV_T['FEH_CALIB'] = RV_T['FEH'] - np.poly1d(coeff_rv)(
    np.log10(RV_T['TEFF'] / teff_ref) / logteff_scale)
print('RV', coeff_rv[::-1], 'SP', coeff_sp[::-1])

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
    curfeh = T["FEH"]
    if cnt % 3 < 2:
        plt.hist2d(COMP['fe_h'][cur_sel],
                   curfeh[cur_sel],
                   range=[[-3, 1], [-3, 1]],
                   bins=[50, 50],
                   norm=maco.PowerNorm(gamma=0.5))
        plt.plot([-3, 1], [-3, 1], color='black', linestyle='--')
        plt.gci().set_rasterized(True)
        plt.text(-2, 0.5, ['GALAH', 'APOGEE'][cnt % 3], color='white')
    else:
        plt.plot(COMP['fe_h'][cur_sel], curfeh[cur_sel], '.')
        plt.xlim(-4, 0)
        plt.ylim(-4, 0)
        plt.plot([-4, 1], [-4, 1], color='black', linestyle='--')
        plt.text(-3, -.5, 'SAGA')

    plt.xlabel('[Fe/H]$_{survey}$ [dex]')

    # plt.title(titl)
    if cnt % 3 == 0:
        plt.annotate(titl, (0.1, .9), xycoords='axes fraction', color='white')
        plt.ylabel(r'[Fe/H]$_{DESI}$ [dex]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0.04, hspace=0.25)
plt.savefig('plots/feh_compar.pdf')

for var_name in ['FEH', 'FEH_CALIB']:
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
            curfeh = T[var_name]

            cur_sel = cur_sel & betw(T['TEFF'], minteff,
                                     maxteff)  # & (T['FEH_ERR'] < .1)
            delt = curfeh[cur_sel] - COMP['fe_h'][cur_sel]
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
        postfix = {'FEH': '', 'FEH_CALIB': r'$_{\rm calibrated}$'}[var_name]
        plt.xlabel(r'$\delta$ [Fe/H]' + postfix + ' [dex]')
        if cnt == 0:
            plt.legend()
    plt.tight_layout()
    # plt.subplots_adjust(wspace=0.04, hspace=0.25)
    plt.savefig({
        'FEH': 'plots/feh_compar_delta.pdf',
        'FEH_CALIB': 'plots/feh_compar_delta_calibrated.pdf'
    }[var_name])
