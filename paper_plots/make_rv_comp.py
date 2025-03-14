import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import plot_preamb as pp
import crossmatcher
import scipy.optimize

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


def combiner(*args):
    # fill nans
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
    aind = main_sel & np.isfinite(X + Y) & betw(teff, 4500, 7000)
    X, Y = X[aind], Y[aind]
    R1 = scipy.optimize.minimize(func, [0, 0, 0], args=(X, Y))
    return R1.x


main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
cnt = 0
# cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
#    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main')  # & (RV_T['SN_R'] > 10)


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']
HOST = open('WSDB', 'r').read()
if False:
    D_GA = crossmatcher.doit(
        'galah_dr4.allstar',
        ra,
        dec,
        'fe_h,teff,logg,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp,rv_comp_1',
        host=HOST,
        db='wsdb',
        asDict=True)
else:
    D_GA = crossmatcher.doit_by_key(
        'galah_dr4.allstar',
        G_T['SOURCE_ID'],
        'fe_h,teff,logg,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp,rv_comp_1',
        key_col='gaiadr3_source_id',
        host=HOST,
        db='wsdb',
        asDict=True)
D_GA['rv_comp_1'][D_GA['flag_sp'] != 0] = np.nan

if False:
    D_AP = crossmatcher.doit('apogee_dr17.allstar',
                             ra,
                             dec,
                             '''vhelio_avg,starflag''',
                             host=HOST,
                             db='wsdb',
                             asDict=True)
else:
    D_AP = crossmatcher.doit_by_key('apogee_dr17.allstar',
                                    G_T['SOURCE_ID'],
                                    '''vhelio_avg,starflag''',
                                    key_col='gaiaedr3_source_id',
                                    host=HOST,
                                    db='wsdb',
                                    asDict=True)
D_AP['vhelio_avg'][D_AP['starflag'] != 0] = np.nan

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .7))
cnt = 0
for jj, program in enumerate(['bright', 'backup', 'dark']):
    plt.subplot(3, 1, cnt + 1)
    xind = (RV_T['SURVEY'] == 'main') & (
        RV_T['PROGRAM'] == program) & main_sel & (RV_T['VRAD_ERR'] < 1)
    rvs = [G_T['RADIAL_VELOCITY'], D_AP['vhelio_avg'], D_GA['rv_comp_1']]
    labs = ['Gaia', 'APOGEE', 'GALAH']
    if program == 'dark':
        aind = (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'bright') & (
            RV_T['VRAD_ERR'] < 1) & main_sel
        lookup = dict(zip(RV_T['TARGETID'][aind], RV_T['VRAD'][aind]))
        rvs = [np.array([lookup.get(x, np.nan) for x in RV_T['TARGETID']])]
        labs = ['DESI,bright']
    for ii, xv in enumerate(rvs):
        cur_xind = xind
        lab = labs[ii]
        if lab == 'Gaia':
            cur_xind = xind & (G_T['RADIAL_VELOCITY_ERROR'] < 5)
        delt = (RV_T['VRAD'] - xv)[cur_xind]
        print(program, lab, np.nanmedian(delt))
        if np.isfinite(delt).sum() < 50:
            delt = [-1000]
        if program != 'dark':
            col = ['#1f77b4', '#ff7f0e', '#2ca02c'][ii]
        else:
            col = '#d62728'
        plt.hist(delt,
                 range=[-15, 15],
                 bins=60,
                 histtype='step',
                 density=True,
                 color=col,
                 label=lab)
    cnt += 1
    if cnt == 1:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
    if cnt in [1, 3]:
        plt.legend()
    plt.text(-15, [.15, .08, .1][jj], 'survey=main\n' + 'program=' + program)
    plt.axvline(0, linestyle='--', color='black')
plt.xlabel('RV$_{DESI}$ - RV$_{ref}$ [km/s]')
plt.tight_layout()
plt.subplots_adjust(hspace=0.)
plt.savefig('plots/rv_comp.pdf')
