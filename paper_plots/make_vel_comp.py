import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp
import crossmatcher
from matplotlib.colors import TABLEAU_COLORS
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
    aind = main_sel & np.isfinite(X + Y) & betw(teff, 4500, 7000)
    X, Y = X[aind], Y[aind]
    R1 = scipy.optimize.minimize(func, [0, 0, 0], args=(X, Y))
    return R1.x


main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
cnt = 0
# cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
#    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['SN_R'] > 10)


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']
HOST = open('WSDB', 'r').read()
D_GA = crossmatcher.doit(
    'galah_dr4.allstar',
    ra,
    dec,
    'fe_h,teff,logg,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp,rv_comp_1',
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

minteff, maxteff = 4500, 7000
from matplotlib.colors import TABLEAU_COLORS, same_color

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .7))
cnt = 0
for jj, program in enumerate(['bright', 'backup']):
    plt.subplot(2, 1, cnt + 1)
    xind = (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == program) & main_sel
    for ii, xv in enumerate(
        [G_T['RADIAL_VELOCITY'], D_AP['vhelio_avg'], D_GA['rv_comp_1']]):
        delt = (RV_T['VRAD'] - xv)[xind]
        lab = ['Gaia', 'APOGEE', 'GALAH'][ii]
        print(program, lab, np.nanmedian(delt))
        if np.isfinite(delt).sum() < 50:
            delt = [-1000]
        plt.hist(delt,
                 range=[-15, 15],
                 bins=60,
                 histtype='step',
                 density=True,
                 color=['#1f77b4', '#ff7f0e', '#2ca02c'][ii],
                 label=lab)
    cnt += 1
    if cnt == 1:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
        plt.legend()
    plt.text(-15, [.25, .13][jj], 'survey=' + program)
plt.xlabel('RV$_{DESI}$ - RV$_{ref}$ [km/s]')
plt.tight_layout()
plt.subplots_adjust(hspace=0.)
plt.savefig('plots/rv_comp.pdf')
