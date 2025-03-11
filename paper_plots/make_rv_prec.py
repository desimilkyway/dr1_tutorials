import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp

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
SC_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'SCORES',
                         mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .7))
cnt = 0
cnt = 0
cur_sel = (main_sel & (RV_T['SURVEY'] == 'main') &
           (RV_T['PROGRAM'] == 'bright') & (RV_T['VSINI'] < 30)
           & betw(RV_T['TEFF'], 4500, 7000)
           & betw(SC_T['TSNR2_LRG'] * 12.15, 180, 220))
zmag = 22.5 - 2.5 * np.log10(FM_T['FLUX_Z'])

zgrid = np.linspace(15.5, 19, 100)
feh_mult = -0.22

for i in range(2):

    plt.subplot(1, 2, i + 1)
    feh = {0: 0, 1: -2}[i]
    cur_sel1 = cur_sel & betw(RV_T['FEH'] - feh, -0.05, 0.05)
    plt.hist2d(
        zmag[cur_sel1],
        np.log10(RV_T['VRAD_ERR'][cur_sel1]),
        range=[[15.5, 19.5], [-.75, 1.2]],
        bins=[60, 60],
        # norm=maco.PowerNorm(gamma=.5)
    )
    plt.plot(zgrid,
             feh * feh_mult - 0.5 + (zgrid - 16) * 0.3,
             color='red',
             linestyle='--',
             label=r'$y= -0.5 + 0.3 \cdot (z-16)$' + '\n' +
             '$ %g [Fe/H]$' % feh_mult)
    # norm=maco.LogNorm()

    plt.xlabel('z [mag]')
    plt.title('[Fe/H]=%d' % feh)
    plt.xlim(15.5, 19.5)

    plt.gci().set_rasterized(True)
    if i == 0:
        plt.ylabel(r'$\log_{10} [\sigma_{RV}/(1 {\rm km/s})]$ ')
    else:
        plt.legend(numpoints=1, loc='lower right')
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())

# plt.subplot(122)
# plt.title('[Fe/H]=-2')
# cur_sel2 = cur_sel & betw(RV_T['FEH'], -2, -1.9)
# plt.xlabel('z [mag]')
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, left=.15)

plt.savefig('plots/rv_prec.pdf')
