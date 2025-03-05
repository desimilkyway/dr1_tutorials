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

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .7))
cnt = 0
cnt = 0
cur_sel0 = (
    main_sel & (RV_T['SURVEY'] == 'main') &
    # (    RV_T['PROGRAM'] == 'bright') &
    (RV_T['SN_R'] > 10))

for cnt in range(2):
    T = [RV_T, SP_T][cnt]
    if cnt == 0:
        cur_sel = cur_sel0
    else:
        cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1')

    plt.subplot(1, 2, cnt + 1)
    plt.hist2d(np.log10(T['TEFF'][cur_sel]),
               T['LOGG'][cur_sel],
               range=[[3.45, 4.2], [-0.5, 6]],
               bins=[100, 100],
               norm=maco.LogNorm())

    plt.gci().set_rasterized(True)
    plt.xlim(4.19, 3.45)
    plt.ylim(6, -0.5)
    plt.xlabel('$\log_{10}$ T$_{eff}$ ')

    plt.text(4, 0, ['RVS', 'SP'][cnt])
    if cnt == 0:
        plt.ylabel('$\log$ g')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, hspace=0.01, top=.99)

plt.savefig('plots/logg_teff.pdf')
