import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp

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

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .75))
cnt = 0
cnt = 0
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)

for cnt in range(2):
    T = [RV_T, SP_T][cnt]
    if cnt == 0:
        cur_sel = cur_sel0
    else:
        cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1')
    plt.subplot(1, 2, cnt + 1)
    plt.hist2d((T['FEH'][cur_sel]),
               T['ALPHAFE'][cur_sel],
               range=[[-5, 1], [-1, 1.3]],
               bins=[100, 100],
               norm=maco.LogNorm())

    plt.gca().set_rasterized(True)
    plt.xlabel('[Fe/H] ')

    plt.title(['RVS', 'SP'][cnt])
    if cnt == 0:
        plt.ylabel(r'$\alpha$/Fe')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.01)

plt.savefig('plots/feh_alpha.pdf')
