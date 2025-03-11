import astropy.table as atpy
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
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .65))
cnt = 0
cnt = 0
cur_sel0 = (
    main_sel & (RV_T['SURVEY'] == 'main') &
    #            (RV_T['PROGRAM'] == 'bright') &
    (RV_T['SN_R'] > 10))


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


import matplotlib.gridspec as gridspec

gs = gridspec.GridSpec(nrows=1,
                       ncols=3,
                       width_ratios=[1, 1, 0.05],
                       wspace=0.05)
for cnt in range(2):
    T = [RV_T, SP_T][cnt]
    if cnt == 0:
        cur_sel = cur_sel0 & (RV_T['VSINI'] < 30)
    else:
        cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1')
    cur_sel = cur_sel & betw(T['TEFF'], 4500, 7000)
    #plt.subplot(1, 2, cnt + 1)
    plt.subplot(gs[cnt])
    im = plt.hist2d((T['FEH'][cur_sel]),
                    T['ALPHAFE'][cur_sel],
                    range=[[-5, 1], [-.5, 1.3]],
                    bins=[100, 100],
                    norm=maco.PowerNorm(gamma=.5, vmax=10000))

    plt.gci().set_rasterized(True)
    plt.xlabel('[Fe/H] ')

    #plt.title(['RVS', 'SP'][cnt])
    plt.text(-2, 1.21, ['RVS', 'SP'][cnt], color='white')
    if cnt == 0:
        plt.ylabel(r'[$\alpha$/Fe]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
ax = plt.subplot(gs[2])
plt.colorbar(im[-1], cax=ax)

# plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.01, top=.98, bottom=.13)

plt.savefig('plots/feh_alpha.pdf')
