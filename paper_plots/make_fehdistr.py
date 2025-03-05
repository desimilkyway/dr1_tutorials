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
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .75))
cnt = 0
cnt = 0
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


cur_sel = cur_sel0 & betw(RV_T['TEFF'], 4500, 7000) & (RV_T['VSINI'] < 30)
plt.hist(RV_T['FEH'][cur_sel], range=[-5, .5], label='RVS', alpha=.5, bins=100)
plt.xlabel('[Fe/H]')
cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1') & betw(
    SP_T['TEFF'], 4500, 7000)
plt.hist(SP_T['FEH'][cur_sel], range=[-5, .5], label='SP', alpha=.5, bins=100)
plt.gca().set_yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plots/feh_distr.pdf')
