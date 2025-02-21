import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as maco
import plot_preamb as pp

pp.run()

RV_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
FM_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'FIBERMAP',
                         mask_invalid=False)
G_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                        'GAIA',
                        mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
fig = plt.figure(figsize=(3.37 * 2, 3.37 * .75))
cnt = 0
cnt = 0
for survey, program in [('sv3', 'bright'), ('main', 'bright'),
                        ('main', 'dark')]:
    cur_sel = main_sel & (RV_T['SURVEY'] == survey) & (RV_T['PROGRAM']
                                                       == program)
    plt.subplot(1, 3, cnt + 1)
    plt.hist2d(-2.5 *
               np.log10(FM_T['FLUX_G'][cur_sel] / FM_T['FLUX_R'][cur_sel]),
               22.5 - 2.5 * np.log10(FM_T['FLUX_R'][cur_sel]),
               bins=[100, 100],
               range=[[-0.49, 2], [15.5, 21.5]])
    plt.gca().set_rasterized(True)
    cnt += 1
    plt.title(f'survey, program: {survey},{program}')
    plt.xlabel('g-r [mag]')
    plt.ylim(21.5, 15.5)
    if cnt == 1:
        plt.ylabel('r [mag]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.01)

plt.savefig('plots/cmd.pdf')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * 1))

cur_sel = main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'backup')
plt.hist2d(G_T['BP_RP'][cur_sel],
           G_T['PHOT_G_MEAN_MAG'][cur_sel],
           bins=[100, 100],
           range=[[0.01, 3], [11, 20]])
plt.gca().set_rasterized(True)
plt.ylim(20, 11)
plt.title('survey, program: main,backup')
plt.xlabel('BP-RP [mag]')
plt.ylabel('G [mag]')
plt.tight_layout()

plt.savefig('plots/cmd_backup.pdf')
