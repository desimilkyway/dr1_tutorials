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
fig = plt.figure(figsize=(3.37 * 2, 3.37 * 1.))
cnt = 0
bitmasks = [('MWS_BROAD', 0), ('MWS_WD', 1), ('MWS_NEARBY', 2), ('MWS_BHB', 6),
            ('MWS_MAIN_BLUE', 8), ('MWS_MAIN_RED', 11)]

survey, program = ('main', 'bright')

for curt, bit in bitmasks:
    objtype_sel = (FM_T['MWS_TARGET'] & (2**bit)) > 0
    cur_sel = objtype_sel & main_sel & (RV_T['SURVEY'] == survey) & (
        RV_T['PROGRAM'] == program)
    plt.subplot(2, 3, cnt + 1)
    plt.hist2d(-2.5 *
               np.log10(FM_T['FLUX_G'][cur_sel] / FM_T['FLUX_R'][cur_sel]),
               22.5 - 2.5 * np.log10(FM_T['FLUX_R'][cur_sel]),
               bins=[100, 100],
               range=[[-0.49, 2], [15.51, 20.5]])
    plt.gca().set_rasterized(True)
    if cnt % 3 == 0:
        plt.ylabel('r [mag]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    cnt += 1
    # plt.title(f'survey, program: {survey},{program}')
    plt.ylim(20.5, 15.51)
    plt.text(0, 15.9, f'{curt}, bitmask {2**bit}', color='white')
    plt.xlabel('g-r [mag]')
plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.01)

plt.savefig('plots/targ_classes.pdf')
