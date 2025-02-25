import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp

pp.run()

RV_T = atpy.Table().read('data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
fig = plt.figure(figsize=(3.37 * 2, 3.37 * 1.3))
cnt = 0
for survey, program in [('sv3', 'bright'), ('main', 'bright'),
                        ('main', 'dark'), ('main', 'backup')]:
    cur_sel = main_sel & (RV_T['SURVEY'] == survey) & (RV_T['PROGRAM']
                                                       == program)
    plt.subplot(2, 2, cnt + 1)
    plt.hist2d(
        RV_T['TARGET_RA'][cur_sel],
        RV_T['TARGET_DEC'][cur_sel],
        bins=[360, 120],
        range=[[0, 360], [-30, 90]],
        weights=1. / np.cos(np.deg2rad(RV_T['TARGET_DEC'][cur_sel])),
    )
    plt.gca().set_rasterized(True)
    if cnt < 2:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
    else:
        plt.xlabel(r'$\alpha$ [deg]')
    cnt += 1
    #plt.title(f'survey, program: {survey},{program}')
    plt.text(300, 85, f'survey, program: {survey},{program}', color='white')
    plt.xlim(360, 0)
    plt.ylabel(r'$\delta$ [deg]')
    plt.colorbar()
plt.tight_layout()
plt.subplots_adjust(wspace=0.12, hspace=0.01)

plt.savefig('plots/density.pdf')
