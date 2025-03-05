from desitarget import targetmask as tm
import astropy.table as atpy
import re

T0 = atpy.Table().read('../data/mwsall-pix-iron.fits', mask_invalid=False)
TF0 = atpy.Table().read('../data/mwsall-pix-iron.fits',
                        'FIBERMAP',
                        mask_invalid=False)


def fixer(x):
    return x.replace('_', '\\_')


for survey in ['main']:
    mm = tm.mws_mask
    for program in ['dark', 'bright', 'backup']:
        sel = (T0['SURVEY'] == survey) & (T0['PROGRAM'] == program)

        # print(T0['SURVEY'], T0['PROGRAM'])
        T = T0[sel]
        TF = TF0[sel]
        if sel.sum() == 0:
            continue

        print('\\hline', survey, program, ' & \\ \\hline')
        SS = survey.upper()
        for n in mm.names():
            if re.match('.*NORTH.*', n) is not None:
                continue
            if re.match('.*SOUTH.*', n) is not None:
                continue
            xind = (TF['MWS_TARGET'] & mm[n].mask) > 0
            if xind.sum() == 0:

                continue
            print(fixer(n), '& ', format(xind.sum(), ','), '\\\\')
