import glob
import astropy.table as atpy
# import crossmatcher
import duckdb
import sqlutilpy as sqlutil
import scipy.stats
import matplotlib.pyplot as plt
from idlplotInd import plot, oplot
import numpy as np

fs = glob.glob('../rv_variability/rvtabs_iron/*exp*fits')
tabs = []
for f in fs:
    T = atpy.Table().read(f, 'RVTAB', mask_invalid=False)
    T1 = atpy.Table().read(f, 'FIBERMAP', mask_invalid=False)
    T2 = atpy.Table().read(f, 'GAIA', mask_invalid=False)
    T.convert_bytestring_to_unicode()
    T1.convert_bytestring_to_unicode()
    T2.convert_bytestring_to_unicode()
    D = {}
    XF = f.replace('.fits', '').split('-')
    survey = XF[-2]
    program = XF[-1]
    # D['G'] = T2['PHOT_G_MEAN_MAG']
    # D['BP'] = T2['PHOT_BP_MEAN_MAG']
    # D['RP'] = T2['PHOT_RP_MEAN_MAG']
    # D['PARALLAX'] = T2['PARALLAX']
    # D['PARALLAX_ERROR'] = T2['PARALLAX_ERROR']
    # D['RUWE'] = T2['RUWE'].filled(0)
    # D['source_id'] = T2['SOURCE_ID'].filled(-1)
    D['survey'] = [survey] * len(T)
    D['program'] = [program] * len(T)
    for k in [
            'TARGETID', 'VRAD', 'VRAD_ERR', 'RVS_WARN', 'FEH', 'LOGG', 'TEFF',
            'TARGET_RA', 'TARGET_DEC', 'RR_SPECTYPE', 'SN_R'
    ]:
        D[k] = T[k]
    D['MJD'] = T1['MJD']
    D['TILEID'] = T1['TILEID']
    tabs.append(atpy.Table(D))
tabs = atpy.vstack(tabs)
tab_p = tabs.to_pandas()
conn1 = duckdb.connect(':memory:')
fname = 'data/mwsall-pix-iron.fits'

Tstack = atpy.Table().read(fname, mask_invalid=False)
Tstack2 = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
Tstack3 = atpy.Table().read(fname, 'GAIA', mask_invalid=False)
Tstack.convert_bytestring_to_unicode()
Tstack2.convert_bytestring_to_unicode()
Tstack3.convert_bytestring_to_unicode()
for k in [
        'SCND_TARGET', 'MWS_TARGET', 'SV1_MWS_TARGET', 'SV2_MWS_TARGET',
        'SV3_MWS_TARGET', 'SCND_TARGET', 'SV1_SCND_TARGET', 'SV1_SCND_TARGET',
        'SV1_SCND_TARGET', 'FLUX_G', 'FLUX_R'
]:
    try:
        Tstack[k] = Tstack2[k].filled(0)
    except AttributeError:
        Tstack[k] = Tstack2[k]
Tstack['source_id'] = Tstack3['SOURCE_ID'].filled(-1)
for k in [
        'PMRA', 'PMDEC', 'PARALLAX', 'PARALLAX_ERROR', 'RUWE',
        'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'DEC_ERROR',
        'ASTROMETRIC_N_GOOD_OBS_AL'
]:
    Tstack[k] = Tstack3[k]
Tstack = Tstack.to_pandas()

conn1.register('repeats0', tab_p)
conn1.register('stack0', Tstack)

conn1.execute('create table repeats as select * from repeats0')
conn1.execute('create table stack as select * from stack0')
conn1.execute('create index qq1 on stack(source_id)')
conn1.execute('create index qq2 on stack(targetid)')

PAIRS = sqlutil.get('''
select r1.mjd as mjd1,
r2.mjd as mjd2,
r1.vrad as vrad1,
r2.vrad as vrad2,
r1.vrad_err as vrad_err1,
r2.vrad_err as vrad_err2,
r1.survey as survey,
r1.program as program,
r1.rvs_warn as rvs_warn1,
r2.rvs_warn as rvs_warn2,
from
repeats as r1,
repeats as r2
where r1.targetid = r2.targetid
and r1.survey=r2.survey
and r1.program=r2.program
and r1.mjd<r2.mjd
and r1.rr_spectype='STAR'
''',
                    conn=conn1,
                    driver='duckdb',
                    asDict=True)

comb_err, delt = ((PAIRS['vrad_err1']**2 + PAIRS['vrad_err2']**2) /
                  2)**.5, (PAIRS['vrad1'] - PAIRS['vrad2']) / np.sqrt(2)

xgrid = np.linspace(-1, 1.5, 100)


def funcer(x):
    return (scipy.stats.scoreatpercentile(x, 84) -
            scipy.stats.scoreatpercentile(x, 16)) / 2.


plt.clf()
survey = 'main'
for i, prog in enumerate(['dark', 'bright', 'backup']):
    sel1 = (PAIRS['program'] == prog) & (PAIRS['survey'] == survey) & (
        PAIRS['rvs_warn1'] == 0) & (PAIRS['rvs_warn2'] == 0)
    SS = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                      delt[sel1],
                                      funcer,
                                      range=[-1, 1.5],
                                      bins=30)
    SC = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                      delt[sel1],
                                      'count',
                                      range=[-1, 1.5],
                                      bins=30)

    plt.subplot(1, 3, i + 1)
    plot(SS.bin_edges[:-1] + .5 * np.diff(SS.bin_edges), (SS.statistic),
         ps=3,
         ylog=True,
         yr=[.5, 30],
         xr=[-1, 1.5],
         xtitle=r'$\log_{10}(\sqrt{(\sigma_1^2+\sigma_2^2)/2}/(1 km/s))$',
         ytitle=r'$StdDev(V_1-V_2)/\sqrt{2} [km/s]$',
         noerase=True,
         title=f'{survey}/{prog}',
         ind=SC.statistic > 100)
    floor = {'dark': 1.3, 'bright': .7, 'backup': 2}[prog]
    oplot(xgrid,
          np.sqrt(10**(2 * xgrid) + floor**2),
          label='floor %s km/s' % floor)
    plt.legend()
plt.gcf().set_size_inches(10, 5)
plt.tight_layout()
plt.savefig('repeats_main.pdf')
