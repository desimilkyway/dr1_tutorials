import glob
import astropy.table as atpy
import duckdb
import sqlutilpy as sqlutil
import scipy.stats
import matplotlib.pyplot as plt
from idlplotInd import plot, oplot
import numpy as np
import plot_preamb as pp

pp.run()
fs = glob.glob('../../rv_variability/rvtabs_iron/*exp*fits')
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
fname = '../data/mwsall-pix-iron.fits'

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
bins = 20
floor_dict = {'dark': 1.2, 'bright': .7, 'backup': 2}
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .9))

for i, prog in enumerate(['dark', 'bright', 'backup']):
    sel1 = (PAIRS['program'] == prog) & (PAIRS['survey'] == survey) & (
        PAIRS['rvs_warn1'] == 0) & (PAIRS['rvs_warn2'] == 0)
    floor = floor_dict[prog]
    erange = [-1, 1.5]
    sel2 = sel1 & (np.abs(PAIRS['mjd1'] - PAIRS['mjd2']) > 1)
    SS = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                      delt[sel1],
                                      funcer,
                                      range=erange,
                                      bins=bins)
    SC = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                      delt[sel1],
                                      'count',
                                      range=erange,
                                      bins=bins)
    SS2 = scipy.stats.binned_statistic(np.log10(comb_err[sel2]),
                                       delt[sel2],
                                       funcer,
                                       range=erange,
                                       bins=bins)
    SC2 = scipy.stats.binned_statistic(np.log10(comb_err[sel2]),
                                       delt[sel2],
                                       'count',
                                       range=erange,
                                       bins=bins)

    plt.subplot(3, 1, i + 1)
    plot(10**(SS.bin_edges[:-1] + .5 * np.diff(SS.bin_edges)), (SS.statistic),
         ps=3,
         ylog=True,
         xlog=True,
         yr=[.5, 30],
         xr=[.11, 30],
         noerase=True,
         ind=SC.statistic > 100)
    oplot(10**(SS2.bin_edges[:-1] + .5 * np.diff(SS2.bin_edges)),
          (SS2.statistic),
          ps=3,
          color='grey',
          ind=SC2.statistic > 100)
    if i == 1:
        plt.ylabel(r'$\frac{1}{\sqrt{2}}$ StdDev($V_1-V_2$) [km/s]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    if i == 2:
        plt.xlabel(r'$\sqrt{\frac{\sigma_1^2+\sigma_2^2}{2}}$ [km/s]')
    plt.text(.2, 10, f'Survey, program: {survey},{prog}')
    oplot(10**xgrid,
          np.sqrt(10**(2 * xgrid) + floor**2),
          label='floor %s km/s' % floor)
    plt.legend()
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('plots/repeats_accuracy.pdf')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .8))
for i, prog in enumerate(['dark', 'bright', 'backup']):
    sel1 = (PAIRS['program'] == prog) & (PAIRS['survey'] == survey) & (
        PAIRS['rvs_warn1'] == 0) & (PAIRS['rvs_warn2'] == 0)
    plt.subplot(3, 1, i + 1)
    floor = floor_dict[prog]
    xdelt = delt / (comb_err**2 + floor**2)**.5
    xr = [-6, 6]
    nbins = 100
    plt.hist(xdelt[sel1], range=xr, bins=nbins, histtype='step', color='black')

    xgrid = np.linspace(-6, 6, 1000)
    plt.plot(xgrid,
             scipy.stats.norm(0, 1).pdf(xgrid) * sel1.sum() * (xr[1] - xr[0]) /
             nbins,
             color='red')
    plt.text(.6,
             .7,
             f'Survey, program:\n {survey},{prog}',
             transform=plt.gca().transAxes)

    # if i == 0:
    plt.ylabel(r'N/bin')
    if i < 2:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
    if i == 2:
        plt.xlabel(r'$\delta_{RV}/\sigma_{RV,calib}$')
    # plt.legend()
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('plots/repeats_delta.pdf')
