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
    D['survey'] = [survey] * len(T)
    D['program'] = [program] * len(T)
    for k in [
            'TARGETID', 'VRAD', 'VRAD_ERR', 'RVS_WARN', 'FEH', 'LOGG', 'TEFF',
            'TARGET_RA', 'TARGET_DEC', 'RR_SPECTYPE', 'SN_R'
    ]:
        D[k] = T[k]
    D['MJD'] = T1['MJD']
    D['TILEID'] = T1['TILEID']
    D['SOURCE_ID'] = T2['SOURCE_ID']
    tabs.append(atpy.Table(D))
tabs = atpy.vstack(tabs)
tab_p = tabs.to_pandas()
conn1 = duckdb.connect(':memory:')
conn1.register('repeats0', tab_p)
conn1.execute('create table repeats as select * from repeats0')
sid, cnt, dmjd = sqlutil.get(
    'select source_id , count(*),max(mjd)-min(mjd) from repeats as r where source_id>0 group by SOURCE_ID',
    conn=conn1,
    driver='duckdb')
fig = plt.figure(figsize=(3.37, 3.37 * .75))
for i, limit in enumerate([1, 10, 100]):
    plt.hist(cnt[dmjd < limit],
             label='<%d day' % limit,
             range=[-0.5, 50.5],
             bins=51,
             histtype='step',
             alpha=0.5,
             color=['blue', 'green', 'red'][i])
plt.xlabel('N observations')
plt.legend()
plt.gca().set_yscale('log')
plt.tight_layout()
plt.savefig('plots/nrepeats.pdf')
