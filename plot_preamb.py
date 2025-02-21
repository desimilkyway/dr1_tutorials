import matplotlib.pyplot as plt
"""

To get proper advantage of this package you need to make sure that you
don't rescale your figures when ingesting them in latex.
That means setting correct figure sizes when you make them.
Either by using something like
plt.figure(1, figsize=(3.,3))
before figure creation
or by using
plt.gcf().set_size_inches(3, 5)
after.

The with of half-column MNRAS figures is ~ 3.37 inches,
and for full width figures I assume 3.37x2
"""


def run(fontsize=7):
    """ This needs to be run before doing the plotting
    I.e.
    import plot_preamb as PP
    PP.run()
    """
    plt.rc(
        'font', **{
            'size': fontsize,
            'sans-serif': ['Helvetica', 'Helvetica Neue'],
            'family': 'sans-serif'
        })
    plt.rc('legend', **{'fontsize': fontsize})
    plt.rc("text.latex",
           preamble="\\usepackage{helvet}\\usepackage[T1]{fontenc}"
           "\\usepackage{sfmath}")
    plt.rc("text", usetex=True)
    plt.rc('xtick.minor', visible=True)
    plt.rc('ytick.minor', visible=True)
    plt.rc('xtick', direction='in')
    plt.rc('ytick', direction='in')
    plt.rc('xtick', top=True)
    plt.rc('ytick', right=True)
    plt.rc('ps', usedistiller='xpdf')
    plt.rc('savefig', **{'dpi': 300})
