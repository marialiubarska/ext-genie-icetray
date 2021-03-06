#!/usr/bin/env python 

from icecube import dataio
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from optparse import OptionParser
from os import path

# input
usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-i", "--infile", default="",
                  dest="INFILE", help="Input .i3 file")
parser.add_option("-o", "--outpdf", default="test",
                  dest="OUTPDF", help="Output .pdf figure name")


(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

# reading data from file 
with dataio.I3File(options.INFILE) as fl:
    
    fr = fl.pop_daq()
    wd = fr['I3MCWeightDict']

    wdict = {}
    for k in wd.keys():
        wdict[k] = []
    print wdict.keys()

    for frame in fl:
        for k in wdict.keys():
            wdict[k].append( frame['I3MCWeightDict'][k] )

# plotting
def plot_val(ax, v, title):
    
    ax.set_title(title)
    h = ax.hist(v, 30, histtype='step', linewidth=2.)

    print k+':'
    print 'hist:', h[0]
    print 'bins:', h[1]

# figure out number of subplots
l = len(wdict)
n_sp1 = int(np.sqrt(l))
n_sp2 = int(l/n_sp1) + 1

# plot everything 
fig, axs = plt.subplots(n_sp1, n_sp2, figsize=(25,15))
fig.suptitle(options.INFILE)

for i in range(l):
    
    ax = axs.flat[i]
    k = wdict.keys()[i]
    plot_val(ax=ax, v=wdict[k], title=k)

plt.savefig(options.OUTPDF+'.pdf', dpi='figure')


