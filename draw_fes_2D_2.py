#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",  # serif font via preamble
         ])
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

parser = argparse.ArgumentParser()
parser.add_argument("pmf", help = "specify the PMF file")
parser.add_argument("-o", "--output", default = "output.png", help = "specify the PNG output image file")
parser.add_argument("--xtitle", default = "X", help = "title along X axis")
parser.add_argument("--ytitle", default = "Y", help = "title along Y axis")
parser.add_argument("--levels", default = 25, type = int, help = "number of levels")
args = parser.parse_args()

def plotfes(pmffilename, pngfilename, xtitle, ytitle, level):
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    z = z - np.min(z)
    z = np.clip(z, 0, 16.0)
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    fig = plt.figure()
    cf = plt.contourf(xi, yi, zi, level, cmap='turbo')
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    #plt.title('Free Energy Surface')
    ax = plt.gca()
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_major_locator(plt.MaxNLocator(6, prune='lower'))
    clb = plt.colorbar()
    clb.ax.set_title("kcal/mol", fontsize=20)
    clb.ax.xaxis.get_major_formatter()._usetex = False
    clb.ax.yaxis.get_major_formatter()._usetex = False
    # clbticks = [float(i.get_position()[1]) for i in clb.ax.get_yticklabels()]
    # clbticksstr = ['{:.1f}'.format(i) for i in clbticks]
    # print(clbticksstr)
    # clb.ax.set_yticklabels(clbticksstr, fontsize = 20)
    plt.savefig(pngfilename, dpi=400, bbox_inches = 'tight', transparent=True)
    return

plotfes(args.pmf, args.output, xtitle = args.xtitle, ytitle = args.ytitle, level = args.levels)
#plotfes('trialanine_apath_cpp+1.czar.pmf', 'apath_2d.png', r'$\xi (s)$', r'$\xi (z)$')
#plotfes('trialanine_gpath_2d.czar.pmf', 'gpath_2d.png', r'$\xi (s)$', r'$\xi (z)$')
