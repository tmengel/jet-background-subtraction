#!/usr/bin/env python3

'''
This file contains the style settings for the paper plots
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# # Set the font
style = {'axes.labelsize': 12,
            'axes.linewidth' : 1.5,
            'font.size': 12,
            'font.family': 'Times New Roman',
            'mathtext.fontset': 'stix',
            'legend.fontsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'text.usetex': False,
            'lines.linewidth': 1,
            'lines.linestyle': ' ',
            'lines.markersize' : 6,
            'lines.markeredgewidth' : 1,
            'xtick.major.size' : 5,
            'xtick.minor.size' : 3,
            'xtick.major.width' : 2,
            'xtick.minor.width' : 1,
            'xtick.direction' : 'in',
            'ytick.major.size' : 5,
            'ytick.minor.size' : 3,
            'ytick.major.width' : 2,
            'ytick.minor.width' : 1,
            'ytick.direction' : 'in',
            'xtick.minor.visible' : True,
            'ytick.minor.visible' : True,
            'savefig.transparent': True,
            'errorbar.capsize': 1.5,
            }
# style = {
#     'axes.facecolor':'white',
#     'axes.edgecolor': 'black',
#     'figure.facecolor' : 'white',
#     'figure.edgecolor' : 'white',
#     'figure.frameon' : False,
#     'figure.subplot.left' : 0.1,
#     'figure.subplot.right' : 0.95,
#     'figure.subplot.bottom' : 0.1,
#     'figure.subplot.top' : 0.95,
    
#     'axes.labelpad' : 1.4,
#     'axes.labelcolor' : 'black',
#     'axes.labelsize': 12,
#     'axes.linewidth' : 1.0,
#     'axes.titlesize' : 'medium',
#     'xaxis.labellocation': 'right',
#     'yaxis.labellocation': 'top',
#     'axes.titlepad': 0.0,
    
#     # Helvetica
#     'text.color': 'black',
#     'font.family': 'Times New Roman',
#     'mathtext.fontset': 'stix',
#     'font.size': 12,
#     'legend.fontsize': 12,
#     'xtick.labelsize': 12,
#     'ytick.labelsize': 12,
#     'text.usetex': False,
#     'mathtext.it' : 'sans:italic',
#     'mathtext.default' : 'regular',

#     # ticks
#     'xtick.major.size' : 7,
#     'xtick.minor.size' : 5,
#     'xtick.major.width' : 1,
#     'xtick.minor.width' : 1,
#     'xtick.direction' : 'in',
#     'ytick.major.size' : 7,
#     'ytick.minor.size' : 5,
#     'ytick.major.width' : 1,
#     'ytick.minor.width' : 1,
#     'ytick.direction' : 'in',
#     'ytick.right' : True,
#     'xtick.top' : True,
#     'xtick.minor.visible' : True,
#     'ytick.minor.visible' : True,
#     'xtick.major.pad': 2,
#     'ytick.major.pad': 2,

#     'lines.linewidth': 1.0,
#     'lines.linestyle': '-',
#     'lines.markersize' : 5,
#     'lines.marker' : 'o',
#     'lines.markeredgewidth' : 1,

#     'savefig.transparent': False,
#     'errorbar.capsize': 0,

#     # legend
#     'legend.frameon': False,
#     'legend.loc': 'best',
#     'legend.numpoints': 1,
#     'legend.scatterpoints': 1,
#     'legend.fontsize': 12,
#     'legend.handlelength': 1.5,
#     'legend.handletextpad': 0.5,
#     'legend.labelspacing': 0.5,
#     'legend.borderpad': 0.4,
#     'legend.borderaxespad': 0.5,
#     'legend.columnspacing': 0.7,
#     'legend.markerscale': 1.0,
#     'legend.shadow': False

#     }


nn_style = {'axes.labelsize': 12,
            'axes.linewidth' : 1.5,
            'font.size': 12,
            'font.family': 'Times New Roman',
            'mathtext.fontset': 'stix',
            'legend.fontsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'text.usetex': False,
            'lines.linewidth': 1,
            'lines.linestyle': ' ',
            'lines.markersize' : 6,
            'lines.markeredgewidth' : 1,
            'xtick.major.size' : 5,
            'xtick.minor.size' : 0,
            'xtick.major.width' : 2,
            'xtick.minor.width' : 0,
            'xtick.direction' : 'in',
            'ytick.major.size' : 5,
            'ytick.minor.size' : 0,
            'ytick.major.width' : 2,
            'ytick.minor.width' : 0,
            'ytick.direction' : 'in',
            'xtick.minor.visible' : True,
            'ytick.minor.visible' : True,
            'savefig.transparent': True,
            'errorbar.capsize': 1.5,
            }


plt.rcParams.update(style)

def GetFigSize(fig_width_pt=510.0, denominator=3.0):
    inches_per_pt = 1.0/72.27 
    golden_mean = (np.sqrt(5)-1.0)/denominator
    fig_width = fig_width_pt*inches_per_pt
    fig_height = fig_width*golden_mean
    return [fig_width, fig_height]

area = dict(color='green', marker='o', markerfacecolor='green', markeredgecolor='green', alpha=1.0, label='Area Method')
dnn = dict(color='blue',  marker='s', markerfacecolor='blue', markeredgecolor='blue', alpha=1.0, label='Neural Network')
mult = dict(color='red',  marker='P', markerfacecolor='red', markeredgecolor='red', alpha=1.0, label='Multiplicity Method')
snn = dict(color='black',  marker='^', markerfacecolor='black', markeredgecolor='black', alpha=0.7, label='Shallow Neural Network')

