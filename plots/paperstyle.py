#!/usr/bin/env python3

'''
This file contains the style settings for the paper plots
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Set the font
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

area = dict(color='green', marker='o', markerfacecolor='green', markeredgecolor='green', alpha=0.7, label='Area Method')
dnn = dict(color='blue',  marker='s', markerfacecolor='blue', markeredgecolor='blue', alpha=0.7, label='Deep Neural Network')
mult = dict(color='red',  marker='P', markerfacecolor='red', markeredgecolor='red', alpha=0.7, label='Multiplicity Method')
snn = dict(color='black',  marker='^', markerfacecolor='black', markeredgecolor='black', alpha=0.7, label='Shallow Neural Network')

