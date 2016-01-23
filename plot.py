#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# Just a figure and one subplot
f, axarr = plt.subplots(1,2, sharey=True)
axarr[0].plot(x, y)
axarr[0].set_title('Simple plot')
axarr[0].set_xlim((0,7))
axarr[0].set_xlabel('x')
axarr[0].set_ylabel(r'sin $x^{2}$')
axarr[0].set(adjustable='box-forced',aspect='equal')
axarr[1].scatter(x, y)
axarr[1].set_title('Scatter plot')
axarr[1].set_xlim((0,7))
axarr[1].set_ylim((-1,1))
axarr[1].set_xlabel('x')
# f.subplots_adjust(wspace=0)

# Show the plot
plt.show()
