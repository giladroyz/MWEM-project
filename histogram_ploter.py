# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 22:33:33 2018

@author: user
"""

import matplotlib.pyplot as plt
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
from IPython.display import Image

def main():

	real = []

	with open('nltcs_histogram_real.txt', 'r') as f:
		for a in f:
			real.append(float(a))
		plt.hist(real, rwidth = 0.1)
		plt.show()

	with open('nltcs_histogram.txt', 'r') as f:
		for a in f:
			real.append(float(a))
		plt.hist(real, rwidth = 0.1)
		plt.show()

if __name__ == '__main__':
    main()