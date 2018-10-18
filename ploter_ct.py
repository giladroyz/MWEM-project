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

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def plot_experiment(file_name):
    
    files_titiles = {'optimization.txt' : 'Adult: age x hours - Optimization'}
    
    name_without_ext = file_name.split('.')[-2]

    epsilons = [float(x)/20 for x in range(1,20)]
    print(epsilons)
    mwem = {}
    
    with open(file_name) as f:
        line = f.readline().split(',')
        number_of_queries = int(line[0])
        number_iterations = int(line[1])
        number_of_tests = int(line[0])

        for i in range(len(epsilons)**2):
            line = f.readline()
            e1,e2 = tuple([float(x) for x in line.split(',')])
            line = f.readline()
            mwem[(e1,e2)] = float(line.rstrip())

    #mwem_mean = [np.mean([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]
    #svd_mean = [np.mean([svd[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]

    #mwem_std = [np.std([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]
    #svd_std = [np.std([svd[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]


    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.array(epsilons)
    Y = np.array(epsilons)
    X, Y = np.meshgrid(X, Y)
    Z = np.array([[mwem[(epsilons[i],epsilons[j])] for j in range(len(epsilons))] for i in range(len(epsilons))])
    #Z /= Z.max()
    Z = np.log2(Z+2)
    Z /= Z.max()

    print(Z)

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

    
def main():
    
    files = ['optimization.txt']
    for file in files:
        plot_experiment(file)
    
if __name__ == '__main__':
    main()

#epsilons = [0.0125, 0.025, 0.1, 0.5]
#mwem = [553313.9292, 209149.7025, 98145.71721, 94824.54253]
#svd = [4.60E+08, 1.15E+08, 7.18E+06, 287399.3726]