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

def plot_experiment(file_name):
    
    files_titiles = {'nltcs_result_new.txt' : 'NLTCS'}
    
    name_without_ext = file_name.split('.')[-2]

    epsilons = []
    mwem = {}
    svd = {}
    
    with open(file_name) as f:
        number_of_epsilons = int(f.readline())
        number_of_tests = int(f.readline())
        f.readline()

        epsilons_strings = f.readline().split(',')

        epsilons = [float(x) for x in epsilons_strings[:-1]]

        # read data of mwem
        for i in range(number_of_tests):
            test_epsilons = [float(x) for x in f.readline().split(',')[:-1]]
            for j, eps in enumerate(epsilons):
                mwem[(i,eps)] = test_epsilons[j]

    mwem_mean = [np.mean([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]

    mwem_std = [np.std([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]

    
    data_mwem = go.Scatter(
        x=epsilons,
        y=mwem_mean,
        name='mwem(T=10)',
        error_y=dict(
            type='data',
            array=mwem_std,
            visible=True
        )
    )

    layout = go.Layout(
    title=files_titiles[file_name],
    xaxis=dict(
        title='epsilon',
        autorange=True
        #fixedrange=True
    ),
    yaxis=dict(
        title='relative entropy',
        #type='log',
        autorange=True,
        exponentformat='E'
    )
    )

    data = [data_mwem]
    fig = go.Figure(data=data, layout=layout)

    #py.image.save_as(fig, filename=name_without_ext + '.png')

    plotly.offline.plot(fig, filename=name_without_ext)
    #Image(name_without_ext + '.png')
    
def main():
    
    files = ['nltcs_result_new.txt']
    for file in files:
        plot_experiment(file)
    
if __name__ == '__main__':
    main()

#epsilons = [0.0125, 0.025, 0.1, 0.5]
#mwem = [553313.9292, 209149.7025, 98145.71721, 94824.54253]
#svd = [4.60E+08, 1.15E+08, 7.18E+06, 287399.3726]