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
    
    files_titiles = {'Adult_age_hour.txt' : 'Adult: age x hours', 
                     'Adult_capitalloss.txt' : 'Adult: capital loss', 
                     'Transfusion_monetary.txt' : 'Transfusion: momentary', 
                     'Transfusion_recency_frequency.txt' : 'Transfusion: recency x frequency'}
    
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

        f.readline()


        # read data of svd
        for i in range(number_of_tests):
            test_epsilons = [float(x) for x in f.readline().split(',')[:-1]]
            for j, eps in enumerate(epsilons):
                svd[(i,eps)] = test_epsilons[j]

        f.readline()

    mwem_mean = [np.mean([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]
    svd_mean = [np.mean([svd[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]

    mwem_std = [np.std([mwem[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]
    svd_std = [np.std([svd[(i,eps)] for i in range(number_of_tests)]) for eps in epsilons]

    
    data_mwem = go.Scatter(
        x=epsilons,
        y=mwem_mean,
        name='mwem(T=10)',
        error_y=dict(
        	name='std of mwem',
            type='data',
            array=mwem_std,
            visible=True
        )
    )

    data_svd = go.Scatter(
        x=epsilons,
        y=svd_mean,
        name='SVD Lower Bound',
        error_y=dict(
        	name='std of svd',
            type='data',
            array=svd_std,
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
    	title='mean square error',
        type='log',
        autorange=True,
        exponentformat='E'
    )
    )

    data = [data_mwem, data_svd]
    fig = go.Figure(data=data, layout=layout)

    #py.image.save_as(fig, filename=name_without_ext + '.png')

    plotly.offline.plot(fig, filename=name_without_ext)
    #Image(name_without_ext + '.png')
    
def main():
    
    files = ['Adult_age_hour.txt', 'Adult_capitalloss.txt', 'Transfusion_monetary.txt', 'Transfusion_recency_frequency.txt']
    for file in files:
        plot_experiment(file)
    
if __name__ == '__main__':
    main()

#epsilons = [0.0125, 0.025, 0.1, 0.5]
#mwem = [553313.9292, 209149.7025, 98145.71721, 94824.54253]
#svd = [4.60E+08, 1.15E+08, 7.18E+06, 287399.3726]