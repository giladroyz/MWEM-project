# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 22:33:33 2018

@author: user
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_experiment(file_name):
    
    epsilons = [0.0125, 0.025, 0.5, 0.1]
    files_titiles = {'Adult_age_hour.txt' : 'Adult: age x hours', 
                     'Adult_capitalloss.txt' : 'Adult: capital loss', 
                     'Transfusion_monetary.txt' : 'Transfusion: momentary', 
                     'Transfusion_recency_frequency.txt' : 'Transfusion: recency x frequency'}
    
    mwem = []
    svd = []
    
    with open(file_name) as f:
        for i in range(4):
            mwem.append(float(f.readline()))
            svd.append(float(f.readline()))
            f.readline()
            f.readline()
            
    #temp = mwem[2]
    #mwem[2] = mwem[3]
    #mwem[3] = temp
    
    #temp = svd[2]
    #svd[2] = svd[3]
    #svd[3] = temp
        
    plt.yscale('log')
    plt.title(files_titiles[file_name])
    plt.ylabel('mean square error')
    plt.xlabel('blue = mwem(T=10) , orange = SVD Lower Bound')
    plt.plot(epsilons, mwem)
    plt.plot(epsilons, svd)
    plt.show()
    
def main():
    
    files = ['Adult_age_hour.txt', 'Adult_capitalloss.txt', 'Transfusion_monetary.txt', 'Transfusion_recency_frequency.txt']
    for file in files:
        plot_experiment(file)
    
if __name__ == '__main__':
    main()

#epsilons = [0.0125, 0.025, 0.1, 0.5]
#mwem = [553313.9292, 209149.7025, 98145.71721, 94824.54253]
#svd = [4.60E+08, 1.15E+08, 7.18E+06, 287399.3726]