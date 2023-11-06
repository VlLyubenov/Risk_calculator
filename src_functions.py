# -*- coding: utf-8 -*-

# -- Sheet --

from typing import Union, Any

from numpy import ndarray, dtype, floating
from numpy._typing import _64Bit


def var_calc(values, percent):


    def partition(array,low,high):
    
        pivot = array[high]
        i = low - 1

        for j in range (low, high):
            if array[j] <= pivot:
                i = i + 1
                (array[i], array[j]) = (array[j], array[i])
        (array[i + 1], array[high]) = (array[high], array[i + 1])
        return  i + 1

    def quicksort(array, low, high):
        if low < high:
            pi = partition(array, low, high)
            quicksort(array, low, pi - 1)
            quicksort(array, pi + 1, high)
    
    data = values
    size = len(data)
    quicksort(data, 0, size - 1)

   
    percentile = size * percent
    percentile = round(percentile)
    if percentile < 1:
        percentile = 1
    return data[percentile]


def plot_var(values,percent):

    import matplotlib.pyplot as plt

    def var_calc(var_values, var_percent):


        def partition(array,low,high):
        
            pivot = array[high]
            i = low - 1

            for j in range (low, high):
                if array[j] <= pivot:
                    i = i + 1
                    (array[i], array[j]) = (array[j], array[i])
            (array[i + 1], array[high]) = (array[high], array[i + 1])
            return  i + 1

        def quicksort(array, low, high):
            if low < high:
                pi = partition(array, low, high)
                quicksort(array, low, pi - 1)
                quicksort(array, pi + 1, high)
    
        data = var_values
        size = len(data)
        quicksort(data, 0, size - 1)

    
        percentile = size * var_percent
        percentile = round(percentile)
        if percentile < 1:
            percentile = 1
        return data[percentile]
    
   
    var_value = var_calc(values, percent)

    #VaR plot
    plt.hist(values,len(values))
    plt.axvline(var_value, color='red', linestyle='dashed', linewidth=2)
    plt.title("VaR")

def cvar_calc(values, percent):


    def var_calc(var_values, var_percent):


        def partition(array,low,high):
        
            pivot = array[high]
            i = low - 1

            for j in range (low, high):
                if array[j] <= pivot:
                    i = i + 1
                    (array[i], array[j]) = (array[j], array[i])
            (array[i + 1], array[high]) = (array[high], array[i + 1])
            return  i + 1

        def quicksort(array, low, high):
            if low < high:
                pi = partition(array, low, high)
                quicksort(array, low, pi - 1)
                quicksort(array, pi + 1, high)
        
        data = var_values
        size = len(data)
        quicksort(data, 0, size - 1)

    
        percentile = size * var_percent
        percentile = round(percentile)
        if percentile < 1:
            percentile = 1
        return data[percentile]
    
    v = var_calc(values, percent)

    cvar_list = []
    for j in values:
        if j < v:
            cvar_list.append(j)
    
    cvar_sum = 0
    for j in cvar_list:
        cvar_sum = j + cvar_sum

    cvar = cvar_sum / len(cvar_list)
    return cvar



def value_at_risk_analysis(values, percent):
    print("VaR: {}".format(var_calc(values, percent))) 
    print("CVaR: {}".format(cvar_calc(values, percent)))
    print("Standard Deviation: {}".format(get_sd(values)))
    import matplotlib.pyplot as plt
    return plot_var(values,percent)

#Statistics

def get_sd(values):
    from statistics import stdev
    return stdev(values)

def get_mean(values):
    from statistics import mean
    return mean(values)

def get_median(values):
    from statistics import median
    return median(values)

def get_variance(values):
    from statistics import variance
    return variance(values)

def get_skewness(values):
    from scipy.stats import skew
    return skew(values)

def get_kurtosis(values):
    from scipy.stats import kurtosis
    return kurtosis(values)

# min(values)
# max(values)

def get_lower_quantile(values):
    from numpy import quantile
    return quantile(values, 0.25)

def get_heigher_quantile(values):
    from numpy import quantile
    return quantile(values, 0.75)



# normal distribution var

def norm_var(values, percent):

    from scipy.stats import norm
    cutoff = norm.ppf(percent, get_mean(values), get_sd(values))

    return cutoff

# monte carlo var

def monte_carlo(mean, sd, n_simulations):

    from numpy.random import normal

    array_rands = normal(mean, sd, n_simulations)

    return array_rands


def monte_carlo_var(var_percent, mean, sd, n_simulations):
    
    values = monte_carlo(mean, sd, n_simulations)
    return var_calc(values, var_percent)

def monte_carlo_cvar(var_percent, mean, sd, n_simulations):
    
    values = monte_carlo(mean, sd, n_simulations)
    return cvar_calc(values, var_percent)

    

def sd_bound_lower(values):
    sd = get_sd(values)

    sd_bound_array = []
    for i in values:
        sd_bound = i - sd
        sd_bound_array.append(sd_bound)

    return sd_bound_array

def sd_bound_upper(values):
    sd = get_sd(values)

    sd_bound_array = []
    for i in values:
        sd_bound = i + sd
        sd_bound_array.append(sd_bound)

    return sd_bound_array

# GARCH VaR

# Require arch packege!!!

#Doesnt work!

def get_garch_var(values, percent):
    
    import pandas as pd
    from arch import arch_model
    from arch.__future__ import reindexing
    from numpy.random import normal

    garch= arch_model(values , vol="Garch", p=1, o=0, q=1, dist="Normal")
    result= garch.fit(update_freq=1)
    forecasts = result.forecast()
    cond_mean = forecasts.mean
    cond_var = forecasts.variance
    m = float(cond_mean.values) - 100
    v = float(cond_var.values)

    distribution = normal(m, v , len(values))
    var_calc(distribution,percent)



# -- Sheet 2 --

import pandas as pd

def get_SPreturns(csv_file):
   
    SP = pd.read_csv(csv_file)

    SPprice = []
    SPprice = SP["Price"].str.replace(',','')
    SPprice = SPprice.apply(float)
    returns = SPprice.pct_change()
    returns = returns.values
    returns = list(returns)
    returns.pop(0)
    return returns


returns = list(get_SPreturns("/data/notebook_files/S&P 500 Historical Data (1).csv"))

value_at_risk_analysis(returns, 0.01)

get_sd(returns)

norm_var(returns, 0.05)

sd_bound_lower(returns)

import pandas as pd
from arch import arch_model
from arch.__future__ import reindexing
from numpy.random import normal


full_percent_returns = []
for i in returns:
    r = (i +1)*100
    full_percent_returns.append(r)

garch= arch_model(full_percent_returns , vol="Garch", p=1, o=0, q=1, dist="Normal")
result= garch.fit(update_freq=1)
forecasts = result.forecast()
cond_mean = forecasts.mean
cond_var = forecasts.variance
m = float(cond_mean.values) - 100
v = float(cond_var.values)

distribution = normal(m, v , len(returns))
var_calc(distribution,0.05)

print(garch_array)



# -- multivatate_cases --



get_input()

