import numpy as np
import sympy
from sympy import *
import matplotlib.pylab as plot
import matplotlib.pyplot as plt

no_range=100
def diff_func(equation,variable,parms):
    """finds partial derivative divided by equation (R0) for particular variable"""
    dR0dv=diff(equation,variable)
    dR0dv_proportion=dR0dv/equation
    for i in dR0dv_proportion.free_symbols: 
        dR0dv_proportion=dR0dv_proportion.subs(i,parms[str(i)])
    return(dR0dv_proportion)

def loop_diff(parms,equation,abs_req):
    """loop through parameters in equation(R0), call diff_func and store importance of each parameter"""
    importance_dictionary={}
    for k in equation.free_symbols:
        if abs_req==True:
            importance_dictionary[str(k)]=abs(diff_func(equation,k,parms))
        else:
            importance_dictionary[str(k)]=diff_func(equation,k,parms)
    return(importance_dictionary)

def plot_parm_imp(parms,equation,abs_req):
    """plot for one set of variables"""
    rank=loop_diff(parms,equation,abs_req)
 
    variable= [i for i in rank]
    value=[rank[str(i)] for i in rank]
    plot.bar(variable,height=value,width=0.8)



def varying_parameters(parms_min,parms_median,parms_max,equation,abs_req):
    """repeat loop_diff for a range of parameters and keep output"""
    rank_list=[]
    parms_temp=parms_median
    for k in equation.free_symbols:
       # if k!=cr:
        # if k!=q:
        param_value_range=np.linspace(parms_min[str(k)],parms_max[str(k)],no_range)
        for l in param_value_range:
            parms_temp[str(k)]=l
            rank_list.append([k,l,loop_diff(parms_temp,equation,abs_req)])
    return(rank_list)

def varying_parameters_ranks(parms_min,parms_median,parms_max,equation):
    """repeat loop_diff for a range of parameters and keep output- ranks """
    rank_list=[]
    parms_temp=parms_median
    for k in equation.free_symbols:
     #   if k!=cr:
       #     if k!=q:
        param_value_range=np.linspace(parms_min[str(k)],parms_max[str(k)],no_range)
        for l in param_value_range:
            parms_temp[str(k)]=l
            argh=loop_diff(parms_temp,equation,abs_req=True)
            argh={key: rank for rank, key in enumerate(sorted(argh, key=argh.get, reverse=False), 1)}
            rank_list.append([k,l,argh])
    return(rank_list)


def plot_value_ranges(a):
    """plots the mean and standard error of dR0/dt based on the results of the sensitivity analysis"""
    value_dics=[i[2] for i in a] #just values 
    parameters=[i for i in value_dics[0]] # just parameters
    parameters.sort()

    just_values=[value_dics[l][s] for s in parameters for l in range(len(parameters*no_range))] #all possible values
    parameters=np.array(parameters)
    no_parameters=len(parameters)
    no_each_parameter=len(just_values)/no_parameters
    no_each_parameter=int(no_each_parameter)
    array_values=np.array(just_values)    
    to_fill = np.reshape(just_values, (no_parameters,no_each_parameter))
    m=np.matrix(to_fill,dtype=np.float64)
    means=np.array(m.mean(1))
    errors=np.array(m.std(1))
    array_means=np.repeat(means,no_each_parameter)
    array_parameters=np.repeat(parameters, no_each_parameter)
    array_errors=np.repeat(a=errors,repeats=no_each_parameter)
    plot.bar(array_parameters,array_means,yerr=array_errors)
    plot.xlabel("Parameter")
    

sigma, mu, cr, epsilon, d, h, gamma, p = var('sigma, mu, cr, epsilon, d, h, gamma, p ',real=True, positive = true)
q = var('q',real=True, positive = False)
R0 = (sigma / (sigma + mu)) * cr * epsilon *  (( (-1 / q) * (1 - exp( q * d)))/((- 1 / q) * (1 - exp( q * d)) + h)) / ((mu + gamma)*( 1 / ( 1 - p )))


parms_median ={"sigma": 0.68 , "mu": 2.06e-5, "cr": 18.54 , "epsilon": 0.05, "d": 4/24, "h": 0.25/24, "gamma": 0.25, "p": 0.005, "q":-39.87 }
parms_max ={"sigma": parms_median["sigma"]*1.5, "mu": parms_median["mu"]*1.5, "cr": parms_median["cr"]*1.5 , "epsilon": parms_median["epsilon"]*1.5, "d": parms_median["d"]*1.5, "h": parms_median["h"]*1.5, "gamma": parms_median["gamma"]*1.5, "p": parms_median["p"]*1.5, "q":parms_median["q"]*1.5}
parms_min ={"sigma": parms_median["sigma"]*0.5, "mu": parms_median["mu"]*0.5, "cr": parms_median["cr"]*0.5 , "epsilon": parms_median["epsilon"]*0.5, "d": parms_median["d"]*0.5, "h": parms_median["h"]*0.5, "gamma": parms_median["gamma"]*0.5, "p": parms_median["p"]*0.5, "q":parms_median["q"]*0.5}

a1=varying_parameters(parms_min,parms_median,parms_max,R0,abs_req=False)

plotting=plot.figure()
plot_value_ranges(a1)
plot.ylabel("psi")
plotting.savefig("../../Writeup/sensitivity.pdf")
#d(R0)/d(Parameter) /R0
