import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time


parser = argparse.ArgumentParser(description='Takes in a distribution of permuted statistic values, a nominal statistic value, plots them, and gets a p-value.') 


parser.add_argument('datafile', help='Name of permutation values data file')
parser.add_argument('ksvalue', help='Nominal statistic value from data')
parser.add_argument('testname', help='Name of test statistic')
parser.add_argument('--show_plots', dest='show_plots', default='False')

args = parser.parse_args()

datafile = args.datafile
ksvalue = float(args.ksvalue)
show_plots = args.show_plots
testname = args.testname


with open(datafile) as datafile: 
    lines=datafile.readlines()[1:]
    ksdist = [float(line) for line in lines]

print('Finished reading file.')
print('File: ' + str(datafile))
print('This is the ' + testname + ' test statistic.')
print('Nominal: ' + str(ksvalue) + '\n')


j = 0
for i in ksdist: 
    if i >= ksvalue: 
        j = j+1


pvalue = j/len(ksdist)

print('>>> The p value is ' + str(pvalue) + ' <<< \n')

if show_plots in ['True', 'y', '1']: 

    plt.figure()

    text = 'Nominal ' +testname+' statistic: ' + str(ksvalue) +'\n'+ 'p-value: ' + str(pvalue) 

    plt.hist(ksdist, bins=120, log=True) #int(len(ksdist)/15), log=False)
    plt.title(testname + ' distribution')
    plt.axvline(x=ksvalue, color='r', linewidth=1, ymax=0.75)
    plt.text(0.5,0.9,text,transform=plt.gca().transAxes)

    plt.show()






