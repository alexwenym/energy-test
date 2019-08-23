import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time


parser = argparse.ArgumentParser(description='') 


parser.add_argument('datafile1', help='Name of p values data file 1')
parser.add_argument('testname', help='Name of test statistic')

args = parser.parse_args()

datafile1 = args.datafile1
testname = args.testname


with open(datafile1) as datafile1: 
    lines=datafile1.readlines()[1:]
    ksdist = [float(line) for line in lines]

print('Finished reading file.')
print('File: ' + str(datafile1))
print('This is the ' + testname + ' test statistic.')

sigma1 = 0
sigma2 = 0
sigma3 = 0

significant = []; 
significant2 = [];

for i in ksdist: 
    if i < 0.15865525: 
        sigma1 = sigma1+1
    if i < 0.02275013: 
        significant.append(i)
        sigma2 = sigma2+1 
    if i < 0.00134989: 
        significant2.append(i)
        sigma3 = sigma3+1 

print('1 sigma: ' + str(sigma1/len(ksdist)))
print('2 sigma: ' + str(sigma2/len(ksdist)))
print('3 sigma: ' + str(sigma3/len(ksdist)))

average = sum(significant) / len(significant)
average2 = sum(significant2) / len(significant2)

print('there are ' + str(len(significant)) + ' events >2sigma')
print('average p value of >2sigma events is ' + str(average))

print('there are ' + str(len(significant2)) + ' events >3sigma')
print('average p value of >3sigma events is ' + str(average2))



plt.figure() 

logbins = np.logspace(np.log10(0.00001),np.log10(0.5),50)

plt.hist(ksdist, bins=logbins, range=[0,0.1], log=True) #int(len(ksdist)/15), log=False)
plt.title(testname + ' distribution')

plt.show()






