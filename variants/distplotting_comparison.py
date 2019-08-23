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
parser.add_argument('datafile2', help='Name of p values data file 2')
parser.add_argument('testname', help='Name of test statistic')

args = parser.parse_args()

datafile1 = args.datafile1
datafile2 = args.datafile2
testname = args.testname


with open(datafile1) as datafile1: 
    lines=datafile1.readlines()[1:]
    ksdist1 = [float(line) for line in lines]

with open(datafile2) as datafile2: 
    lines=datafile2.readlines()[1:]
    ksdist2 = [float(line) for line in lines]

print('Finished reading file.')
print('File: ' + str(datafile1) + ' and ' + str(datafile2))
print('This is the ' + testname + ' test statistic.')

sigma1 = 0
sigma2 = 0
sigma3 = 0

for i in ksdist1: 
    if i < 0.15865525: 
        sigma1 = sigma1+1
    if i < 0.02275013: 
        sigma2 = sigma2+1 
    if i < 0.00134989: 
        sigma3 = sigma3+1 

print('datafile1 1 sigma: ' + str(sigma1))
print('datafile1 2 sigma: ' + str(sigma2))
print('datafile1 3 sigma: ' + str(sigma3))

sigma1 = 0
sigma2 = 0
sigma3 = 0

for i in ksdist2: 
    if i < 0.15865525: 
        sigma1 = sigma1+1
    if i < 0.02275013: 
        sigma2 = sigma2+1 
    if i < 0.00134989: 
        sigma3 = sigma3+1 

print('datafile2 1 sigma: ' + str(sigma1))
print('datafile2 2 sigma: ' + str(sigma2))
print('datafile2 3 sigma: ' + str(sigma3))

diffcounter = 0
for i in range(len(ksdist1)): 
    if ksdist1[i] < ksdist2[i]:
        diffcounter = diffcounter+1 

print('There are ' + str(diffcounter) + ' events in sample2 which are greater than their counterpart in sample1')

diff = np.subtract(ksdist1, ksdist2)

print('datafile1 minus datafile2')
print(diff)

plt.figure()

plt.hist(diff, bins=1000)

plt.show()


#plt.figure() 

#logbins = np.logspace(np.log10(0.00001),np.log10(0.5),50)

#plt.hist(ksdist, bins=logbins, range=[0,0.1], log=True) #int(len(ksdist)/15), log=False)
#plt.title(testname + ' distribution')

#plt.show()






