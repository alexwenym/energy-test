import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time

datafile1 = '500eventdata.txt'
datafile2 = '1000eventdata.txt'
datafile3 = '2000eventdata.txt'


#500 event samples
with open(datafile1) as datafile:

    lines=datafile.readlines()
    
    
    nominal1 = float(lines[0][2:])
    ksdist1 = [float(line) for line in lines[1:]]
    
    nominal1=nominal1*1
    ksdist1=ksdist1*1

#1000 event samples
with open(datafile2) as datafile: 

    lines=datafile.readlines()

    nominal2 = 1*float(lines[0][2:])
    ksdist2 = [1*float(line) for line in lines[1:]]


#2000 event samples
with open(datafile3) as datafile: 

    lines=datafile.readlines()

    nominal3 = 1*float(lines[0][2:])
    ksdist3 = [1*float(line) for line in lines[1:]]


print('Finished reading files.')
#print('File: ' + str(datafile))
#print('This is the ' + testname + ' test statistic.')
#print('Nominal: ' + str(ksvalue) + '\n')


j = 0
for i in ksdist1: 
    if i >= nominal1: 
        j = j+1
k = 0 
for i in ksdist2: 
    if i >= nominal2: 
        k = k+1

l = 0 
for i in ksdist3: 
    if i >= nominal3:
        l = l+1


pvalue1 = j/len(ksdist1)
pvalue2 = k/len(ksdist2)
pvalue3 = l/len(ksdist3)

print('>>> The p value for the 500 event run is ' + str(pvalue1) + ' <<< \n')
print('>>> The p value for the 1000 event run is ' + str(pvalue2) + ' <<< \n')
print('>>> The p value for the 2000 event run is ' + str(pvalue3) + ' <<< \n')


plt.figure()

#text = 'Nominal ' +testname+' statistic: ' + str(ksvalue) +'\n'+ 'p-value: ' + str(pvalue) 

plt.hist(ksdist1, bins=int(len(ksdist1)/100), log=False, label='500 events, npermuations=10k', alpha=0.5, density=True)
plt.hist(ksdist2, bins=int(len(ksdist2)/100), log=False, label='1000 events, npermutations=10k', alpha=0.4, density=True)
plt.hist(ksdist3, bins=int(len(ksdist3)/16), log=False, label='2000 events, npermutations=1k', alpha=0.3, density=True)
plt.title('ks distributions, NORMALIZED')
plt.legend(loc='upper right')
plt.axvline(x=nominal1, color='blue', linewidth=1, ymax=0.75)
plt.axvline(x=nominal2, color='orange', linewidth=1, ymax=0.75)
plt.axvline(x=nominal3, color='green', linewidth=1, ymax=0.75)
#plt.text(0.5,0.9,text,transform=plt.gca().transAxes)

plt.show()






