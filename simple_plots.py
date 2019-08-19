import numpy as np
from numpy import asarray
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time
from mpl_toolkits.mplot3d import Axes3D
 
parser = argparse.ArgumentParser(description='Compares distance distributions used for the Energy Test by using the Chi2 test.')

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#size = int(sys.argv[3])

parser.add_argument('file1', help='Name of first sample file')
parser.add_argument('file2', help='Name of second sample file')
#parser.add_argument('distsize', help='Number of distances to subsample for each test. Since the total number of distances is huge, this is a faster way to it without using all distances (order of 10e10 distances for 10e5 events.)')


args = parser.parse_args()
file1 = args.file1
file2 = args.file2


#print('Using a subsampled size of ' + str(distsize) + ' from the distance distribution')

matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000


def open_files(file1, file2): 
    '''
    Opens the two files and outputs them as two 
    N by 3 arrays with each row containing
    m12, m13, and m23. 
    '''
    
    #print('Reading ' + file1 + ' as file1 and ' + file2 + ' as file2.')

    with open(file1) as file1:
        lines = file1.readlines()
        s1m12 = [line.split()[0] for line in lines]
        s1m13 = [line.split()[1] for line in lines]
        s1m23 = [line.split()[2] for line in lines]

    with open(file2) as file2:
        lines = file2.readlines()
        s2m12 = [line.split()[0] for line in lines]
        s2m13 = [line.split()[1] for line in lines]
        s2m23 = [line.split()[2] for line in lines]


    s1m12 = [float(i) for i in s1m12]
    s1m13 = [float(i) for i in s1m13]
    s1m23 = [float(i) for i in s1m23]
    s1 = np.column_stack((s1m12, s1m13, s1m23))

    s2m12 = [float(i) for i in s2m12]
    s2m13 = [float(i) for i in s2m13]
    s2m23 = [float(i) for i in s2m23]
    s2 = np.column_stack((s2m12, s2m13, s2m23))

    #print('Finished reading files.')

    return s1, s2


def get_permuted_samples(s1, s2):
    '''
    Given two input samples, returns two output ones of the same shape
    that have shuffled contents; ie the two output ones should have the same
    distributions, even if the two input ones were distinct.
    '''

    s12 = np.concatenate((s1, s2), axis=0)

    np.random.shuffle(s12)

    s1out = s12[0:len(s1)]
    s2out = s12[len(s1):len(s12)]

    return s1out, s2out





# Function definitions end here
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
# Function calling starts here

s1, s2 = open_files(file1, file2)

fig1 = plt.figure()
ax1 = plt.axes()
ax1.scatter(s1[:,0], s1[:,1], s=0.1)

fig2 = plt.figure()
ax2 = plt.axes()
ax2.scatter(s2[:,0], s2[:,1], s=0.1)

plt.show()


#chi2values = []

#chi2valuenominal = chi2test(s1, s2, nbins)

#s1distances, s2distances = distances.get_distances(s1, s2)

#ksvaluenominal=5
#print('\n >>> the nominal chi2 value is ' + str(chi2valuenominal) + ' <<< \n')

#print('Running the Chi2 test ' + str(npermutations) + ' time(s)') 

#for i in range(npermutations): 

#    if i % 50 == 0: 
#        print('Reached permutation #' + str(i))

#    s1, s2 = get_permuted_samples(s1, s2)
    
#    chi2value = chi2test(s1, s2, nbins)

    #print('the ks value is ' + str(ksvalue))
    
#    chi2values.append(chi2value)

#print('The ks value is ' + str(sum(ksvalues)/len(ksvalues)))

#plt.figure()
#plt.hist(ksvalues, bins=25)
#plt.show()

#chi2values = np.array(chi2values)

#j=0
#for i in chi2values: 
#    if i > chi2valuenominal: 
#        j=j+1

#if len(chi2values) != 0:
#    pvalue = j/len(chi2values)
#else: 
#    pvalue = 'NaN'

#np.savetxt(outfile, chi2values, delimiter=' ', header=str(chi2valuenominal)+' p='+str(pvalue) )
print('Script Ended.')

#plt.figure()

#plt.subplot(1,2,1)
#plt.hist(s1distances, bins=nbins, log=False, alpha=0.7)
#plt.title('sample1 distance distribution')

#plt.subplot(1,2,2)
#plt.hist(s2distances, bins=nbins, log=False, alpha=0.7)
#plt.title('sample2 distance distribution')

#plt.show()
