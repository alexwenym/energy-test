import numpy as np
from numpy import asarray
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time
import distances

 
parser = argparse.ArgumentParser(description='Plots and compares distance distributions used for the Energy Test by using the Kolmogorov-Smirnov test.')

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#size = int(sys.argv[3])

parser.add_argument('file1', help='Name of first sample file')
parser.add_argument('file2', help='Name of second sample file')
#parser.add_argument('distsize', help='Number of distances to subsample for each test. Since the total number of distances is huge, this is a faster way to it without using all distances (order of 10e10 distances for 10e5 events.)')
parser.add_argument('npermutations', help='Number of permutations.`')
parser.add_argument('outfile', help='name of file to write output KS permutation distribution to.')
parser.add_argument('--show_plots', dest='featuretrue', action='store_true')
parser.add_argument('--numbins', dest='numbins', default=50)
parser.set_defaults(featuretrue=False)


args = parser.parse_args()
file1 = args.file1
file2 = args.file2
#distsize = int(args.distsize)
npermutations = int(args.npermutations)
#show_plots = args.show_plots
featuretrue = args.featuretrue
numbins = int(args.numbins) 
outfile = args.outfile


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



def ks_2samp(data1, data2):
    
    data1, data2 = map(asarray, (data1, data2))
    n1 = data1.shape[0]
    n2 = data2.shape[0]

    n1 = len(data1)
    n2 = len(data2)

    data1 = np.sort(data1)
    data2 = np.sort(data2)

    data_all = np.concatenate([data1,data2])

    cdf1 = np.searchsorted(data1,data_all,side='right')/(1.0*n1)
    cdf2 = (np.searchsorted(data2,data_all,side='right'))/(1.0*n2)

    d = np.max(np.absolute(cdf1-cdf2)) 

    return d




'''
def calculate_distance(v1, v2): 

    dij = v1-v2

    dij = np.sqrt(dij.dot(dij))

    return dij
'''

'''
def get_distances(s1, s2): 
    
    #Gets a specified number sample of distances. 
    #input: sample is an array, from which to sample from (in this case s1/s2) 
    #number is the number of distances to output. For a given sample of size S,
    #the number of unique distances within the sample is S*(S-1)/2, or S Choose 2. For 50k events, this is on the order of 10e9
    
    da1 = []
    da2 = []

    assert len(s1) == len(s2)
    
    start=time.time()
    #print('Starting loop to sample distances.')

        #indexes = np.random.choice(len(sample), 2, replace=False) # <-- SO SLOW
        #indexes = [i, i + 1]
        #L = len(sample)
        #indexes = np.random.randint(0, L, size=2) # <-- MUCH FASTER (>10x)
        #while indexes[0] == indexes[1]: 
        #    indexes = np.random.randint(0, L, size=2)

    for j in range(len(s1)): 

        #if j % 100 == 0: 
        #    print('outer loop reached event' + str(j))
            
        for k in range(j+1,len(s1)): 
    
            u1 = s1[j]
            u2 = s1[k]
            
            v1 = s2[j]
            v2 = s2[k]

            d1 = distances.calculate_distance(u1, u2)
            d2 = distances.calculate_distance(v1, v2)

            #distances1 = np.append(distances1, [distance1], axis=0)
            da1.append(d1)
            da2.append(d2)

    end=time.time()
    #print('Finished getting distances. That took ' + str(end-start) + ' seconds.')
    #print('the number of distances is ' + str(len(da1)))

    return da1, da2
'''


def kstest(s1, s2):

    s1distances, s2distances = distances.get_distances(s1, s2)

    ks = ks_2samp(s1distances, s2distances)

    #print('the p value from the K-S 2 sample test is ' + str(p))
    #print('the K-S test statistic is ' + str(ks))

    return ks


# Function definitions end here
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
# Function calling starts here

s1, s2 = open_files(file1, file2)

ksvalues = []

ksvaluenominal = kstest(s1, s2)

s1distances, s2distances = distances.get_distances(s1, s2)

#ksvaluenominal=5
print('\n >>> the nominal ks value is ' + str(ksvaluenominal) + ' <<< \n')

print('Running the KS test ' + str(npermutations) + ' time(s)') 

for i in range(npermutations): 

    if i % 50 == 0: 
        print('Reached permutation #' + str(i))

    s1, s2 = get_permuted_samples(s1, s2)
    
    ksvalue = kstest(s1, s2)

    #print('the ks value is ' + str(ksvalue))
    
    ksvalues.append(ksvalue)

#print('The ks value is ' + str(sum(ksvalues)/len(ksvalues)))

#plt.figure()
#plt.hist(ksvalues, bins=25)
#plt.show()

ksvalues = np.array(ksvalues)

np.savetxt(outfile, ksvalues, delimiter=' ', header=str(ksvaluenominal))
print('Script Ended. Files successfully saved.')

plt.figure()

plt.subplot(1,2,1)
plt.hist(s1distances, bins=numbins, log=False, alpha=0.7)
plt.title('sample1 distance distribution')

plt.subplot(1,2,2)
plt.hist(s2distances, bins=numbins, log=False, alpha=0.7)
plt.title('sample2 distance distribution')

plt.show()
