import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time

 
parser = argparse.ArgumentParser(description='Plots and compares distance distributions for the Energy Test')

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#size = int(sys.argv[3])

parser.add_argument('file1', help='Name of first sample file')
parser.add_argument('file2', help='Name of second sample file')
parser.add_argument('distsize', help='Number of distances to subsample for each test. Since the total number of distances is huge, this is a faster way to it without using all distances (order of 10e10 distances for 10e5 events.)')
parser.add_argument('--sampsize', dest='sampsize', default=1, help='Number of times to sample distances, ie. the number of data points you want, ie, the number of K-S test statistics you want.')
parser.add_argument('--permute', dest='permute', default='False')
parser.add_argument('--show_plots', dest='show_plots', default='False')
parser.add_argument('--numbins', dest='numbins', default=50)

args = parser.parse_args()
file1 = args.file1
file2 = args.file2
distsize = int(args.distsize)
sampsize = int(args.sampsize)
permute = args.permute
show_plots = args.show_plots
numbins = int(args.numbins) 


print('Using a subsampled size of ' + str(distsize) + ' from the distance distribution')



matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000


def open_files(file1, file2): 
    '''
    Opens the two files and outputs them as two 
    N by 3 arrays with each row containing
    m12, m13, and m23. 
    '''
    
    print('Reading ' + file1 + ' as file1 and ' + file2 + ' as file2.')

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

    print('Finished reading files.')

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
    


def calculate_distance(v1, v2): 

    a = v1[0] - v2[0]
    b = v1[1] - v2[1]
    c = v1[2] - v2[2]

    dij = np.sqrt(a**2 + b**2 + c**2)

    return dij


def get_distances(sample): 
    '''
    Gets a specified number sample of distances. 
    input: sample is an array, from which to sample from (in this case s1/s2) 
    number is the number of distances to output. For a given sample of size S,
    the number of unique distances within the sample is S*(S-1)/2, or S Choose 2. For 50k events, this is on the order of 10e9
    '''
    distances = []
    
    start=time.time()
    print('Starting loop to sample distances.')

        #indexes = np.random.choice(len(sample), 2, replace=False) # <-- SO SLOW
        #indexes = [i, i + 1]
        #L = len(sample)
        #indexes = np.random.randint(0, L, size=2) # <-- MUCH FASTER (>10x)
        #while indexes[0] == indexes[1]: 
        #    indexes = np.random.randint(0, L, size=2)

    for j in range(len(sample)): 

        #if j % 100 == 0: 
        #    print('outer loop reached event' + str(j))
            
        for k in range(j+1,len(sample)): 
    
            v1 = sample[j]
            v2 = sample[k]

            distance = calculate_distance(v1, v2)

            distances.append(distance)

    end=time.time()
    print('Finished sampling distances. That took ' + str(end-start) + ' seconds.')
    print('the number of distances is ' + str(len(distances)))

    return distances


def kstest(s1, s2, permute):

    #print('Running K-S test...')    

    if permute in ['True', 'y', '1']: 
        #print('WITH permutation')
        #print('Permuting samples...')
        s1, s2 = get_permuted_samples(s1, s2)
        #print('Samples are permuted. Expect similar distributions.')
    #else: 
        #print('WITHOUT permutation')


    s1distances = get_distances(s1)
    s2distances = get_distances(s2)

    ks, p = scipy.stats.ks_2samp(s1distances, s2distances)



    #print('the p value from the K-S 2 sample test is ' + str(p))
    #print('the K-S test statistic is ' + str(ks))

    return ks, p


# Function definitions end here
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
# Function calling starts here

s1, s2 = open_files(file1, file2)

ksvalues = []

print('Running the KS test ' + str(sampsize) + ' time(s), with ' + str(distsize) + ' distances each time.')

if permute in ['True', 'y', '1']:
    print('WITH permutation.')
else: 
    print('WITHOUT permutation')


for i in range(sampsize): 

    if i % 20 == 0: 
        print('Reached test #' + str(i))

    ksvalue, p = kstest(s1, s2, permute)

    print('the ks value is ' + str(ksvalue))
    
    ksvalues.append(ksvalue)

#print('The ks value is ' + str(sum(ksvalues)/len(ksvalues)))

#plt.figure()
#plt.hist(ksvalues, bins=25)
#plt.show()

ksvalues = np.array(ksvalues)

#np.savetxt('kspermuted.txt', ksvalues, delimiter=' ')
#print('Files successfully saved.')

if show_plots in ['True', 'y', '1']: 
    plt.figure()

    plt.subplot(1,2,1)
    plt.hist(s1distances, bins=numbins, log=False)
    plt.title('sample1 distance distribution')

    plt.subplot(1,2,2)
    plt.hist(s2distances, bins=numbins, log=False)
    plt.title('sample2 distance distribution')

    plt.show()
