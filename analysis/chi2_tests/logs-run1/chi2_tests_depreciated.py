import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time

 
parser = argparse.ArgumentParser(description='Obtains the p-value for a pair oof distance distribution samples using the chi2 test statistic.')

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#size = int(sys.argv[3])

parser.add_argument('file1', help='Name of first sample file')
parser.add_argument('file2', help='Name of second sample file')
parser.add_argument('--nbins', dest='nbins', default=10, help='Number of bins to use for the chi2 test')
parser.add_argument('--npermutations', dest='npermutations', default=1, help='Number of permutations. Keep in mind this is very time-intensive.')
parser.add_argument('--show_plots', dest='show_plots', default='False')
parser.add_argument('--run_nominal', dest='run_nominal', default='False')

args = parser.parse_args()
file1 = args.file1
file2 = args.file2
npermutations = int(args.npermutations)
show_plots = args.show_plots
nbins = int(args.nbins) 
run_nominal = args.run_nominal


matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000

print('WARNING: this assumes that your sample sizes are equal. If they are not, manually introduce an alpha scaling factor in the chi2test() function.')


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
    Gets array of distances for the given sample. 
    input: sample is an array, from which to get distances from (in this case s1/s2) 
    the number of unique distances within the sample is S*(S-1)/2, or S Choose 2. For 50k event samples, this is on the order of 10e9 (will probably crash)
    '''
    distances = []
    
    start=time.time()

        #indexes = np.random.choice(len(sample), 2, replace=False) # <-- SO SLOW
        #indexes = [i, i + 1]
        #L = len(sample)
        #indexes = np.random.randint(0, L, size=2) # <-- MUCH FASTER (>10x)
        #while indexes[0] == indexes[1]: 
        #    indexes = np.random.randint(0, L, size=2)

    for j in range(len(sample)): 

        #if j % 100 == 0: 
            #print('outer loop reached event' + str(j))
            
        for k in range(j+1,len(sample)): 
    
            v1 = sample[j]
            v2 = sample[k]

            distance = calculate_distance(v1, v2)

            distances.append(distance)

    end=time.time()
    #print('Finished getting distances in one sample. That took ' + str(end-start) + ' seconds.')
    #print('the number of distances is ' + str(len(distances)))

    return distances


def chi2test(s1, s2, nbins, nominal):    
    
    if nominal == True: 
        pass
    else:
        s1, s2 = get_permuted_samples(s1, s2)

    s1distances = get_distances(s1)
    s2distances = get_distances(s2)
    

    Sarray = []

    hS1 = np.histogram(s1distances, bins=nbins, density=False, range=(0,np.sqrt(3)))
    hS2 = np.histogram(s2distances, bins=nbins, density=False, range=(0,np.sqrt(3)))

    for i in range(nbins): 
        S1 = hS1[0][i]
        S2 = hS2[0][i]

        #print('reached ' + str(i) + ' with value ' + str(S1) + ' and ' + str(S2))

        if S1 == 0 and S2 == 0:
            pass
        else:
            S = ((S1-S2)**2)/(S1+S2)
            Sarray.append(S)

    statistic = sum(Sarray)

    #print('the p value from the K-S 2 sample test is ' + str(p))
    #print('the K-S test statistic is ' + str(ks))

    return statistic


# Function definitions end here
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
# Function calling starts here

s1, s2 = open_files(file1, file2)

chi2list = []

print('Running the chi2 test with ' + str(npermutations) + ' permutations.')

if run_nominal == 'True': 
    nominalchi2 = chi2test(s1, s2, nbins, nominal=True)
    print('The nominal chi2 value is ' + str(nominalchi2))


for i in range(npermutations): 

    if i % 20 == 0: 
        print('Reached test #' + str(i))

    chi2value = chi2test(s1, s2, nbins, nominal=False)

    #print('the chi2 value is ' + str(ksvalue))
    
    chi2list.append(chi2value)

#print('The ks value is ' + str(sum(ksvalues)/len(ksvalues)))

#plt.figure()
#plt.hist(ksvalues, bins=25)
#plt.show()

chi2list = np.array(chi2list)

np.savetxt('chi2permuted_'+str(nbins)+'bins.txt', chi2list, delimiter=' ')
print('Files successfully saved.')

if show_plots in ['True', 'y', '1']: 
    plt.figure()

    plt.subplot(1,2,1)
    plt.hist(s1distances, bins=numbins, log=False)
    plt.title('sample1 distance distribution')

    plt.subplot(1,2,2)
    plt.hist(s2distances, bins=numbins, log=False)
    plt.title('sample2 distance distribution')

    plt.show()
