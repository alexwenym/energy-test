import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random 

file1 = str(sys.argv[1])
file2 = str(sys.argv[2])
size = int(sys.argv[3])

print('Using a subsampled size of ' + str(size) + ' from the distance distribution')

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


def get_distance_sample(sample, number): 
    '''
    Gets a specified number sample of distances. 
    input: sample is an array, from which to sample from (in this case s1/s2). 
    number is the number of distances to output. For a given sample of size S,
    the number of unique distances within the sample is S*(S-1)/2, or S Choose 2.
    '''
    distances = []

    for i in range(number):

        indexes = np.random.choice(len(sample), 2, replace=False)
    
        v1 = sample[indexes[0]]
        v2 = sample[indexes[1]]

        distance = calculate_distance(v1, v2)

        distances.append(distance)

    return distances



s1, s2 = open_files(file1, file2)

s1, s2 = get_permuted_samples(s1, s2)

s1distances = get_distance_sample(s1, size)
s2distances = get_distance_sample(s2, size)

ks, p = scipy.stats.ks_2samp(s1distances, s2distances)

print('the p value from the ks 2 sample test is ' + str(p))
print('the ks test statistic is ' + str(ks))




plt.figure()

plt.subplot(1,2,1)
plt.hist(s1distances, bins=200)
plt.title('sample1 distance distribution')

plt.subplot(1,2,2)
plt.hist(s2distances, bins=200)
plt.title('sample2 distance distribution')

plt.show()



'''
p_array = []
chi2_array = []
r=20

for i in range(r): 

    a, b = chi2_test('sample1.txt', 'sample2.txt', i, i)

    print('reached ' + str(i))

    p_array.append(b)
    chi2_array.append(a)


plt.figure()

plt.subplot(1,2,1)
plt.plot(list(range(r)), p_array) 
plt.title("p value vs bins on each axis")

plt.subplot(1,2,2)
plt.plot(list(range(r)), chi2_array)
plt.title("chi2 statistic vs. bins on each axis")

plt.show()
'''


#data=[float(i) for i in data]
#print(type(data[1][1]))
#plot = plt.scatter(s1_m12, s1_m13, cmap="Purples")
#plt.pcolormesh(data)
#plt.show()
