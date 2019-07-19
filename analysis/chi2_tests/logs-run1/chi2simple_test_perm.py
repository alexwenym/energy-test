import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random 

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#xbins = int(sys.argv[3])
#ybins = int(sys.argv[4])

matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000



def get_chi2(s1, s2, xbins, ybins):
    '''
    Input: s1 and s2 are L by 2 arrays with each row containing a pair of invariant 
    masses for a decay into a three body system. 
    xbins and ybins refer to the number of bins along each axis. 
    Output: the chi-square test statistic for this particular sample. 
    '''

    # designate sample1 as X and sample2 as Xbar
    # therefore alpha is N(s1)/N(s2)

    alpha = len(s1)/len(s2)

    Sarray = np.array([]).reshape(0,1)

    for i in range(xbins): 

        s1xcut = [entry for entry in s1 if i*(1/xbins) <= entry[0] < (i+1)*(1/xbins)]
        s2xcut = [entry for entry in s2 if i*(1/xbins) <= entry[0] < (i+1)*(1/xbins)]
    
        for j in range(ybins): 
    
            s1ycut = [entry for entry in s1xcut if j*(1/ybins) <= entry[1] < (j+1)*(1/ybins)]
            s2ycut = [entry for entry in s2xcut if j*(1/ybins) <= entry[1] < (j+1)*(1/ybins)]

            s1ycut = np.array(s1ycut)
            s2ycut = np.array(s2ycut)
               
            p = (len(s1ycut)-(alpha*len(s2ycut)))
            q = np.sqrt(len(s1ycut)+ (len(s2ycut)*alpha**2))
        
            if q != 0: 

                #plt.scatter(s1ycut[:,0], s1ycut[:,1])
                #plt.show()

                S = p/q
        
                Sarray = np.append(Sarray, np.array(S).reshape(1,1), axis=0)
    
    S2 = np.square(Sarray)

    chi2 = np.sum(S2)

    #print('the Chi2 value is ' + str(chi2)) 
    
    #p = scipy.stats.chi2.sf(chi2, xbins*ybins)

    #print('there are ' + str(xbins*ybins) + ' dof')

    #print('the p value (from scipy) is ' + str(p))
    
    return chi2



def permutation(file1, file2, xbins, ybins, iterations):
    '''
    Gets a distribution of chi-square values to mimic the null hypothesis. 
    Input: two input files containing two different sets of samples, the number of
    bins along each axis, and the number of iterations - the number of data
    points you need in the final distribution. Could be costly to run. Need to optimize

    Output: list of chi-square values that mimic the null hypothesis, as well as 
    the chi-square value that comes from the differentiated dataset. 
    Consequently there is an estimate of the p-value.
    '''

    print('Doing permutations. Will do ' + str(iterations))
    print('Reading file ' + file1 + ' as file1 and file ' + file2 + ' as file2.')

    with open(file1) as file1:
        lines = file1.readlines()
        s1m12 = [line.split()[0] for line in lines]
        s1m13 = [line.split()[1] for line in lines]

    with open(file2) as file2:
        lines = file2.readlines()
        s2m12 = [line.split()[0] for line in lines]
        s2m13 = [line.split()[1] for line in lines]


    s1m12 = [float(i) for i in s1m12]
    s1m13 = [float(i) for i in s1m13]
    s1 = np.column_stack((s1m12, s1m13))

    s2m12 = [float(i) for i in s2m12]
    s2m13 = [float(i) for i in s2m13]
    s2 = np.column_stack((s2m12, s2m13))

    true_chi2 = get_chi2(s1,s2,xbins,ybins)
    print('The real chi2 value is ' + str(true_chi2))
    scipy_p = scipy.stats.chi2.sf(true_chi2, xbins*ybins)
    print('the p value as estimated by scipy.stats is ' + str(scipy_p))


    s12 = np.concatenate((s1, s2), axis=0)
    chi2_dist = []

    print('reached permutation loop')

    for i in range(iterations):

        if i % 20 == 0: 
            print('reached ' + str(i) + ' events')

        np.random.shuffle(s12)

        s1 = s12[0:len(s1)]
        s2 = s12[len(s1):len(s12)]

        chi2 = get_chi2(s1, s2, xbins, ybins)

        chi2_dist.append(chi2)

    
    return chi2_dist, true_chi2


dist = permutation('sample1.txt', 'sample2.txt', 3, 3, 10000)

counter = 0 
for i in dist[0]: 
    if i > dist[1]: 
        counter = counter+1
 
p_value = counter/10000

print('the p value from the permutation method is ' + str(p_value))


plt.figure()
plt.hist(dist[0], bins=50)
plt.title("chi2 distribution")
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
