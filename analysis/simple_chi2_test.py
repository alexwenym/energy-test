import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats

#file1 = str(sys.argv[1])
#file2 = str(sys.argv[2])
#xbins = int(sys.argv[3])
#ybins = int(sys.argv[4])

matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000



def chi2_test(file1, file2, xbins, ybins):

    print('Reading file ' + file1 + ' as file1 and file ' + file2 + ' as file2.')
    print('Using ' + str(xbins) + ' xbins and ' + str(ybins) + ' ybins for a total of ' + str(xbins*ybins) + ' bins.')



    # designate sample1 as X and sample2 as Xbar
    # therefore alpha is N(s1)/N(s2)

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

    print(s1)
    print(s1[1][0])

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

    print('the Chi2 value is ' + str(chi2)) 
    
    p = scipy.stats.chi2.sf(chi2, xbins*ybins)

    print('there are ' + str(xbins*ybins) + ' dof')

    print('the p value is ' + str(p))
    
    return chi2, p 




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



#data=[float(i) for i in data]
#print(type(data[1][1]))
#plot = plt.scatter(s1_m12, s1_m13, cmap="Purples")
#plt.pcolormesh(data)
#plt.show()
