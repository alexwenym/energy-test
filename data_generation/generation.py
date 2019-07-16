import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy.stats
import numpy.random
import argparse
import time
from mpl_toolkits.mplot3d import Axes3D


parser = argparse.ArgumentParser(description='Generates a simple custom data set, one uniform and the other with some degree of ununiformity. Meant as a simple set of toy data for analysis by experimental statistical tests.')

parser.add_argument('nevents', help='Number of events to generate in each sample. Same for both samples.')
parser.add_argument('outfile1', help='Name of file to write first output sample to.')
parser.add_argument('outfile2', help='Name of file to write second output sample to.') 
parser.add_argument('eccratio', help='Eccentricity ratio - ratio of sample 2 events inside a smaller area to the number of total events. Crude way of controlling the similarity between the two samples.')

args = parser.parse_args()
outfile1 = args.outfile1
outfile2 = args.outfile2
nevents = int(args.nevents)
eccratio = float(args.eccratio)

print('\nWarning: NOT optimized. I\'m not intending to generate more than a few thousand events per sample. If heavy lifting is required, I suggest optimizing the generation to be faster. \n')

print('Generating ' + str(nevents) + ' events.')

def generate_uniform(nevents): 

    l= np.random.uniform(0,1,size=(nevents,3))

    z = np.zeros((nevents, 2))

    l = np.concatenate((l,z), axis=1)

    return l

def generate_eccentric(nevents, eccratio): 
    
    normalnum = int(nevents*(1-eccratio)) 
    eccnum = int(nevents*eccratio)
    if normalnum+eccnum > nevents: 
        eccnum = eccnum-1
    if normalnum+eccnum < nevents: 
        eccnum = eccnum+1
    assert normalnum+eccnum == nevents

    l = np.random.uniform(0,1,size=(normalnum, 3))
    s = np.random.uniform(0,0.2,size=(eccnum, 3))

    ls = np.concatenate((l,s), axis=0)

    z = np.zeros((nevents, 2))

    ls = np.concatenate((ls, z), axis=1)

    np.random.shuffle(ls)
    #print(len(ls))

    return ls

    
s1 = generate_uniform(nevents)

s2 = generate_eccentric(nevents, eccratio)

np.savetxt(outfile1, s1, delimiter=' ')
np.savetxt(outfile2, s2, delimiter=' ')

#fig = plt.figure()

#ax = plt.axes(projection='3d')

#ax.scatter(s2[:,0], s2[:,1], s2[:,2])

#plt.show()



