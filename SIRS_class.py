import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.pylab as plt
import random
import sys
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation


#class to simulate and perform measurements on the SIRS model

class SIRS(object):

#Initialises a square array of user defined size
    def __init__(self,dimension,p1,p2,p3):
        self.array = np.random.choice([-1,1,0],size=(dimension,dimension))
        self.dimension = dimension
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.N = int(dimension**2)

#Instance method which imposes periodic boundary conditions on a specific index
    def PBC(self,index):
        max = self.dimension - 1
        if index > max:
            index = 0
        elif index < 0:
            index = max
        return index

#Instance method which returns the 4 nearest neighbours of a given index in the SIRS array employing periodic boundary conditions.
    def NNs(self,a,b):
        nn1 = self.array[self.PBC(a-1),self.PBC(b)]
        nn2 = self.array[self.PBC(a+1),self.PBC(b)]
        nn3 = self.array[self.PBC(a),self.PBC(b-1)]
        nn4 = self.array[self.PBC(a),self.PBC(b+1)]

        return [nn1,nn2,nn3,nn4]

#Static method which trials a given probability
    @staticmethod
    def Probability(p):
        r = np.random.uniform(0,1)
        if p > r: return 1
        else: return 0
    
#Instance method which imposes the rules of the SIRS model on a given index in the array.
    def Rules(self,a,b):
        nn_list = self.NNs(a,b)
        #Dict item which contains the number of each type of nns of an index
        c = Counter(nn_list)

        if self.array[a,b] == 1 and c[0]>=1 and SIRS.Probability(self.p1)==True:
            self.array[a,b] = 0

        elif self.array[a,b] == 0 and SIRS.Probability(self.p2) == True:
            self.array[a,b] = -1

        elif self.array[a,b] == -1 and SIRS.Probability(self.p3)==True:
            self.array[a,b] = 1

#Instance method which randomly samples the array N times and simulates the SIRS model
    def Sweep(self):
        for i in range(self.N):
            a = int(np.random.uniform(0,self.dimension))
            b = int(np.random.uniform(0,self.dimension))
            self.Rules(a,b)

#Instance methods which returns the number of infected cells in the array
    def InFraction(self):
        inf = np.count_nonzero(self.array == 0)
        #print("Number of infected: "+str(inf))
        return inf

#Instance method which takes a list of samples of infected cells and the sqaured of the list. Returns the variance/N and avg number of infected cells/N
    def Variance(self,inf,infsq):
        avg_inf = np.average(inf)
        avg_infsq = np.average(infsq)
        var = (avg_infsq - avg_inf**2)/(self.N)
        return var, avg_inf/self.N

#Static method which generates a contour plot for lists x and y and an array z. Dimensions must agree
    def Contour(x,y,z,title):
        X,Y = np.meshgrid(x,y)
        plt.contourf(X,Y,z,30)
        plt.xlabel("P1 (S->I)")
        plt.ylabel("P3 (R->S)")
        plt.title("Phase diagram of " +str(title))
        plt.colorbar();
        plt.show()

    def Set_Fraction(self,fraction):
        indices = []
        L = list(range(self.dimension))
        while len(indices)<int(self.N*fraction):
            x = random.choice(L)
            y = random.choice(L)
            item = [x,y]
            if item not in indices:
                indices.append(item)
        for i in indices:
            self.array[i[0],i[1]] = 2
