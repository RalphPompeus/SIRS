import numpy as np
import sys
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation
from datetime import date
from SIRS_class import SIRS

######################
#p3 is on x-axis, p1 is on y-axis
#######################



#Defines global variables: number of iterations to collect data for and user defined dimensions and 3 probabilities
iterations = 10000
dimension = int(input('please enter a lattice dimension, N x N: '))
p1 = float(input('please enter a value for p1: '))
p2 = float(input('please enter a value for p2: '))
p3 = float(input('please enter a value for p3: '))
#p1,p2,p3 = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
task = input('please choose a task: viz, phase, immune: ')
task = task.lower()

# for waves of infection: p1 = 0.8, p2 = 0.095, p3 = 0.01
# for absorbing states: p1 = 0., p2 = 0.1, p3 = 1.
# for dynamical equilibrium: p1 = 0.5, p2 = 0.5, p3 = 0.5


#Animates the simulation of the SIRS model by calling on the SIRS class.
if task == 'viz':
    
    A = SIRS(dimension,p1,p2,p3)
    #A.Set_Fraction(0.3)
    
    def UpdatePlot(*args):
        image.set_array(A.array)
        A.Sweep()
        return image,

    SIRS_Image = plt.figure()
    image = plt.imshow(A.array,animated=True)
    model = FuncAnimation(SIRS_Image,UpdatePlot,interval=50,blit=True)
    plt.show()

#Simulates the SIRS model for a range of probabilities of S->I and I->R. Samples data and generates graphs of the average infected fraction/N and the variance in the infected fraction per N.
elif task == 'phase':
    
    #Creates arrays of probabilities to generate contour plots and creates empty square arrays for avg inf/N and avg var/N. Also lists for the same variables for writing to data file.
    
    steps = 2
    p1_list = np.linspace(p1,1,steps,endpoint=True)
    p3_list = np.linspace(p1,1,steps,endpoint=True)
    
    
    
    
    avg_inf_array = np.zeros(shape=(len(p1_list),len(p3_list)))
    avg_var_array = np.zeros(shape=(len(p1_list),len(p3_list)))
    avg_inf_list = []
    inf_var_list = []
    
    f = open('DataFile.txt','w')
    
    #For loops to iterate through probabilities
    for i in range(11):
        for j in range(11):
            
            #Reinitialises an instance of SIRS class to simulate with updated probabilities
            A = SIRS(dimension,i/10.,p2,j/10.)
            print (A.p1,A.p3)
            inf_list=[]
            infsq_list = []
            
            #Performs a predefined number of sweeps for each set of probs
            for t in range(iterations):
                A.Sweep()
                #If sufficient sweeps completed so equilibrium is reached, array is sampled for number of infected cells.
                if t>=100 and t%10 ==0:
                    inf = A.InFraction()
                    inf_list.append(inf)
                    infsq_list.append(inf**2.)
        
            #Calculates the variance and average infected fraction. Writes to data lists.
            var, inf_fraction = A.Variance(inf_list,infsq_list)
            avg_inf_list.append(inf_fraction)
            inf_var_list.append(var)
            
            #Assigns the calculated values to an index the respective arrays for contour plotting
            avg_inf_array[j,i] = inf_fraction
            avg_var_array[j,i] = var

            #Writes to text file
            f.write("{0:1.2}".format(A.p1)+" " +"{0:1.2}".format(A.p3)+" "+str(inf_fraction)+" "+str(var)+"\n")

    f.close()

    #Generates coloured contour plots
    SIRS.Contour(p1_list,p3_list,avg_inf_array,'Normalised number of infected cells as a function of p1 and p3')
    SIRS.Contour(p1_list,p3_list,avg_var_array,'Variance in normalised number of infected cells as a function of p1 and p3')

#Generates the graph for the average fraction of infected cells vs fraction of immune cells
elif task == 'immune':
    avg_inf_list = []
    immune_list = []
    error_list = []
    f = open('ImmunityFile.txt','w')

    #Iterates through lists of immune fractions
    for i in np.arange(0,1,0.01):
        print(i)
        inf_list = []
        #error_list = []
        A = SIRS(dimension,p1,p2,p3)
        A.Set_Fraction(i)
        for t in range(iterations):
            A.Sweep()
            if t>= 100 and t%10 ==0:
                inf = A.InFraction()
                inf_list.append(inf)
    
        if inf_list[-1]==0:
            avg_inf_list.append(0)
            
        else:
            avg_inf_list.append((np.average(inf_list))/A.N)
        immune_list.append(i)
        #print(immune_list)
        
        #error calculations
        #std_dev = (float((np.std(inf_list)/A.N))/float(np.sqrt(len(inf_list)))) 
        std_dev = (float((np.std(inf_list)/A.N)))    
        #std_error = np.std((len(inf_list)/dimension**2)/len(inf_list))**0.5
        print(std_dev)
        
        error_list.append(std_dev)
        
        
        f.write(str(i)+" "+str(np.average(inf_list))+"\n")
    plt.plot(immune_list, avg_inf_list)
    plt.errorbar(immune_list,avg_inf_list, yerr = error_list, fmt = 'none', capsize=4, capthick=1)
    plt.xlabel("Fraction of Permanently Immune Cells")
    plt.ylabel("Average Fraction of Infected Cells")
    plt.title("Plot of the Fraction of Infected cells vs Fraction of Immune cells")
    plt.savefig('immune_fraction_plot_50x50_'+str(date.today()),format='pdf')
    plt.show()
