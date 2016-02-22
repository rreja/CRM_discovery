from optparse import OptionParser , IndentedHelpFormatter
import sys, os, re, operator, math
from random import shuffle
from collections import OrderedDict
from operator import add
import numpy as np
from pylab import *
import scipy.stats
from math import erf, sqrt


def  process_file_9(chromHmm_output):
    for direc in os.listdir(chromHmm_output):
        direc_path = os.path.join(chromHmm_output,direc)
        if not os.path.isdir(direc_path):
            continue
        
        for cdtDir in os.listdir(direc_path):
            if not cdtDir.startswith("_CDT"):
                continue
    
            xdata = []
            ydata = []
            significant = {}
            iterations = 1000
            observed = {}
            z = []
            data = {}
            data_to_shuffle = {}
            outfile = os.path.join(os.path.join(direc_path,cdtDir),"Means_vs_std_plot.svg")
            outfile2 = os.path.join(os.path.join(direc_path,cdtDir),"distribution_of_stds.txt")
            for fname in os.listdir(os.path.join(direc_path,cdtDir)):
                
                if not (fname.endswith('.cdt')):
                    continue
                
                in1 = open(os.path.join(os.path.join(direc_path,cdtDir),fname))
                label = fname.split("_")[0]
                print "Processing "+fname
                for line in in1:
                    if line.startswith("Uniqe") or line.startswith("EWEIGHT") or line.startswith("gene"):
                        tmp = line.rstrip().split("\t")[2:]
                        X = [int(x) for x in tmp]
                        xmin = min(X)
                        xmax = max(X)
                        Y = [0]*len(X)
                        continue
                    tmplist = line.rstrip().split("\t")[2:]
                    newList = [float(x) for x in tmplist]
                    Y = map(add,Y,newList)
                
                Y = movingaverage(Y,50)
                X = movingaverage(X,50)
                length = len(Y)
                
                X = [int(x) for x in X]
                xmin = min(X)
                xmax = max(X)
                #print xmin,xmax
                Y = [float(x)/max(Y) for x in Y]
                
                data_to_shuffle[label] = Y
                
                # Get the peak
                maximum = max(Y)
                # Find the corresponding index associated with the peak
                index_of_closest_value =  min(range(len(Y)), key=lambda i: abs(Y[i]-maximum))
                # Change the index to relative X-coordinate value.
                rel_value = index_of_closest_value + xmin
                std = calculate_std(rel_value,xmin,xmax,Y)
                #print std,label
                xdata.append(rel_value)
                ydata.append(std)
                z.append(label)
                observed[label] = std
                
            print ydata
            print z
            
            # Scatter plots
            colors = np.random.rand(50)
            fig, ax = plt.subplots()
            ax.scatter(xdata, ydata,c=colors,alpha=0.5)
            ax.set_xlim(-200,200)
            for i, txt in enumerate(z):
                ax.annotate(txt, (xdata[i],ydata[i]))
            savefig(outfile)
            
            
            # Shuffle original Y values
            for label,Y in data_to_shuffle.items():
                stds, count = create_shuffle_vectors(Y,xmin,xmax,iterations,ydata[z.index(label)],label)
                data[label] = stds
                #significant[label] = float(count)/iterations
                z_score = (observed[label] - np.mean(stds))/np.std(stds)
                #significant[label] = scipy.stats.norm.sf(abs(z_score))
                significant[label] = z2p(z_score)
            
            # create file for stds
            out2 = open(outfile2,"w")
            # writing header
            l = "fname"
            for k in z:
                l = l+"\t"+k
            out2.write(l+"\n")
            
            for j in range(0,iterations):
                line = "count"+str(j)
                for keys in z:
                    line = line+"\t"+str(data[keys][j])
                out2.write(line+"\n")
            
            print "The significant factors around this reference point are:"
            for k,v in significant.items():
                if v <= 0.05:
                    print k+"\t"+str(v)
    
            

def z2p(z):
    return 0.5 * (1 + erf(z / sqrt(2)))
   
            
def create_shuffle_vectors(copy_Y,xmin,xmax,iterations,obs,label):
    stds = []
    count = 1
    for j in range(0,iterations):
        shuffle(copy_Y)
        maximum = max(copy_Y)
        index_of_closest_value =  min(range(len(copy_Y)), key=lambda i: abs(copy_Y[i]-maximum))
        rel_value = index_of_closest_value + xmin
        std_exp = calculate_std(rel_value,xmin,xmax,copy_Y)
        if std_exp <= obs:
            count = count + 1
        stds.append(std_exp)
    
    return(stds,count)
        
        
            
            
def calculate_std(ref_point,xmin,xmax,array):
    variance = []
    weights = []
    for j in range(xmin,xmax):
        index = j + xmax
        variance.append(array[index]*(j-ref_point)*(j-ref_point))
        weights.append(array[index])
    
    #std = math.pow(sum(variance)/len(variance),0.5)
    std = math.pow(sum(variance)/sum(weights),0.5)
    return(std)

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')
    
    


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Plot_of_mean_vs_std_for_all_factors.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM output directory')
    #parser.add_option('-w', action='store', type='int', dest='window', default = 50,
    #                  help='Window size of moving average., Default=50')
    
    
    (options, args) = parser.parse_args()
    process_file_9(options.output)
    
    
    
if __name__ == "__main__":
    run() 