from optparse import OptionParser , IndentedHelpFormatter
import sys, os, re, operator, math
from collections import OrderedDict
from operator import add
import numpy as np
from pylab import *


def  process_file(cdtDir,options):
    xdata = []
    ydata = []
    z = []
    #data = defaultdict(list)
    outfile = os.path.join(cdtDir,"Means_vs_std_plot.svg")
    for fname in os.listdir(cdtDir):
        data = {}
        if not (fname.endswith('.cdt') or fname.endswith('.txt')):
            continue
        
        in1 = open(os.path.join(cdtDir,fname))
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
        
        Y = movingaverage(Y,options.window)
        X = movingaverage(X,options.window)
        X = [int(x) for x in X]
        xmin = min(X)
        xmax = max(X)
        print xmin,xmax
        Y = [float(x)/max(Y) for x in Y]
        # Get the peak
        maximum = max(Y)
        # Find the corresponding index associated with the peak
        index_of_closest_value =  min(range(len(Y)), key=lambda i: abs(Y[i]-maximum))
        # Change the index to relative X-coordinate value.
        rel_value = index_of_closest_value + xmin
        std = calculate_std(rel_value,xmin,xmax,Y)
        #print rel_value,std,label
        xdata.append(rel_value)
        ydata.append(std)
        z.append(label)
        
    
    # Scatter plots
    colors = np.random.rand(50)
    fig, ax = plt.subplots()
    ax.scatter(xdata, ydata,c=colors,alpha=0.5)
    for i, txt in enumerate(z):
        ax.annotate(txt, (xdata[i],ydata[i]))
    savefig(outfile)
        
        
        
            
            
def calculate_std(ref_point,xmin,xmax,array):
    variance = []
    for j in range(xmin,xmax):
        index = j + xmax
        variance.append(array[index]*(j-ref_point)*(j-ref_point))
    
    std = math.pow(sum(variance)/len(variance),0.5)
    return(std)

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')
    
    


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python sort_CDT_by_given_file.py  <path_to_CDT_folder>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='window', default = 50,
                      help='Window size of moving average., Default=5')
    #parser.add_option('-i', action='store', type='string', dest='sortbyFile',
    #                  help='Input file which will decide the sorting order')
    #parser.add_option('-o', action='store', type='int', dest='criteria',default=0,
    #                  help='Sorting criteria, 0=> sort by first column of any file, 1=> sort by given chr:start-end order in the file, 2 => last column of gff')
    #parser.add_option('-n', action='store', type='int', dest='extra',default=0,
    #                  help='Extra N lines to add to the end of CDT')
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    #outdir = os.path.join(os.path.dirname(args[0]),"sorted")
    #if not os.path.exists(outdir): os.makedirs(outdir)
    
     
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    
    process_file(args[0],options)
    
    
    
if __name__ == "__main__":
    run() 