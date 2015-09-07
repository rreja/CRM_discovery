import sys, os, re, difflib, csv
from operator import add
from itertools import izip, cycle, tee
from optparse import OptionParser , IndentedHelpFormatter
from pylab import *
import numpy as np
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as plticker
from collections import defaultdict
from scipy.stats.stats import pearsonr
import scipy, pylab
import scipy.cluster.hierarchy as sch


data = {}
def  process_files(options):
    for fname in os.listdir(options.cdtDir):
        if not fname.endswith(".cdt"):
            continue
        label = fname.split("_")[0]
        print "processing "+fname
    
        # Read the cdt and sum all rows
        in0 = open(os.path.join(options.cdtDir,fname),"rt")
        for line in in0:
            if line.startswith("Uniqe"):
                tmp = line.rstrip().split("\t")[2:]
                X = [int(x) for x in tmp]
                xmin = min(X)
                xmax = max(X)
                Y1 = [0]*len(X)
                continue
            if line.startswith("EWEIGHT"):
                continue
            tmplist = line.rstrip().split("\t")[2:]
            newList = [float(x) for x in tmplist]
            Y1 = map(add,Y1,newList)
            
        in0.close()
        data[label] = Y1

    
    # Compute correlation to create NxN matrix
    coeff_matrix = []
    new_label = []
    for k1,v1 in data.items():
        per_row_list = []
        new_label.append(k1)
        for k2,v2 in data.items():
            coeff, pval = pearsonr(v1,v2)
            per_row_list.append(coeff)
        coeff_matrix.append(per_row_list)
    
    outfile = os.path.join(options.cdtDir,"correlation_heatmap.png")
    plot_dendogram(coeff_matrix,new_label,outfile)
            

def plot_dendogram(D,labels,outfile):
    # Convert matrix to numpy matrix
    D = np.matrix(D)
    fig = pylab.figure()
    axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
    Y = sch.linkage(D, method='centroid')
    Z = sch.dendrogram(Y, orientation='right')
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    
    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
    index = Z['leaves']
    D = D[index,:]
    D = D[:,index]
    
    # Convert the newely clustered indexes to new xticklabels.
    xtick_labels = []
    for i in index:
        for idx,vals in enumerate(labels):
            if i == idx:
                xtick_labels.append(vals)
    
    # Plot the matrix.    
    im = axmatrix.matshow(D, aspect='auto', origin='lower')
    
    axmatrix.set_xticks(range(len(labels)))
    axmatrix.set_yticks([])
    
    # Set the labels for xticks
    axmatrix.set_xticklabels(xtick_labels,rotation=45)
    
    
    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
    pylab.colorbar(im, cax=axcolor)
    
    # Display and save figure.
    fig.savefig(outfile,bbox_inches='tight')
    
    
    








usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Create_correlation_heatmap_for_ref.py /usr/local/folder_containing_CDT_files
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-i', action='store', type='string', dest='cdtDir',
                      help='Directory containing CDT files.')
    #parser.add_option('-j', action='store', type='string', dest='CDT2', 
    #                  help='Second CDT file')

    (options, args) = parser.parse_args()
    process_files(options)

    
    
if __name__ == "__main__":
    run() 