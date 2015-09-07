import sys, os
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
from scipy.stats.stats import pearsonr
import scipy, pylab
import scipy.cluster.hierarchy as sch

def process_file(options,outfile,outtmp):
    print "Calculating the Pearson correlation between all the columns."
    data = defaultdict(list)
    labels = []
    input = open(options.matrixFile,"rt")
    for line in input:
        if line.startswith("state") or line.startswith("genes"):
            #labels = line.rstrip().split("\t")[1:]
            continue
        cols = line.rstrip().split("\t")
        data[int(cols[0])] = map(float,cols[1:])
        labels.append(cols[0])
        #no_of_factors = len(map(float,cols[1:]))
    
    
    no_of_genes = len(data.keys())
    coeff_matrix = []
    for i in range(1,no_of_genes+1):
        per_row_list = []
        for j in range(1,no_of_genes+1):
            list1 = []
            list2 = []
            list1 = data[i]
            list2 = data[j]
            #print list1, list2
            coeff, pval = pearsonr(list1,list2)
            per_row_list.append(coeff)
        coeff_matrix.append(per_row_list)
    
    #In case you want to write down the correlation matrix to a file uncomment this region.
    #count = 0
    #out = open(outtmp,"w")
    ##out.write("ID\tRap1\tIfh1\tSfp1\tHmo1\tFhl1\tTup1\tAbf1\tGcn4\tSua7\n")
    #for l in coeff_matrix:
    #    line = labels[count]
    #    count = count + 1
    #    for m in l:
    #        line = line+"\t"+str(m)
    #    out.write(line+"\n")
    
    print "Clustering the correlation matrix and plotting it."
    plot_dendogram(coeff_matrix,labels,outfile)
        
    


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
    #print D
    im = axmatrix.matshow(D, aspect='auto', origin='lower')
    
    axmatrix.set_xticks(range(len(labels)))
    axmatrix.set_yticks([])
    
    # Set the labels for xticks
    axmatrix.set_xticklabels(xtick_labels)
    
    
    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
    pylab.colorbar(im, cax=axcolor)
    
    # Display and save figure.
    fig.savefig(outfile,bbox_inches='tight')



usage = '''
input_paths may be:
- The matrix file only.

############################################
# IMPORTANT NOTES:
# 1) The colorbar might appear black for certain cases. It's not due to the script,
it is due to the matplotlib versions. Try it on Rohit's Mac book pro, works fine
there.
# 2) The input file should be in the following format:
ID      Factor1     Factor2 ....
Gene1   120         230
Gene2   12          20
...
############################################
example usages:
python compute_correlation_and_plot_heatmap.py -m  <matrix_file_in_the_above_mentioned_format.>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
   
    #parser.add_option('--consensus', action='store', type='string', dest='cons',default='G[AT]ATA',
    #                  help = 'Enter your consensus sequence from the MEME motif, multiple bases at one position are enclosed by []')
    parser.add_option('-m', action='store', type='string', dest='matrixFile',
                      help = 'Data matrix for all the factors.')
    #parser.add_option('-i', action='store', type='string', dest='idxDir',
    #                  help = 'Directory containing idx files.')
    #parser.add_option('-g', action='store', type='string', dest='genFile',
    #                  help='Chromosome length file')
   
    (options, args) = parser.parse_args()
    

    outfile = os.path.join(os.path.dirname(options.matrixFile),"heatmap_from_matrix.svg")
    outtmp = os.path.join(os.path.dirname(options.matrixFile),"heatmap_from_matrix.txt")
    process_file(options,outfile,outtmp)
    
      
    
if __name__ == "__main__":
    run() 