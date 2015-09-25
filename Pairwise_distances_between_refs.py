import sys, os
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import pybedtools
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import numpy as np


def process_file(options):
       
    # Read the directory containing different states
    for direc in os.listdir(options.output):
        direc_path = os.path.join(options.output,direc)
        if not os.path.isdir(direc_path):
            continue
        
        Files = {}
        # Get information about the motif
        Tss_file = os.path.join(direc_path,"enriched_TSS.gff")
        Files["TSS"] = Tss_file
        outfile = os.path.join(direc_path,"pairwise-distance_comparison.png")
        
        for fname in os.listdir(os.path.join(direc_path,"_candidate_refs")):
            if not fname.endswith(".gff"):
                continue
            ref_label = fname.split("-")[0]
            # Create directory for motifs
            
            ref_file = os.path.join(os.path.join(direc_path,"_candidate_refs"),fname)
            Files[ref_label] = ref_file
        
        f,ax = plt.subplots(len(Files.keys()),len(Files.keys()), sharex = 'col')
        
        # Compare pair-wise distances
        count1 = 0
        for label1,file1 in Files.items():
            count2 = 0
            for label2, file2 in Files.items():
                if label1 == label2:
                    ax[count1,count2].set_yticks([])
                    count2 = count2 + 1
                    continue
                X,Y = get_distances(file1,file2,options,direc_path)
                plot_graph(ax,X,Y,count1,count2,label1,label2)
                count2 = count2 + 1
            count1 = count1 + 1
        
        
        plt.tight_layout()
        plt.savefig(outfile)        
        sys.exit(1)


def get_distances(file1,file2,options,outdir):
    tmpout = os.path.join(outdir,"tmp.txt")
    sloped_ref = pybedtools.BedTool(file1).slop(g=options.genFile,s=True,l=500,r=500)
    sloped_ref.intersect(file2,wao=True).saveas(tmpout)
    in0 = open(tmpout,"rt")
    dist = []
    for line in in0:
        cols = line.rstrip().split("\t")
        if cols[9] == ".":
            continue
        mid_point_ref = int((int(cols[3]) + int(cols[4]))/2)
        mid_point_test = int((int(cols[12]) + int(cols[13]))/2)
        if cols[6] == "+":
            dist.append(mid_point_test - mid_point_ref)
        else:
            dist.append((-1)*(mid_point_test - mid_point_ref))
    
    # Calculate X,Y from dist list:
    dicti = {}
    X = []
    Y = []
    
    for k in dist:
        if k in dicti:
            dicti[k] = dicti[k]+1
        else:
            dicti[k] = 1
    
    for i in range(-500,500):
        X.append(i)
        if i in dicti:
            Y.append(dicti[i])
        else:
            Y.append(0)
            
    return(X,Y)
    

def plot_graph(ax,X,Y,count1,count2,xlab,ylab):
    
    #Y = smoothListGaussian(Y,20)
    #X = smoothListGaussian(X,20)
    Y = movingaverage(Y,50)
    X = movingaverage(X,50)
    ax[count1,count2].plot(X,Y,lw=1.0)
    ax[count1,count2].set_xlim(-300,300)
    list1 =  ax[count1,count2].get_yticks().tolist()
    list2 =  ax[count1,count2].get_xticks().tolist()
    ax[count1,count2].set_yticks([min(list1),max(list1)])
    ax[count1,count2].set_xticks([-200,0,200])
    ax[count1,count2].set_title(ylab)
    ax[count1,count2].set_xlabel(xlab)
    

    
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')    

    
def smoothListGaussian(list,degree):  

     window=degree*2-1  
     weight=np.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(np.exp((4*(frac))**2))  
         weightGauss.append(gauss)  
     weight=np.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  
     return smoothed  



usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Pairwise_distances_between_refs.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM output directory')
    parser.add_option('-g', action='store', type='string', dest='genFile',
                      help='File containing chromosome lengths.')
    
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
    
if __name__ == "__main__":
    run() 