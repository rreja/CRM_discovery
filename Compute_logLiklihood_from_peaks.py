import sys, os, math, pybedtools
from optparse import OptionParser , IndentedHelpFormatter
from operator import itemgetter
from collections import defaultdict
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.mlab as mlab
from matplotlib import pyplot as plt
from scipy.stats.stats import pearsonr,spearmanr
import numpy as np
import statsmodels.formula.api as sm 




def process_file(options):
    
    LogLik_TSS = 0
    # Slop the TSS file.
    sloped = pybedtools.BedTool(options.ref).slop(g=options.genFile,s=True,l=options.dist,r=options.dist)
    # Read through the peak directory.
    for fname in os.listdir(options.peakDir):
        if not fname.endswith(".gff"):
            continue
        # Get the label for the file.
        print "Comparing main file =  "+fname
        
        # Create intersection output file.
        tmpIntersect = os.path.join(options.peakDir,"tmpIntersect.txt")
        # Take intersection with TSS
        sloped.intersect(os.path.join(options.peakDir,fname),wao=True).saveas(tmpIntersect)
        X_distances = get_distances(tmpIntersect)
        os.system("rm "+tmpIntersect)
        
        for fname2 in os.listdir(options.peakDir):
            if not fname2.endswith(".gff"):
                continue
            if fname2 == fname:
                continue
            
            print "with "+fname2
            newtmpIntersect = os.path.join(options.peakDir,"newtmpIntersect.txt")
            sloped.intersect(os.path.join(options.peakDir,fname2),wao=True).saveas(newtmpIntersect)
            Y_distances = get_distances(newtmpIntersect)
            os.system("rm "+newtmpIntersect)
            
            # Get X,Y values from defaultdictionary items.
            X,Y = convert_to_list(X_distances,Y_distances)
            result = sm.OLS(Y,X).fit()
            LogLik_TSS = LogLik_TSS + result.llf
            
            
    print "Lok liklihood value wrt TSS = "+str(LogLik_TSS)
    
    #### Now computing the log liklihood for each factor ####
    
                
    
    for fname in os.listdir(options.peakDir):
        if not fname.endswith(".gff"):
            continue
        print ""
        print "Computing the log liklihood for =  "+fname
        logL = 0
        #sloped = os.path.join(options.peakDir,"sloped.txt")
        sloped_peak = pybedtools.BedTool(os.path.join(options.peakDir,fname)).slop(g=options.genFile,b=options.dist)
        #pybedtools.BedTool(os.path.join(options.peakDir,fname)).slop(g=options.genFile,b=options.dist).saveas(sloped)
        
        for fname2 in os.listdir(options.peakDir):
            if not fname2.endswith(".gff"):
                continue
            if fname2 == fname:
                continue
            
            print "first file= "+fname2
            newtmpIntersect2 = os.path.join(options.peakDir,"newtmpIntersect2.txt")
            sloped_peak.intersect(os.path.join(options.peakDir,fname2),wao=True).saveas(newtmpIntersect2)
            X_distances = get_distances(newtmpIntersect2)
            os.system("rm "+newtmpIntersect2)
            
            for fname3 in os.listdir(options.peakDir):
                if not fname3.endswith(".gff"):
                    continue
                if (fname3 == fname2 or fname3 == fname):
                    continue
                
                print "Processing "+fname3
                newtmpIntersect3 = os.path.join(options.peakDir,"newtmpIntersect3.txt")
                sloped_peak.intersect(os.path.join(options.peakDir,fname3),wao=True).saveas(newtmpIntersect3)
                Y_distances = get_distances(newtmpIntersect3)
                os.system("rm "+newtmpIntersect3)
                
                # Get X,Y values from defaultdictionary items.
                X,Y = convert_to_list(X_distances,Y_distances)
                result = sm.OLS(Y,X).fit()
                logL = logL + result.llf
        
        print "Logliklihood value for "+fname+" = "+str(logL)




    
        
def convert_to_list(X_distances,Y_distances):
    X = []
    Y = []
    for k,v in X_distances.items():
        for m in v:
            for n in Y_distances[k]:
                X.append(m)
                Y.append(n)
                
                
    return(X,Y)

def get_distances(infile):
    distances = defaultdict(list)
    in0 = open(infile,"rt")
    for line in in0:
        cols = line.rstrip().split("\t")
        original_TSS = int(cols[3]) + 500
        if cols[9] == ".":
            distances[cols[8]].append(0)
        if (cols[6] == "+" or cols[6] == "."):
            distances[cols[8]].append(int(cols[12]) - original_TSS)
        elif cols[6] == "-":
            distances[cols[8]].append((-1)*(int(cols[12]) - original_TSS))
    
    return(distances)
        
        
    
   
    
    


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python plot_scatter_plot.py --f2=filename.txt
'''.lstrip()

 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-r', action='store', type='string', dest='ref',
                      help='Reference coordinate, for instance TSS/+1 nuc/TES, UNIQUE LAST COLUMN REQUIRED.')
    parser.add_option('-i', action='store', type='string', dest='peakDir',
                      help='Directory containing peak calls (UNIQUE LAST COLUMN REQUIRED) in gff format.')
    parser.add_option('-d', action='store', type='int', dest='dist',default=500,
                      help='Distance to go from TSS/peaks to look for another peak.')
    parser.add_option('-g', action='store', type='string', dest='genFile',
                      help='Chromosome length file')
    
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
    

    
if __name__ == "__main__":
    run() 