import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee
from compute_euclidean_distance import get_fullvectors



def seed_lookup(idxData,filehash,options,mean,std,lookup_regions):
    # Go through 5 iterations of seed lookup
    for i in range(5):
        for locus in lookup_regions:
            # Get all possible 100bp windows and the corresponding tags in those window. Store in vec1 which is dictionary, keys are window-offset, values is a list of 220 bins.
            vec1 = get_fullvectors(locus[0],idxData,filehash,options)
            for k,v in vec1.items():
                sys.exit(1)
                
                
                
    
   

def get_fullvectors(key,idxData,filehash,options):
    region_dict = {}
    taglist = []
    finallist = []
    count = 0
    chrom = key.split(":")[0]
    start = key.split(":")[1]
    end = key.split(":")[2]
    windows = sliceIterator(range(int(start),int(end)),options.windowLength)
    for i in windows:
        for key,val in filehash.items():
            for j in i:
                # idx if formed by fileNo:chr+start, where filenNo is stored in filehash. Each entry in filehash corresponds to one tab/tag file.
                idx = str(key)+":"+chrom+":"+str(j)
                if idx in idxData:
                    taglist.append(idxData[idx])
                else:
                    taglist.append(0)
            taglist = bindata(taglist,options)
            finallist = finallist + taglist
            taglist = []
        region_dict[count] = finallist
        count = count+1
        finallist = []
    return(region_dict)

def bindata(taglist,options):
    bintaglist = []
    summation = 0
    tmplist = range(0,len(taglist)+1,options.bins)
    for elem , next_elem in pairwise(tmplist):
        for j in taglist[elem:next_elem:1]:
            summation = summation+j
        bintaglist.append(summation)
        summation = 0
    return(bintaglist)

def pairwise(seq):
    a, b = tee(seq)
    next(b)
    return izip(a, b)

def sliceIterator(lst, sliceLen):
    for i in range(len(lst) - sliceLen + 1):
        yield lst[i:i + sliceLen]


usage = '''
input_paths may be:
Not to be used as as independent script. To be called by 3_create_seed_motif.py

example usages:
python search_seed_in_enriched_loci.py
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
    
  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-w', action='store', type='int', dest='window',default=100,
    #                  help='Size of the moving window.Default=100')
    #parser.add_option('-b', action='store', type='int', dest='bins',default=5,
    #                  help='Sub bin within a window. This will also be used as a sliding distance.Default=5')
    (options, args) = parser.parse_args()
       
    if not args:
        parser.print_help()
        sys.exit(1)
    seed_lookup(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])   
    #compute_distance(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],options,sys.argv[6])
    
if __name__ == "__main__":
    run() 
