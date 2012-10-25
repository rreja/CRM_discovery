import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee


def process_files(idxdir,options):
    count = 1
    idxData = {}
    encData = {}
    allData = {}
    filehash = {}
    taglist = []
    enriched = open(options.enrcfile,"r")
    for line in enriched:
        chrom,junk,junk,start,end,junk,junk,junk,junk  = line.rstrip().split("\t")
        encData[chrom+":"+start+":"+end] = 1
    
    # Reading all idx files in the directory one by one
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                chrom, start,ftag,rtag,ttag = line.rstrip().split("\t")
                idxData[chrom+":"+start] = int(ttag)
                
        # Calculating the tags present in the enrcihed regions.
        for key, val in encData.items():
            l = key.split(":")
            for i in range(int(l[1]),int(l[2])+1):         
                idx = l[0]+":"+str(i)
                if idx in idxData:
                    taglist.append(idxData[idx])
                else:
                    taglist.append(0)
            allData[str(count)+"_"+key] = bindata(taglist,options)
            taglist = []
        filehash[count] = fname
        count = count+1
    
    create_matrix_for_regions(allData,filehash)               

def create_matrix_for_regions(allData,filehash):
    for majorkey, majorval in allData.items():
        key1,val1 = get_vectors(majorkey,allData,filehash)
        for minorkey, minorval in allData.items():
            key2,val2 = get_vectors(minorkey,allData,filehash)
            #print key1,key2,val1, val2
            #sys.exit(1)
            compute_distance(key1,key2,val1,val2)
            

def get_vectors(key,allData,filehash):
    tmpdict = {}
    for k,v in filehash.items():
        keystr = str(k)+"_"+key.split("_")[1]
        tmpdict[keystr] = allData[keystr]
    return(key.split("_")[1],tmpdict)


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
    
def print_dict(dictionary):
    for key, val in dictionary.items():
        print key,val


def pairwise(seq):
    a, b = tee(seq)
    next(b)
    return izip(a, b)

usage = '''
input_paths may be:
- a single file.

example usages:
python 2_create_distance_matrix.py [options] /dir_to_idx_files/
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-e', action='store',  dest='enrcfile',
                      help='File containing the enriched regions.')
    parser.add_option('-w', action='store', type='int', dest='window',default=100,
                      help='Size of the moving window.Default=100')
    parser.add_option('-b', action='store', type='int', dest='bins',default=5,
                      help='Sub bin within a window. This will also be used as a sliding distance.Default=5')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all idx files."
    else:
        process_files(args[0],options)
    
if __name__ == "__main__":
    run() 