import sys, os, operator, random
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee
from collections import defaultdict
from multiprocessing import Process

from compute_euclidean_distance import compute_distance, get_fullvectors, merge_list, sliceIterator

def process_files(idxdir,options,outfile1,outfile2):
    count = 1
    idxData = {}
    encData = {}
    allData = {}
    filehash = {}
    taglist = []
    out1 = open(outfile1,"w")
    out2 = open(outfile2,"w")
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
    
    print_matrix_header(out1,out2,encData)
    create_matrix_for_regions(encData,allData,filehash,options,out1,out2)               

def create_matrix_for_regions(encData,allData,filehash,options,out1,out2):
    known_dist = {}
    known_offset = {}
    jobs = []
    count = 0
    for majorkey, majorval in encData.items():

        val1 = get_vectors(majorkey,allData,filehash)
        #line1,line2 = run_job_per_line(majorkey,val1,encData,allData,filehash,options)
        p = Process(target=run_job_per_line,args=(majorkey,val1,encData,allData,filehash,options))
        p.start()
        jobs.append(p)
        count = count + 1
        if count == 4:
            for p in jobs:
                p.join()
            sys.exit(1)

def run_job_per_line(majorkey,val1,encData,allData,filehash,options):
    print majorkey
    line1 = majorkey
    line2 = majorkey
    for minorkey, minorval in encData.items():
        val2 = get_vectors(minorkey,allData,filehash)
        offset, dist = compute_distance(majorkey,minorkey,val1,val2,options.window,options.bins,filehash)
        line1 = line1+"\t"+str(dist)
        line2 = line2+"\t"+offset
    #return(line1,line2)
    print line1
    #return(line1)
    
def print_matrix_header(out1,out2,encData):
    line = ""
    for k,v in encData.items():
        line = line+"\t"+k
    out1.write(line+"\n")
    out2.write(line+"\n")
        
    
def get_vectors(key,allData,filehash):
    tmpdict = {}
    for k,v in filehash.items():
        #keystr = str(k)+"_"+key.split("_")[1]
        keystr = str(k)+"_"+key
        tmpdict[keystr] = allData[keystr]
    #return(key.split("_")[1],tmpdict)  # return key in the form of chr:start:end
    return(tmpdict)


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
    outfile1 = os.path.join(os.path.dirname(options.enrcfile),"distance_matrix.txt")
    outfile2 = os.path.join(os.path.dirname(options.enrcfile),"offset_matrix.txt")
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all idx files."
    else:
        process_files(args[0],options,outfile1,outfile2)
    
if __name__ == "__main__":
    run() 
