import sys, os, operator, time
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee
from compute_euclidean_distance import get_fullvectors
from math import exp, log
from multiprocessing import Process,Queue
import multiprocessing as mp




def seed_lookup(idxData,filehash,options,mean,std,lookup_regions,bmean,bstd):
    # Go through 5 iterations of seed lookup
    # excluded list includes all regions we do not want to consider either they have no motif or their motif was already found.
    excluded = []
    resample = 0
    # Output_offsets will include all the locus and their offset that had motif. This will be printed in a file to create cluster plot.
    output_offsets = {}
    for iteration in range(5):
        print "Round number: "+str(iteration+1)
        #print mean
        out_q1 = Queue()
        out_q2 = Queue()
        out_q3 = Queue()
        all_window_results = {}
        all_window_vals = {}
        all_window_offsets = {}
        count = 0
        num = (mp.cpu_count()) - 2
        print "Using "+str(num)+" processors in your computer."
        jobs = []
        linecount = 0
        for locus in lookup_regions:
            if locus[0] in excluded:
                continue
            else:
                vec1 = get_fullvectors(locus[0],idxData,filehash,options)
                p = Process(target=run_process_in_parallel,args=(vec1,locus[0],options,mean,std,bmean,bstd,out_q1,out_q2,out_q3))
                jobs.append(p)
                p.start()
                count = count + 1
                if count == num:
                    for i in range(num):
                        all_window_results.update(out_q1.get())
                        all_window_vals.update(out_q2.get())
                        all_window_offsets.update(out_q3.get())
                    for p in jobs:
                        p.join()
                    count = 0
                
        # Accouting for the jobs that were running but the loop had to exit.
        for j in range(count):
            all_window_results.update(out_q1.get())
            all_window_vals.update(out_q2.get())
            all_window_offsets.update(out_q3.get())
        for p in jobs:
            p.join()
        
        all_window_vals = get_top_n(all_window_results,all_window_vals,iteration)
            
        if iteration == 0:
            if len(all_window_vals.keys()) <= 30:
                print "Not a good seed. Less than 30 matches found."
                # Exclude this locus from the further analysis.
                excluded.append(locus[0])
                # variable resample to send a signal to 3_seed...py script to ignore the current seed motif and resample the next one.
                resample = 1
                # Also exclude all the locus that matched with a higher likelihood with this locus.
                for k,v in all_window_vals.items():
                    excluded.append(k)
                break
                
        if len(all_window_vals.keys()) == 0:
            break
        else:
            mean,std,output_offsets,excluded = get_mean_and_std(all_window_vals,all_window_offsets,excluded,output_offsets)
        

    return(resample,excluded,output_offsets,mean)
    
def get_mean_and_std(dicti,offsets_dicti,excluded,output_offsets):
    summed_list = []
    count = 0
    #output_offsets = {}
    n = len(dicti.keys())
    print "Locus with high Likelihood= "+str(n)
    # Calculating mean.
    for k,v in dicti.items():
        # updating offsets that will be printed.
        output_offsets[k] = offsets_dicti[k]
        # updating excluded regions.
        excluded.append(k)
        if count == 0:
            summed_list = v
        else:
            for j in range(len(v)):
                summed_list[j] = summed_list[j] + v[j]
        count = count + 1
    # Adding psudocount of 0.5 to the mean
    newmeanList = [(float(x) + 0.5)/(n+1) for x in summed_list]
    # Calculatind std here.
    sum_of_squares = [0]*len(newmeanList)
    tmpnewList = []
    for k,v in dicti.items():
        for m in range(len(newmeanList)):
            sum_of_squares[m] = sum_of_squares[m] + ((v[m] - newmeanList[m])*(v[m] - newmeanList[m]))
    for i in range(len(newmeanList)):
        ##Adding psudo-count of 0.5 to sigma.
        tmpnewList.append((sum_of_squares[i] + ((0.5-newmeanList[i])**2))/n)
    
    newstdList = [(x**0.5) for x in tmpnewList]
    # Returning all the computed values.
    return(newmeanList,newstdList,output_offsets,excluded)
            
def get_top_n(all_window_results,all_windows_val,iteration):
    # tmpdict will store dict with keys as regions and values as likelihood score.
    tmpdict = {}
    sorted_result =  sorted(all_window_results.iteritems(), key=operator.itemgetter(1),reverse=True)
    # take top 50, for the remaining take whatever has L > 0.
    if iteration == 0:
        for i in sorted_result[:50]:
            if all_window_results[i[0]] > 0:
                tmpdict[i[0]] = all_windows_val[i[0]]
    else:
        for i in sorted_result:
            if all_window_results[i[0]] > 0:
                tmpdict[i[0]] = all_windows_val[i[0]]
    return(tmpdict)
        
                
def run_process_in_parallel(vec1,locus,options,mean,std,bmean,bstd,out_q1,out_q2,out_q3):
    # Store in this dict one window per locus that had the highest liklihood. the key will be chr:start:end:offset, value will be the liklihood score.
    all_aligned_windows = {}
    # Store in this dict the value vector for the highest scoring window.
    all_aligned_windows_val = {}
    # Storing all the aligned offsets in the dictionary, with key as locus and value as offset.
    all_aligned_offsets = {}
    # Initializing highest score.
    highestScore = -1
    highest_match_offset = 0
    highest_match_val = []
    for k,v in vec1.items():
        L = calculate_likelihood(mean,std,bmean,bstd,v,options)
        if highestScore > L:
            continue
        else:
            highestScore = L
            highest_match_offset = k
            highest_match_val = v
            
    all_aligned_windows[locus] = L
    all_aligned_windows_val[locus] = highest_match_val
    all_aligned_offsets[locus] = highest_match_offset
    
    
    out_q1.put(all_aligned_windows)
    out_q2.put(all_aligned_windows_val)
    out_q3.put(all_aligned_offsets)
                
    
def calculate_likelihood(mean,std,bmean,bstd,v,options):
    # Declare proiors, PA=P2B = 0.49 and P2A = 0.01
    bins_per_factor = options.windowLength/options.bins
    pa = 0.49
    p2b = 0.49
    p2a = 0.01
    likelihood = 0
    sum_L = 0
    for i in range(0,len(v),bins_per_factor):
        M = 1
        B = 1
        A = 1 
        for j in range(i,i+20):
            M = M*(exp(-((v[j] - mean[j])**2)/(2*(std[j]**2))))
            B = B*(exp(-((v[j] - bmean[j])**2)/(2*(bstd[j]**2))))
            A = A*(exp(-((1.75 - 0)**2)/(2*(1**2))))

        try:
            NUM = log(M*pa)
        except ValueError:
            NUM = 0
        try:
            DEN = log((B*p2a)+(A*p2b))
        except ValueError:
            DEN = 0
            
        likelihood = NUM - DEN      
        sum_L = sum_L + likelihood
    return(sum_L)
        
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
    #Returnng the region_dict with keys as window offset and values as bined tags.
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
