import sys, os, operator, random, time
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee
from search_seed_in_enriched_loci import seed_lookup


def process_files(idxData,filehash,options,bmean,bstd):
    # The main loop to create "n" seed motifs
    outfile = os.path.join(os.path.dirname(options.distFile),"seed_motif_profiles.txt")
    out = open(outfile,"w")
    for i in range(options.nseed):
        print "Running the "+str(i+1)+" round"
        Dmatx = {}
        Omatx = {}
        count = 0
        highestScore = 0
        avgprofile = []
        standardDeviation = []
        seedKey = ""
        randList = generate_random_numbers(options.sample,options.max)
        #print randList
        # First reading the distance matrix file and sampling 100 rows.
        input1 = open(options.distFile,"rt")
        for line in input1:
            count = count+1
            if count == 1: # Take the first line which is the heder and the order will be restored in it.
                header = line.rstrip().split("\t")
                header.remove('')
            if count in randList:
                tmp = line.rstrip().split("\t")
                Dmatx[tmp[0]] = map_list_to_list(tmp[1:],header,"dist")
            else:
                continue
            # If you have read all the lines in randList, just stop reading the file.
            if count > max(randList):
                break
        # Now reading the offset matrix file and reading the corresponding rows.
        count2 = 0
        input2 = open(options.offFile,"rt")
        for line in input2:
            count2 = count2+1
            if count2 in randList:
                tmp = line.rstrip().split("\t")
                Omatx[tmp[0]] = map_list_to_list(tmp[1:],header,"offset")
            else:
                continue
            if count2 > max(randList):
                break
        # Now for each of the 100 sampled rows, find the top 20 distances and their corresponding offsets.
        for key,val in Dmatx.items():
            #print key
            sorted_val = sorted(val.iteritems(), key=operator.itemgetter(1),reverse=False)
            ## Getting offset for the same key.
            ## Taking top 100 closest loci and looking which "key" window was found overlapping the most.
            offsets = get_offsets(Omatx[key],sorted_val[:100]) # Change the number here to change the 'top x' selection.
            # Get all the enriched windows or continue with the loop if you do not file a row that has a window represented atleast 20 times.
            encWindows = find_most_enriched_window(offsets)
            # iterate through all the enriched windows and create average profile
            if not encWindows == 0:
                ##Binned vectors corresponding to all the enriched windows
                allVectors = create_data_vectors(encWindows,idxData,filehash,options)
                meanVector = get_mean(allVectors,len(allVectors.keys()))
                #print meanVector
                ## Get a standard deviation vector for each bin.
                std  = get_std(allVectors,meanVector,len(allVectors.keys()))
                ## Using the mean and std find out the seed score.
                seedscore = get_seed_score(meanVector,std,options)
                # Sorting in the loop to find out the best seed score and the corresponding average profile.
                if highestScore > seedscore:
                    continue
                else:
                    highestScore = seedscore
                    avgprofile = meanVector
                    standardDeviation = std
                    seedKey = key
            else:
                continue
            ## Run this command from CRM_discovey folder:  time python 3_create_seed_motif.py -m 10 -n 3 -o testdata/output/tmpoffset.txt -d testdata/output/tmpdist.txt ../tags/ -e testdata/output/tmp.enriched.gff 
        ## Write the average seed profiles to a file.
        # Instead of giving all regions for seed lookup, we can give the top 1000/5000 regions from the same row that created the seed motif.
        lookup_regions = sorted(Dmatx[seedKey].iteritems(), key=operator.itemgetter(1),reverse=False)
        # Function imported from another script.
        seed_lookup(idxData,filehash,options,avgprofile,standardDeviation,lookup_regions[:5000],bmean,bstd)
        # Remove this once the script is complete.
        sys.exit(1)
        #Commenting out the section below since we may not be writing the seeds to file anymore.
        #line = seedKey
        #for z in avgprofile:
        #    line = line+"\t"+str(z)
        #out.write(line+"\n")
        #out.flush()
        

def get_seed_score(mean,std,options):
    # First find out how many bins for each factor to consider.
    # len(mean) represents the total bins present.
    windows_per_factor = options.windowLength/options.bins
    half_window = windows_per_factor/2
    seed_score = 0
    for i in range(0,len(mean),windows_per_factor):
        list1 = sorted(mean[i:(i+windows_per_factor)],reverse=True)
        list2 = std[i:(i+windows_per_factor)]
        seed_score = seed_score + seed_score_calculation(list1,list2,half_window,windows_per_factor)
    return(seed_score)
    
def seed_score_calculation(mean,std,half_window,full_window):
    first_half = 0
    second_half = 0
    sum_std = 0
    seedscore = 0
    for j in range(0,half_window):
        first_half = first_half +mean[j]
        sum_std = sum_std + std[j]
    for k in range(half_window,full_window):
        second_half = second_half + mean[k]
        sum_std = sum_std + std[k]
    std_per_factor = float(sum_std)/full_window
    seedscore = float(abs(first_half) - abs(second_half))/std_per_factor
    #print mean
    #print first_half,second_half,std_per_factor,seedscore
    #sys.exit(1)
    return(seedscore)
    
def create_data_vectors(encwindows,idxdata,filehash,options):
    taglist = []
    finallist = []
    allData = {}
    #iterate over all enriched windows and get their tag.
    for k,v in encwindows.items():
        # Get the coordinates for actual window in that locus.
        new_k = get_actual_cord(k,(v.split(":"))[1],options)
        l = new_k.split(":")
        # Iterate over all tag/tab files one by one
        for key,val in filehash.items():
            for i in range(int(l[1]),int(l[2])):
                # idx if formed by fileNo:chr+start, where filenNo is stored in filehash. Each entry in filehash corresponds to one tab/tag file.
                idx = str(key)+":"+l[0]+":"+str(i)
                if idx in idxdata:
                    taglist.append(idxdata[idx])
                else:
                    taglist.append(0)
            taglist = bindata(taglist,options)
            finallist = finallist+taglist
            taglist = []
        ##allData[k] = bindata(taglist,options)
        ##taglist = []
        allData[k] = finallist
        finallist = []
        
    return(allData)   
    
def read_tag_files(idxdir,options):
    idxData = {}
    filehash = {}
    count = 1
    # Reading all idx files in the directory one by one
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                chrom, start,ftag,rtag,ttag = line.rstrip().split("\t")
                # Appending the count string in the key to recognize each tag file seperately in the same idxData dict
                idxData[str(count)+":"+chrom+":"+start] = int(ttag)
            filehash[count] = fname
            count = count+1
    #for a,b in filehash.items():
    #    print a,b
    #sys.exit(1)
    return(idxData,filehash)            
    
def find_most_enriched_window(offsets):
    tmpdict = {}
    returndict = {}
    for k,v in offsets.items():
        tmp = v.split(":")
        if tmp[0] in tmpdict:
            tmpdict[tmp[0]] = tmpdict[tmp[0]]+1
        else:
            tmpdict[tmp[0]] = 1
    sorted_tmpdict = sorted(tmpdict.iteritems(), key=operator.itemgetter(1),reverse=True)
    top_window =  int(sorted_tmpdict[0][0])
    # From the dictionary sepeate the ones containing the top window and put it in returndict
    for k,v in offsets.items():
        tmp = v.split(":")
        if int(tmp[0]) == top_window:
            returndict[k] = v
    # So return the dictionary when you find atleast 20 enriched windows among 100.
    if len(returndict.keys()) < 20:
    #if len(returndict.keys()) <= 1:
        return(0)
    else:
        return(returndict)
            
def get_actual_cord(key,windowNo,options):
    tmp = key.split(":")
    start = int(tmp[1]) + options.bins*int(windowNo)
    end = start + options.windowLength
    return(tmp[0]+":"+str(start)+":"+str(end))
         
def get_offsets(dicti,tuple1):
    offsetDict = {}
    list1 = []
    for i in range(len(tuple1)):
        list1.append(tuple1[i][0])
    for k,v in dicti.items():
        if k in list1:
            offsetDict[k] = v
    return(offsetDict)
              
def map_list_to_list(list1,list2,param):
    dicti = {}
    if param == "dist":
        for i in range(len(list2)):
            dicti[list2[i]] = int(list1[i])
    elif param == "offset":
        for i in range(len(list2)):
            dicti[list2[i]] = list1[i]
    return(dicti)
        
def generate_random_numbers(s,m):
    randlist = []
    randlist = random.sample(range(2,m), s)
    return(randlist)

def get_mean(dicti,n):
    summed_list = []
    count = 0
    for k,v in dicti.items():
        if count == 0:
            summed_list = v
        else:
            for j in range(len(v)):
                summed_list[j] = summed_list[j] + v[j]
        count = count + 1
    # Adding psudocount of 0.5 to the mean
    newList = [(float(x) + 0.5)/(n+1) for x in summed_list]
    return(newList)

def get_std(dicti,meanList,n):
    sum_of_squares = [0]*len(meanList)
    tmpnewList = []
    for k,v in dicti.items():
        for m in range(len(meanList)):
            sum_of_squares[m] = sum_of_squares[m] + ((v[m] - meanList[m])*(v[m] - meanList[m]))
    for i in range(len(meanList)):
        ##Adding psudo-count of 0.5 to sigma.
        tmpnewList.append((sum_of_squares[i] + ((0.5-meanList[i])**2))/n)
    
    ##tmpnewList2 = [float(x)/(n-1) for x in sum_of_squares]
    ##newList2 = [(x**0.5) for x in tmpnewList2]
    ##print newList2
    newList = [(x**0.5) for x in tmpnewList]
    return(newList)
                
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

def calculate_background_mean_std(idxData,filehash,options):
    input = open(options.enrFile,"rt")
    sum_total = 0
    allData = {}
    bmean  = {}
    taglist = []
    bstd = {}
    mean = []
    std = []
    bins_per_factor = options.windowLength/options.bins
    # Initialization
    for k,v in filehash.items():
        allData[k] = []
    for line in input:
        if line.startswith("chrom") or line.startswith("#"):
            continue;
        chrom,junk,junk,start,end,junk,junk,junk,junk = line.rstrip().split("\t")
        # Iterate over all tag/tab files one by one
        for key,val in filehash.items():
            for i in range(int(start),int(end)):
                # idx if formed by fileNo:chr+start, where filenNo is stored in filehash. Each entry in filehash corresponds to one tab/tag file.
                idx = str(key)+":"+chrom+":"+str(i)
                if idx in idxData:
                    taglist.append(idxData[idx])
                else:
                    taglist.append(0)
            taglist = bindata(taglist,options)
            allData[key] = allData[key] + taglist
            taglist = []
    for k,v in allData.items():
        ss = 0
        length = len(v)
        summation = sum(v)
        bmean[k] = float(summation)/length
        for val in v:
            ss = ss + ((val - bmean[k])*(val - bmean[k]))
        bstd[k]  =  (float(ss)/(length-1))**(0.5)
        #print bmean[k],bstd[k]
    for k,v in filehash.items():
        mean = mean + [bmean[k]]*bins_per_factor
        std = std + [bstd[k]]*bins_per_factor
    return(mean, std)
        
        

usage = '''
input_paths may be:
- a single file.

example usages:
python 3_create_seed_motif.py -o <offset_matrix.txt>
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-n', action='store', type='int', dest='sample',default = 100,
                      help='[OPTIONAL]: The number of rows to sample at a time.')
    parser.add_option('-e', action='store', dest='enrFile',
                      help='[MANDATORY]: The file containing the enriched regions.')
    parser.add_option('-m', action='store', type='int', dest='max',default = 5000,
                      help='[OPTIONAL]: Generate random nummbers between 1 and max.')
    parser.add_option('-o', action='store', dest='offFile',
                      help='[MANDATORY]: File containing the offset matrix')
    parser.add_option('-d', action='store', dest='distFile',
                      help='[MANDATORY]: File containing the distance matrix')
    parser.add_option('-s', action='store', type='int', dest='nseed',default = 50,
                      help='[OPTIONAL]: The number of seed motif to generate.')
    parser.add_option('-w', action='store', type='int', dest='windowLength',default=100,
                      help='[OPTIONAL]: Size of the moving window.Should be same as the one used in previous step.Default=100')
    parser.add_option('-b', action='store', type='int', dest='bins',default=5,
                      help='[OPTIONAL]: Sub bin within a window. Should be same as the one used in previous step.Default=5')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    idxData,filehash = read_tag_files(args[0],options)
    bmean,bstd = calculate_background_mean_std(idxData,filehash,options)
    process_files(idxData,filehash,options,bmean,bstd)
    print "All iterations completed. Check your results."
    
if __name__ == "__main__":
    run() 
