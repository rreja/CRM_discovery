import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy


def compute_distance(key1,key2,val1,val2,window,bins,filehash):
    #distance = []
    distance = {}
    binned_window_length = window/bins
    if binned_window_length != 0:
        vec1 = get_vectors(key1,val1,binned_window_length,filehash)
        vec2 = get_vectors(key2,val2,binned_window_length,filehash)
        for ak,av in vec1.items():
            for bk,bv in vec2.items():
                distance[str(ak)+":"+str(bk)] = numpy.linalg.norm(numpy.array(av)-numpy.array(bv))
                #print str(ak)+":"+str(bk), distance[str(ak)+":"+str(bk)]
                
    else:
        print "Your binned window length became zero!"
    sorted_list = sorted(distance, key=distance.get, reverse=True)
    return(sorted_list[0], distance[sorted_list[0]])
    

    
def get_vectors(key,val,window,filehash):
    region_list = []
    region_dict = {}
    count = 0
    for k,v in filehash.items():
        idx = str(k)+"_"+key
        tmp = sliceIterator(val[idx],window)
        if count == 0:
            region_list = tmp
        else:
            region_list = merge_list(region_list,tmp)
        count = count+1
    key = 0    
    for j in region_list:
        region_dict[key] = j
        key = key+1
    #return(region_list)
    return(region_dict)
          
    

def merge_list(list1,list2):
    list3 = []
    for i,j in map(None,list1,list2):
        list3.append(i+j)
    return(list3)

def sliceIterator(lst, sliceLen):
    for i in range(len(lst) - sliceLen + 1):
        yield lst[i:i + sliceLen]


usage = '''
input_paths may be:
- a single file.

example usages:
python compute_euclidean_distance.py 
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
    
  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='window',default=100,
                      help='Size of the moving window.Default=100')
    parser.add_option('-b', action='store', type='int', dest='bins',default=5,
                      help='Sub bin within a window. This will also be used as a sliding distance.Default=5')
    (options, args) = parser.parse_args()
       
    if not args:
        parser.print_help()
        sys.exit(1)
        
    compute_distance(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],options,sys.argv[6])
    
if __name__ == "__main__":
    run() 