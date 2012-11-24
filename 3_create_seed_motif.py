import sys, os, operator, random, time
from optparse import OptionParser , IndentedHelpFormatter

def process_files(infile,options):
    for i in range(1,options.nseed):
        Dmatx = {}
        Omatx = {}
        count = 0
        randList = generate_random_numbers(options.sample,options.max)
        # First reading the distance matrix file and sampling 100 rows.
        input1 = open(infile,"rt")
        for line in input1:
            count = count+1
            if count == 1: # Take the first line which is the heder and the order will be restored in it.
                header = line.rstrip().split("\t")
                header.remove('')
            if count in randList:
                tmp = line.rstrip().split("\t")
                Dmatx[tmp[0]] = map_list_to_list(tmp[1:],header)
            else:
                continue
            if count > max(randList):
                break
        # Now reading the offset matrix file and reading the corresponding rows.
        count2 = 0
        input2 = open(options.offFile,"rt")
        for line in input2:
            count2 = count2+1
            if count2 in randList:
                tmp = line.rstrip().split("\t")
                Omatx[tmp[0]] = map_list_to_list(tmp[1:],header)
            else:
                continue
            if count2 > max(randList):
                break
        
        #for k,v in Omatx.items():
        #    print k+"hello"
        #    for key, val in v.items():
        #        print key,val
        #        sys.exit(1)
        sys.exit(1)
                
                
            
        
    
def map_list_to_list(list1,list2):
    dicti = {}
    for i in range(len(list2)):
        dicti[list2[i]] = list1[i]
    return(dicti)
        


def generate_random_numbers(s,m):
    randlist = []
    randlist = random.sample(range(m), s)
    return(randlist)

 
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
                      help='The number of rows to sample at a time.')
    parser.add_option('-m', action='store', type='int', dest='max',default = 5000,
                      help='Generate random nummbers between 1 and max.')
    parser.add_option('-o', action='store', dest='offFile',
                      help='File containing the offset matrix')
    parser.add_option('-s', action='store', type='int', dest='nseed',default = 50,
                      help='The number of seed motif to generate.')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    process_files(args[0],options)
    
if __name__ == "__main__":
    run() 
