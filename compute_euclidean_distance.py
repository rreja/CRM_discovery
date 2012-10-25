import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter

def compute_distance(key1,key2,val1,val2,options,):
    coor_list1 = parse_key(key1)
    coor_list2 = parse_key(key2)
    




def parse_key(key1):
    region_list = []
    tmplist1 = key1.split("_")
    tmplist2 = tmplist1[1].split(":")
    region_list = range(int(tmplist2[1]),int(tmplist2[2]),options.window)






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