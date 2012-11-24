import sys, os, operator, random, time
from optparse import OptionParser , IndentedHelpFormatter







        
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
    #parser.add_option('-d', action='store',  dest='distFile',
    #                  help='File containing the distance matrix.')
    parser.add_option('-o', action='store', dest='offFile',
                      help='File containing the offset matrix')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        process_files(args[0],options)
    
if __name__ == "__main__":
    run() 
