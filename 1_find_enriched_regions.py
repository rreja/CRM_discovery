import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy as np













usage = '''
input_paths may be:
- a single file.

example usages:
python 1_find_enriched_regions.py
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-p', action='store', type='string', dest='PPdir',
    #                  help='Dir cotaining all the peak-pair file in gff format.')
    parser.add_option('-g', action='store',  dest='reffile',
                      help='File containing chromosome length')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    try:
        import pybedtools
    except ImportError:
        print "You need to install pybedtools."
        sys.exit(1)
     
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all gff files."
    
if __name__ == "__main__":
    run() 