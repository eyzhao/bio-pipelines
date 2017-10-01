''' Usage: clean_germline.py -i INPUTPATH [ -o OUTPUTPATH ]

Options:
    -i INPUTPATH    Path to input germline txt file
    -o OUTPUTPATH   Path to output cleaned germline txt file [ default: None ]
'''

from docopt import docopt

args = docopt(__doc__)

input = open(args['-i']).read()
out_str = input.replace('\nCOSM', '; COSM').replace('\nrs', '; rs')
out_path = args['-o'] if args['-o'] is not None else args['-i'] + '.clean'
output = open(out_path, 'w')
output.write(out_str)
output.close()
