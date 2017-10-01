''' Mirrors POG files into a local directory structure.

    August 19, 2017
    Written by Eric Zhao, BC Genome Sciences Centre, ezhao@bcgsc.ca

    Usage: mirror_files.py -i INPUT -r ROOT -p PREFIX [ -n -c COLNAME ]

    Options:
        -i --input INPUT        Path to input TSV, containing a column which holds the paths.
                                    that has the paths to mirror.

        -r --root ROOT          Root of input paths to remove and replace with PREFIX

        -p --prefix PREFIX      Prefix to mirror files to

        -n --dry-run            Only print the mirroring paths. Don't actually copy.

        -c --col-name COLNAME   Name of the column containing paths. If not provided, then reads
                                    data not as a TSV but as a list of paths, one per line with
                                    no header.
'''

from docopt import docopt
import pandas as pd
import os
import re
import shutil
from pprint import pprint

class decorate:
    BOLD = '\033[1m'
    END = '\033[0m'

def mirror_files_from_tsv(input_path, root, prefix, dry_run=False, col_name=None):
    if not root.endswith('/'):
        root = root + '/'

    if not prefix.endswith('/'):
        prefix = prefix + '/'

    if col_name is None:
        input_tsv = pd.read_csv(input_path, sep='\t', header=None, names=['path'])
        paths = input_tsv['path'].tolist()
    else:
        input_tsv = pd.read_csv(input_path, sep='\t')
        paths = input_tsv[col_name].tolist()

    paths_matching = [path for path in paths if path.startswith(root)]
    paths_nonmatching = [path for path in paths if not path.startswith(root)]
    if paths_nonmatching:
        print('Warning: the following paths did not start with the specified root and will be ignored')
        pprint(paths_nonmatching)

    mirror_paths_no_root = [path[len(root):] for path in paths_matching]
    mirror_paths = [prefix + path for path in mirror_paths_no_root]

    paths_mapping = zip(paths_matching, mirror_paths)
    for old, new in paths_mapping:
        print('{0} {1}==>{2} {3}'.format(old, decorate.BOLD, decorate.END, new))
        if not dry_run:
            status = mirror_file(old, new)
            if status == 0:
                print('done')
        print('')

def mirror_file(old_path, new_path):
    os.makedirs(os.path.dirname(new_path), exist_ok = True)
    shutil.copy2(old_path, new_path)
    return 0

if __name__ == "__main__":
    args = docopt(__doc__, version="pogfinder v0.0")

    mirror_files_from_tsv(
        args['--input'],
        args['--root'],
        args['--prefix'],
        args['--dry-run'],
        args['--col-name']
    )

