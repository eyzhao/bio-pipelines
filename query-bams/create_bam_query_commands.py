''' create_bam_query_commands.py - Runs bam_query on the appropriate POG files

    Usage: create_bam_query_commands.py [ -l LOGFILE ]
'''

from docopt import docopt
import os
import glob
import re

def main(log):
    allpaths = glob.glob('data/*/*.germline.xls.txt.clean') \
            + glob.glob('data/*/*.germline.pathogenic.xls.txt.clean') \
            + glob.glob('data/*/*.bam')
    patient_dirs = glob.glob('data/*')

    pathsets = {pd: [path for path in allpaths if path.startswith(pd)]
            for pd in patient_dirs}

    pog_id_re = re.compile('(POG\d+)')

    for (pd, ps) in pathsets.items():
        germline = [p for p in ps if p.endswith('germline.xls.txt.clean')]
        pathogenic = [p for p in ps if p.endswith('pathogenic.xls.txt.clean')]
        bam = [p for p in ps if p.endswith('.bam')]
        patient_id_search = pog_id_re.search(pd)

        if patient_id_search is not None:
            patient_id = patient_id_search.group(1)
            if germline and bam:
                process_path_set(patient_id, pd, germline, None, bam)
            if pathogenic and bam:
                process_path_set(patient_id, pd, None, pathogenic, bam)

    logfile = open(log, 'w')
    logfile.write('Successfully created BAM query commands')
    logfile.close()

def process_path_set(patient_id, patient_dir, germline, pathogenic, bam):
    if germline is None and pathogenic is None:
        PathSetError('Only one of germline or pathogenic can be None')
    elif germline is None:
        depth_dir = patient_dir + '/pathogenic_germline_depth'
        path = pathogenic[0]
    elif pathogenic is None:
        depth_dir = patient_dir + '/germline_depth'
        path = germline[0]
    else:
        PathSetError('One of germline / pathogenic must be None')

    print('mkdir -p {0}'.format(depth_dir))

    for bam_path in bam:
        bam_fn = bam_path.split('/')[-1]
        print('python scripts/query-bams/query_bams.py -1 -a -b {0} -l {1} > {2}'.format(
            bam_path,
            path,
            depth_dir + '/{0}_germline_{1}.txt'.format(patient_id, bam_fn)
        ))


if __name__ == '__main__':
    args = docopt(__doc__, version='v0.0')
    main(log=args['LOGFILE'])
