''' Systematically locates and outputs POG files

    November 2016
    Written by Eric Zhao, BC Genome Sciences Centre, ezhao@bcgsc.ca

    Usage: pogfinder.py (-b | -t | -g | -p | -a | -s | -n | -v | -V | -I | -i | -d | -D | -A) [ -c CASES -l LOGPATH ]

    Options:
        -b --blood-bam              Germline BAM files
        -t --tumour-bam             Tumour BAM files
        -g --germline               Germline excel files
        -p --germline-pathogenic    Pathogenic germline excel files
        -a --germline-annotated     Annotated germline excel files
        -s --loh-seg                LOH seg files
        -n --cnv                    LOH seg files
        -v --somatic-snv            Somatic SNVs - Unannotated
        -V --somatic-snv-ann        Somatic SNVs - Annotated
        -I --indel                  Somatic INDELS
        -i --signature              Somatic Signatures
        -d --drug-target            Drug target analysis summaries
        -D --delly                  DELLY raw structural variant calls
        -A --trans-abyss            Remapped SV contigs from trans-abyss
        -c CASES                    Path to file with POG case IDs, one per line
        -l LOGPATH                  Path to output log file to
'''

from docopt import docopt
import glob
import re

def read_flatfile(flatfile_path):
    flatfile = open(flatfile_path).readlines()
    header = flatfile[0].strip().split('\t')
    flatfile_data = flatfile[1:len(flatfile)]

    table = {header : [line.split('\t')[i] for line in flatfile_data]
             for i, header in enumerate(header)}
    return(table)

def search_flatfiles(column_names):
    pog_re = re.compile('(POG\d+)')
    flatfile_paths = glob.glob('/projects/POG/POG_data/POG*/flatfile/POG???.tab')
    sample_mapping = {pog_re.search(path).group(1) : path for path in flatfile_paths}

    flatfile_combined = {col: [] for col in column_names}
    flatfile_combined['sample'] = []

    for sample in sorted(sample_mapping.keys()):
        flatfile_path = sample_mapping[sample]
        flatfile = read_flatfile(flatfile_path)

        skip_entry = False
        for col in column_names:
            if not col in flatfile:
                print('Warning: Column named "{0}" not found in flatfile {1}'.format(
                    col, flatfile_path)
                )
                skip_entry = True
                break
        if skip_entry:
            continue

        for col in column_names:
            flatfile_combined[col] += flatfile[col]
        flatfile_combined['sample'] += ([sample] * len(flatfile[column_names[0]]))

    return flatfile_combined

def get_bam_paths(tissue):
    flatfiles = search_flatfiles(['sample_prefix', 'merged_bam'])
    print(flatfiles)

    nt_re = re.compile(r'(\d+nt)')

    bam_globs = [glob.glob(path + '/*.bam') for path in flatfiles['merged_bam']]
    bam_paths = [bamglob[0] if bamglob else '' for bamglob in bam_globs]
    nt_strings = [nt_re.search(path).group(1) if nt_re.search(path) is not None else ''
            for path in bam_paths]

    output = [
        [bam_paths[i], "_".join([flatfiles['sample'][i],
                                 flatfiles['sample_prefix'][i],
                                 nt_strings[i]])
        ]
        for i in range(len(flatfiles['sample']))
        if tissue in flatfiles['sample_prefix'][i] and bam_paths[i].strip()
    ]

    output_strings = ['\t'.join(row) for row in output]

    return (output_strings, '')

def get_germline_paths(pathogenic=False, annotated=False):
    ''' Retrieves paths for all germline files
    '''
    germline_annotated_glob = "/projects/tumour_char/pog/germline/small_mutations/snv/POG???/*_hg19/v0.4.0/POG*.germline.annotated.xls"
    germline_pathogenic_glob = "/projects/tumour_char/pog/germline/small_mutations/snv/POG???/*_hg19/v0.4.0/POG*.germline.pathogenic.xls"
    germline_glob = "/projects/tumour_char/pog/germline/small_mutations/snv/POG???/*_hg19/v0.4.0/POG*.germline.xls"

    if pathogenic:
        globstring = germline_pathogenic_glob
    elif annotated:
        globstring = germline_annotated_glob
    else:
        globstring = germline_glob

    paths = glob.glob(globstring)
    log = ''
    return (paths, log)

def get_drug_target_paths():
    drug_target_glob = '/projects/POG/POG_data/POG???/drug_target_analysis/transient_latest_run/biop*blood*/summary/gene_drug_summary.tsv'
    paths = glob.glob(drug_target_glob)
    log = ''
    return(paths, log)

def get_loh_paths():
    loh_glob = '/projects/POG/POG_data/POG*/wgs/*/reviewed/loh/results/apolloh_out_segs.txt'
    paths = glob.glob(loh_glob)
    log = ''
    return(paths, log)

def get_cnv_paths():
    cnv_glob = '/projects/POG/POG_data/POG???/wgs/biop?_t_*_blood?_n_*/reviewed/cnv/w200/w200_segs.txt'
    paths = glob.glob(cnv_glob)
    log = ''
    return(paths, log)

def get_snv_paths(annotated=False):
    if annotated:
        snv_glob_1 = '/projects/POG/POG_data/POG*/wgs/*/*/strelka/*/bwa/results/passed.somatic.snvs.eff.snvs.vcf'
        snv_glob_2 = '/projects/POG/POG_data/POG*/wgs/*/*/strelka/*/bwa/results/passed.somatic.snvs.eff.vcf'
        paths = glob.glob(snv_glob_1) + glob.glob(snv_glob_2)
    else:
        snv_glob = '/projects/POG/POG_data/POG*/wgs/*/*/strelka/*/bwa/results/passed.somatic.snvs.vcf'
        paths = glob.glob(snv_glob)
    log = ''
    return(paths, log)

def get_indel_paths():
    indel_glob_1 = '/projects/POG/POG_data/POG*/wgs/*/*/strelka/*/bwa/results/passed.somatic.indels.eff.snvs.vcf'
    indel_glob_2 = '/projects/POG/POG_data/POG*/wgs/*/*/strelka/*/bwa/results/passed.somatic.indels.eff.vcf'
    paths = glob.glob(indel_glob_1) + glob.glob(indel_glob_2)
    log = ''
    return(paths, log)

def get_signature_paths():
    paths = glob.glob('/projects/tumour_char/pog/somatic/signature/POG*/biopsy_normal/v0.3.0/*nnls.txt')
    return(paths, '')

def get_delly_paths():
    delly_glob_1 = '/projects/POG/POG_data/POG???/sv/delly/delly-?.?.?/*/Somatic_Germline_Quality_tagged.vcf'
    delly_glob_2 = '/projects/POG/POG_data/POG???/sv/delly/delly-?.?.?/*/Somatic_Germline_Quality_tagged.vcf_bak'
    paths = glob.glob(delly_glob_1) + glob.glob(delly_glob_2)
    log = ''
    return(paths, log)

def get_contig_remap_paths():
    paths = glob.glob('/projects/POG/POG_data/POG*/wgs/*/trans-ABySS/fusions/LSR.tsv')
    return(paths, '')

def limit_to_case_list(paths, case_list):
    case_set = set(case_list)

    pog_re = re.compile('(POG\d+)')
    pog_ids = [pog_re.search(path).group(1) for path in paths]

    new_paths = [paths[i] for i, pogid in enumerate(pog_ids) if pogid in case_set]

    return new_paths


def write_to_log(data, path):
    handle = open(path, 'w')
    handle.write(data)
    handle.close()
    return 0

if __name__ == "__main__":
    args = docopt(__doc__, version="pogfinder v0.0")

    if args['--germline']:
        (paths, log) = get_germline_paths()
    elif args['--germline-pathogenic']:
        (paths, log) = get_germline_paths(pathogenic=True)
    elif args['--germline-annotated']:
        (paths, log) = get_germline_paths(annotated=True)
    elif args['--blood-bam']:
        (paths, log) = get_bam_paths(tissue='blood')
    elif args['--tumour-bam']:
        (paths, log) = get_bam_paths(tissue='biop')
    elif args['--loh-seg']:
        (paths, log) = get_loh_paths()
    elif args['--drug-target']:
        (paths, log) = get_drug_target_paths()
    elif args['--cnv']:
        (paths, log) = get_cnv_paths()
    elif args['--somatic-snv']:
        (paths, log) = get_snv_paths(annotated=False)
    elif args['--somatic-snv-ann']:
        (paths, log) = get_snv_paths(annotated=True)
    elif args['--indel']:
        (paths, log) = get_indel_paths()
    elif args['--signature']:
        (paths, log) = get_signature_paths()
    elif args['--delly']:
        (paths, log) = get_delly_paths()
    elif args['--trans-abyss']:
        (paths, log) = get_contig_remap_paths()
    else:
        print('No acceptable file type received.')

    if args['-c'] is not None:
        paths = limit_to_case_list(paths, [id.strip() for id in open(args['-c'])])

    print('\n'.join(paths))

    if args['-l'] is not None:
        write_to_log(log, args['-l'])
