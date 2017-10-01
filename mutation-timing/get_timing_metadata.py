''' get_timing_metadata.py - Searches through directory tree to construct mutation timing config file

    Usage: get_timing_metadata.py

'''

import glob
import re

def paths_to_table(paths):
    tissue_re = re.compile('([bioparchld]+\d+)')
    pog_re = re.compile('(POG\d+)')
    ids = [pog_re.search(path).group(1) for path in paths]
    tissues = [tissue_re.findall(path) for path in paths]
    if len(tissues[0]) == 1:
        table = [{'id': '_'.join([ids[i], tissues[i][0]]), 'path': paths[i]}
            for i in range(len(ids))]
    elif len(tissues[0]) == 2:
        table = [{'id': '_'.join([ids[i], tissues[i][0], tissues[i][1]]), 'path': paths[i]}
            for i in range(len(ids))]
    elif len(tissues[0]) == 0:
        table = [{'id': ids[i], 'path': paths[i]}
            for i in range(len(ids))]
    return(table)

def combine_dicts(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d

def merge_tables(t1, t2, by='id', name1='t1', name2='t2'):
    t1_pre = [combine_dicts({name1 + '_' + key: row[key] for key in row if key != by}, {by: row[by]})
            for row in t1]
    t2_pre = [combine_dicts({name2 + '_' + key: row[key] for key in row if key != by}, {by: row[by]})
            for row in t2]

    merged = []
    for d1 in t1_pre:
        for d2 in t2_pre:
            if d1[by] == d2[by]:
                merged.append(combine_dicts(d1, d2))

    return(merged)


def get_tumor_content_dict(flatfile_path):
    pog_re = re.compile('(POG\d+)')
    handle = open(flatfile_path)
    header = next(handle).strip().split('\t')
    tc_index = header.index('biofx_tc')
    prefix_index = header.index('sample_prefix')

    tc_dict = {}

    for line in handle:
        row = line.strip().split('\t')
        id = pog_re.search(line).group(1) + '_' + row[prefix_index].strip()
        tc = row[tc_index].strip()
        tc_dict[id] = tc

    return(tc_dict)

def has_pathogenic_calls(germline_path):
    content = open(germline_path).read().strip().split('\n')
    return len(content) > 1

def get_timing_metadata():
    vcf_paths = glob.glob('data/POG*/somatic_snv/*snvs.vcf')
    loh_paths = glob.glob('data/POG*/*.apolloh_out_segs.txt')
    germline_paths = glob.glob('data/POG*/POG*.germline.pathogenic.xls.txt')
    germline_paths = [g for g in germline_paths if has_pathogenic_calls(g)]

    vcf_table = paths_to_table(vcf_paths)
    loh_table = paths_to_table(loh_paths)
    germline_table = paths_to_table(germline_paths)
    included_cases = set([d['id'] for d in germline_table])

    merged = merge_tables(vcf_table, loh_table, by='id', name1='vcf', name2='loh')
    included = [merged[i] for i in range(len(merged)) if merged[i]['id'].split('_')[0] in included_cases
            if not 'POG000' in merged[i]['id']]
    input_paths = ['data/{0}/mutation_timing/{1}.Rdata'.format(d['id'].split('_')[0], d['id']) for d in included]
    flatfile_paths = [glob.glob('/projects/POG/POG_data/{0}/flatfile/{0}.tab'.format(d['id'].split('_')[0]))[0] for d in included]
    tc_dicts = [get_tumor_content_dict(f) for f in flatfile_paths]
    tumor_content = [tc_dicts[i]['_'.join(included[i]['id'].split('_')[0:2])] for i in range(len(included))]
    contamination = [str(1.0 - (float(tc)/100.0)) for tc in tumor_content]

    output = ['\t'.join([included[i]['id'],
                               included[i]['loh_path'],
                               included[i]['vcf_path'],
                               input_paths[i],
                               contamination[i]])
                        for i in range(len(included))]

    print('\t'.join(['id', 'loh', 'vcf', 'timing_input', 'contamination']))
    print('\n'.join(output))


if __name__ == "__main__":
    get_timing_metadata()
