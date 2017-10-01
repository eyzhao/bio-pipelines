''' query_bams.py - Queries positions in BAM files indicated line by line

    November 2016
    Written by Eric Zhao - BC Genome Sciences Centre - ezhao@bcgsc.ca / eyzhao.com

    Usage: query_bams.py -l POSITIONFILE [ -1 -a ]

    Options:
        -l POSITIONFILE         Path to a file where first 5 cols are chr, pos, ref, alt, bampath
        -1 --header             File has header - ignore first row
        -a --append             Append position file columns 3 and onwards
'''

from docopt import docopt
import re
import pysam

class TestError(Exception):
    pass

def parse_position(position):
    return {'chr': position.split(':')[0], 'pos': int(position.split(':')[1])}

def query_position(bampath, position, ref, alt):
    locus = parse_position(position)
    samfile = pysam.AlignmentFile(bampath, 'rb')

    pileup_obj = samfile.pileup(locus['chr'], locus['pos']-1, locus['pos'])
    pileup = [pileupcolumn.pileups for pileupcolumn in pileup_obj
            if pileupcolumn.pos == locus['pos'] - 1]
    if len(pileup) == 0:
        return None

    pileup = pileup[0]
    bases = [pileupread.alignment.query_sequence[pileupread.query_position]
            for pileupread in pileup if pileupread.query_position is not None]
    alleles = {allele: len([b for b in bases if b.upper() == allele])
            for allele in ['A', 'C', 'G', 'T']}
    indel = [pileupread.indel for pileupread in pileup]
    indelcount = len([i for i in indel if i != 0])

    if alt == 'NA':
        alt = max(alleles, key=lambda key: alleles[key] if not key == ref else -1)

    if len(ref) > 1 or len(alt) > 1:
        var_allele_count = indelcount
    else:
        var_allele_count = alleles[alt]

    row = [locus['chr'], locus['pos'], len(bases), var_allele_count, ref, alt,
            ] + [alleles[base] for base in ['A', 'C', 'G', 'T']
            ] + [indelcount]

    return([str(item) for item in row])

if __name__ == '__main__':
    args = docopt(__doc__)
    num_re = re.compile(r'\d+')
    handle = open(args['-l'])

    if args['--header']:
        input_header = next(handle).strip().split('\t')

    header = ['chr', 'pos', 'depth', 'var_depth', 'Ref', 'Alt', 'A', 'C', 'G', 'T', 'indel']

    if args['--append']:
        header = header + input_header[4:len(input_header)]

    print('\t'.join(header))

    for line in handle:
        if line.strip():
            row = line.strip().split('\t')
            if not num_re.match(row[1]):
                continue
            chr = row[0]
            pos = int(row[1])
            ref = row[2].upper()
            alt = row[3].upper()
            bampath = row[4]

            out = query_position(bampath, '{0}:{1}'.format(chr, pos), ref, alt)

            if out is not None:
                if args['--append']:
                    additional_column_values = row[4:len(row)]
                    out = out + additional_column_values

                print('\t'.join(out))


