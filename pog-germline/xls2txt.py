''' xls2txt.py - Conversion of excel variant calls to TXT files

    November, 2016
    Written by Eric Zhao, BC Genome Sciences Centre, ezhao@bcgsc.ca
    --

    Usage: xls2txt.py -x EXCELPATH [ -l LOGPATH ]
           xls2txt.py -g GLOB [ -l LOGPATH ]

    Options:
        -x EXCELPATH    Input excel file path - outputs to STDOUT
        -g GLOB         Path string with wildcards to excel files
        -l LOGPATH      Path to output log file (optional)
'''

from docopt import docopt
import glob
import xlrd

def get_header_row(sheet):
    ''' Finds the header row and returns the tuple (header, idx)
        where header is a list and idx is the row number of the header
    '''
    rows = [[sheet.cell(i, j).value for j in range(sheet.ncols)]
            for i in range(sheet.nrows)]

    for i, row in enumerate(rows):
        if "Chr" in row and "Pos" in row:
            return (row, i)

    HeaderError('Could not find Header Row. Is this a valid POG Germline Excel File?')

def convert_excel(excel_path):
    wb = xlrd.open_workbook(excel_path)
    sheet = wb.sheets()[0]
    (header, header_index) = get_header_row(sheet)

    #expression_indices = []
    #for i, value in enumerate(header):
    #    if 'TCGA' in value or 'FC' in value or 'Expression' in value:
    #        expression_indices.append(i)

    col_names = ['Sample Name', 'Chr', 'Pos', 'Ref', 'Alt', 'Variant', 'Gene',
            'dbSNP', 'ClinVar', 'Zygosity in germline', 'Zygosity in tumour', 'Type',
            'Flagged', 'GMAF', 'Histology', 'Mutational landscape', 'Transcript',
            'Complete annotation', 'Patient history', 'Notes'
            ]
    col_indices = [header.index(c) if c in header else None for c in col_names] # + expression_indices

    output_table = [[sheet.cell(row_idx, col_idx).value if col_idx is not None else 'NA' for col_idx in col_indices]
            for row_idx in range(header_index, sheet.nrows)]

    return output_table

def table_to_string(table):
    try:
        output_string = '\n'.join(['\t'.join([''.join([i if ord(i) < 128 else '' for i in unicode(r)]) for r in row]) for row in table]) + '\n'
        return output_string
    except:
        print('ERROR IN THIS LINE:')
        print(row)

def write_to_file(content, output_path):
    handle = open(output_path, 'w')
    handle.write(content)
    handle.close()

def run_multiple(paths, output_dir):
    pass

class HeaderError(Exception):
    pass

if __name__ == "__main__":
    args = docopt(__doc__, version="xls2txt v0.0")

    if args['-x']:
        content = convert_excel(args['-x'])
        print(table_to_string(content))

    elif args['-g']:
        paths = glob.glob(args['-g'])
        out_paths = [path.strip() + '.txt' for path in paths]
        for i, path in enumerate(paths):
            out_path = out_paths[i]
            write_to_file(table_to_string(convert_excel(path)), out_path)

    if args['-l']:
        write_to_file('', args['-l'])
