''' Formats clinical OaSIS treatment regimens output into relational style text files.
    Works on tab-delimited converted version of OaSIS spreadsheet.
    Make sure that the TSV file has unix line endings, or else this script will not work.

    Usage: format.py -i INPUT -p PATIENTOUT -t TXOUT

    Options:
        -i INPUT        Input file - OaSIS dump
        -p PATIENTOUT   Output patient information file path
        -t TXOUT        Output treatment information file path
'''

from docopt import docopt
import pandas as pd


def date_parse_function(array):
    return pd.to_datetime(array, yearfirst = True)

def parse_table(input_path):
    table = pd.read_table(input_path, sep='\t', header=0, index_col=None, date_parser = date_parse_function,
            parse_dates = [4, 9, 10, 15, 21, 22, 27, 33, 34, 39, 45, 46, 51, 57, 58, 63,
                69, 70, 75, 81, 82, 87, 93, 94, 99, 105, 106, 111, 117, 118, 123,
                129, 130, 135, 137, 139, 140, 143, 146, 158, 170, 171, 176])

    ### TREATMENTS TABLE

    cc_indices = [
            [header_item.endswith('_cc' + str(cc_idx)) for header_item in list(table)]
            for cc_idx in range(1,13)
            ]

    cc_tables = [table.iloc[:, indices] for indices in cc_indices]
    pd.options.mode.chained_assignment = None
    for df in cc_tables:
        cc_number = int(list(df)[-1].split('_')[-1].replace('cc', ''))
        df['treatment_number_cc'] = cc_number
        df.columns = ['_'.join(c.split('_')[0:-1]) for c in df.columns]

    treatments_table = pd.concat(cc_tables)
    pog_id = table['gsc_pog_id'][treatments_table.index]
    treatments_table = pd.concat([pog_id, treatments_table], axis=1)

    treatments_table = treatments_table.loc[[
        not b for b in list(pd.isnull(treatments_table['course_begin_on']))
        ]].reset_index(drop=True)

    ### PATIENTS TABLE

    patient_indices = [not '_cc' in header_item for header_item in list(table)]
    print(patient_indices)
    patient_table = table.iloc[:, patient_indices]

    return {'treatment':treatments_table, 'patient':patient_table}

def write_table(table, path):
    table.to_csv(path, sep='\t', index=False)

if __name__ == '__main__':
    args = docopt(__doc__)
    output = parse_table(args['-i'])
    write_table(output['treatment'], args['-t'])
    write_table(output['patient'], args['-p'])
