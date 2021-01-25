import argparse
import os, pandas as pd
import requests
from typing import List
import urllib
from rdkit import Chem
import numpy as np
from itertools import zip_longest

from get_and_clean import clean_bioprofile, get_activity_vectors


def bioassay_post(identifier: str, identifier_list: List, output='csv'):
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/assaysummary/{}'.format(identifier, output)

    # encoded_list = [requests.utils.quote(s) for s in identifier_list]

    #encoded_list = requests.utils.quote(','.join(identifier_list))
    headers = {'Content-Type': 'multipart/form-data'}
    data = {identifier: ','.join(identifier_list)}

    response = requests.post(url, data=data)

    return response


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def chunk_responses(identifier: str, identifier_list: List, chunk_size=100, outfile='assay_data'):

    num_compounds = len(identifier_list)
    counter = 0

    f = open(outfile, 'w')
    header_written = False

    for gp in grouper(identifier_list, chunk_size):
        batch = [cid for cid in gp if cid]
        response = bioassay_post(identifier, batch, output='csv')

        if response.status_code == 200:
            text = response.text
            if not header_written:
                f.write(text)
                header_written = True
            else:
                header = text.split('\n')[0]
                text = text.replace(header, '')
                f.write(text)
            counter = counter + len(batch)
            print("Retrieved data for {} out of {} compounds.".format(counter, num_compounds))
        else:
            f.close()
            os.remove('data/assay_data.csv')
            print("Error: {}".format(response.status_code))
            return
    f.close()


def make_matrix(data_file, min_actives=0):
    """ will turn assay data file into wide (i.e., a matrix) format but perform all filtering in
        long format to save RAM"""

    df = pd.read_csv(data_file, usecols=['AID', 'CID', 'Bioactivity Outcome'])
    df = df.drop_duplicates(subset=['AID', 'CID'])


    df_tmp = df.groupby('AID')['Bioactivity Outcome'].apply(lambda x: (x == 'Active').sum()).reset_index(name='Num_Active')
    df = df.merge(df_tmp)
    print(df)
    df = df[df.Num_Active >= min_actives]

    matrix = df.pivot(index='CID', columns='AID', values='Bioactivity Outcome')
    matrix = matrix.replace('Inactive', -1).replace('Active', 1).replace('Probe', 1).replace('Inconclusive', 0).replace('Unspecified', 0).fillna(0)

    # eliminate
    new_file = data_file.split('.')[0] + '_matrix.csv'
    matrix.to_csv(new_file)




if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Build bioprofile and heatmap from current PubChem data')

    parser.add_argument('-f', '--cidsfile', metavar='df', type=str,
                        help='Name of data file containing cids (as txt)')

    parser.add_argument('-bs', '--batch_size', metavar='batch', type=int,
                        help='Name of data file containing casrns (as csv)')

    args = parser.parse_args()
    datafile_name = args.cidsfile
    batch_size = args.batch_size

    df = pd.read_csv(datafile_name, header=None, sep='\t')

    cids = [str(int(cid)) for cid in df[1].tolist() if not np.isnan(cid)]

    assay_file = datafile_name.split('.')[0] + 'assay_data.csv'
    chunk_responses('cid', cids, chunk_size=batch_size, outfile=assay_file)
    make_matrix(assay_file, min_actives=5)
