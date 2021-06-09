import requests, os
import pandas as pd
from itertools import zip_longest

ACCEPTABLE_IDENTIFERS = ['cid']


def bioassay_post(identifier_list, identifier='cid', output='csv'):
    """ uses PubChem's PUG-Rest service to obtain assay summaries for a target set of chemicals
    via a POST request.

    identifer_list: a list of chemical identifiers, the type of identifier should match with that outlined
    in the 'identifer' arguement with the default being PubChem Compound Identifier (CID).

    identifier: the chemical identifier (e.g., cid, smiles, etc.) of the chemicals in 'identifier_list'. Can be any chemical identifier
    as outlined in the PubChem PUG-Rest documentation.  Default=cid

    output: the output format (e.g., csv, json, etc.) of the

    """

    # convert list of identifers to str
    identifier_list = list(map(str, identifier_list))

    # make the base URL for the PubChem POST Request
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/assaysummary/{}'.format(identifier, output)

    # encoded_list = [requests.utils.quote(s) for s in identifier_list]

    # encoded_list = requests.utils.quote(','.join(identifier_list))
    headers = {'Content-Type': 'multipart/form-data'}
    data = {identifier: ','.join(identifier_list)}

    response = requests.post(url, data=data)

    return response


def grouper(iterable, n, fillvalue=None):
    """ support function for bioprofile, used to create n equaled-sized
    sets from an iterable.
    """

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def bioprofile(identifier_list,
               identifier='cid',
               outfile='bioprofile.csv',
               chunk=False,
               restrict_aids=None,
               ):
    """

    :param identifier_list: a list of identifiers to query
    :param identifier: right now only valid for cids
    :param outfile: name of the bioprofile to cache to
    :param chunk: batch size to query in
    :param restrict_aids: a list of aids to restrict to.  E.g.,
    any aid not in this list will not be in the final bioprofile
    :return:
    """

    if identifier not in ACCEPTABLE_IDENTIFERS:
        raise Exception('Sorry, currently only the following identifers are valid:'.format('\n\t'.join(ACCEPTABLE_IDENTIFERS)))

    num_compounds = len(identifier_list)

    # check to see whether
    # the list should be queried
    # in chunks of data or
    # processed as a whole
    if not chunk:
        chunk_size = num_compounds
    else:
        chunk_size = chunk

    counter = 0

    f = open(outfile, 'w')
    header_written = False

    for gp in grouper(identifier_list, chunk_size):
        batch = [cid for cid in gp if cid]
        response = bioassay_post(batch, identifier=identifier, output='csv')

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
            os.remove(outfile)
            print("Error: {}".format(response.status_code))
            return
    f.close()

    if restrict_aids:
        df = pd.read_csv(outfile)
        df = df[df.AID.isin(list(map(int, restrict_aids)))]
        df.to_csv(outfile, index=False)

def make_matrix(data_file,
                min_actives=0,
                merge_dups_by='max',
                outfile='bioprofile_matrix.csv',
                identifier_list=None
                ):
    """ will turn assay data file into wide (i.e., a matrix) format but performs all filtering in
        long format to save RAM"""

    df = pd.read_csv(data_file, usecols=['AID', 'CID', 'Bioactivity Outcome'])
    df['Bioactivity Outcome'] = df['Bioactivity Outcome'] \
        .replace('Inactive', -1) \
        .replace('Active', 1) \
        .replace('Probe', 1) \
        .replace('Inconclusive', 0) \
        .replace('Unspecified', 0) \
        .fillna(0)

    df['Activity Transformed'] = df.groupby(['CID', 'AID'])['Bioactivity Outcome'].transform(merge_dups_by)

    # for the num actives count, a CID is considered active
    # for a given AID if it has an active response for any of its
    # bioactivity outcomes for that AID.
    df['Bioactivity Outcome Max'] = df.groupby(['CID', 'AID'])['Bioactivity Outcome'].transform('max')

    # take only one response for a
    # CID/AID pair, ie., the transformed
    # bioactivity value
    df = df.drop_duplicates(['CID', 'AID', 'Activity Transformed'])
    # CID/AID/Bioactivity Outcome should be unique
    # just like CID/AID/Bioactivity Outcome Maz
    df_tmp = df.groupby('AID')['Bioactivity Outcome Max'].apply(lambda x: (x == 1).sum()).reset_index(name='Num_Active')
    df = df.merge(df_tmp)
    df = df[df.Num_Active >= min_actives]

    # turn into wide format
    matrix = df.pivot(index='CID', columns='AID', values='Activity Transformed').fillna(0)

    if identifier_list:
        # add compounnds that are
        # not in the matrix to have
        # all zeros, then reindex
        for cmp in identifier_list:
            if cmp not in matrix.index:
                matrix.loc[cmp] = [0]*matrix.shape[1]
        matrix = matrix.loc[identifier_list]
    matrix.to_csv(outfile)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Build bioprofile and heatmap from current PubChem data')

    parser.add_argument('-f', '--file', metavar='df', type=str,
                        help='Name of data file containing identifiers. (as txt)')

    parser.add_argument('-bs', '--batch_size', metavar='batch', type=int,
                        help='Name of data file containing casrns (as csv)')

    parser.add_argument('-ma', '--min_actives', metavar='ma', type=int,
                        help='Minimum number of active responses required for bioassays')

    parser.add_argument('-ra', '--restrict_aids', metavar='ma', type=str,
                        help='A comma-delimited list of aids to limit the bioprofile to.  E.g., '
                             'eliminate aids not in this list')

    args = parser.parse_args()
    datafile_name = args.file
    batch_size = args.batch_size if args.batch_size else 100
    min_actives = args.min_actives if args.min_actives else 0
    restrict_aids = args.restrict_aids.split(',') if args.restrict_aids else None


    # open the csv file and read
    # store identifiers as a list
    target_compounds = []

    csv_file = open(datafile_name, "r")

    for line in csv_file:
        cmp = line.strip()
        target_compounds.append(cmp)

    csv_file.close()

    # collects data in 'long' format
    # and writes to a csv file specified
    # in the outfile parameter
    bioprofile(target_compounds, chunk=batch_size, outfile='bioprofile_long.csv', restrict_aids=restrict_aids)

    # convert to a matrix that is
    # more suitable for modeling
    make_matrix('bioprofile_long.csv', min_actives=min_actives, outfile='bioprofile_matrix.csv', identifier_list=target_compounds)