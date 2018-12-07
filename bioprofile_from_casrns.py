import argparse
import os

from get_and_clean import clean_bioprofile, get_activity_vectors

parser = argparse.ArgumentParser(description='Build bioprofile and heatmap from current PubChem data')
parser.add_argument('c', '--cluster', metavar='c', type=str,
                    help='Are you running on the cluster? (Y/N)')
parser.add_argument('-df', '--datafile', metavar='df', type=str,
                    help='Name of data file containing casrns (as csv)')
parser.add_argument('-ev', '--env_var', metavar='ev', type=str,
                    help='Environmental variable of project directory')
parser.add_argument('ma', '--min_actives', metavar='ma', type=int,
                    help='Minimum number of active responses required for bioassays')

args = parser.parse_args()
datafile_name = args.datafile
directory = os.getenv(args.env_var)

if args.cluster in ['Y', 'y']:
    import matplotlib
    matplotlib.use('Agg')

import pandas as pd
import seaborn as sns

casrns = pd.read_csv(os.path.join(directory, datafile_name), header=None)
activity_vectors, skipped = get_activity_vectors('name', casrns[0])
bioprofile = clean_bioprofile(pd.concat(activity_vectors, axis=1), args.min_actives)
bioprofile.to_csv(os.path.join(directory, f'{datafile_name}_bioprofile.csv'))
heatmap = sns.clustermap(bioprofile, col_cluster=True, cmap='coolwarm', xticklabels=False, yticklabels=False)
heatmap.savefig(os.path.join(directory, f'{datafile_name}_bioprofile_heatmap.png'))

if len(skipped) > 0:
    pd.Series(skipped).to_csv(os.path.join(directory, f'{datafile_name}_skipped_cpds.csv'))
