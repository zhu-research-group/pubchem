import logging
import pandas as pd
import requests

from io import StringIO
from pandas.errors import ParserError
from time import sleep
from urllib.error import HTTPError, URLError


def get_pubchem_activities(url, identifier):
    s = requests.get(url).content
    c = pd.read_csv(StringIO(s.decode('utf-8'))).set_index('AID')
    activities = c.loc[:, 'Bioactivity Outcome']
    activities.rename(f'{identifier}', inplace=True)
    duplicates = activities[activities.index.duplicated()]
    cleaned_activities = activities[~activities.index.duplicated()]

    if duplicates.shape[0] > 0:
        duplicate_aids = duplicates[~duplicates.duplicated(keep='first')].index
        for aid in duplicate_aids:
            aid_activities = activities.loc[aid]
            if len(list(aid_activities.unique())) > 1:
                cleaned_activities.loc[aid] = 0

    return cleaned_activities


def get_activity_vectors(identifier_type, identifier_list):
    activity_vectors = []
    skipped = []
    log = logging.getLogger(__name__)

    for name in identifier_list:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{identifier_type}/{name}/assaysummary/CSV'

        try:
            activity_vectors.append(get_pubchem_activities(url, name))
        except KeyError:
            if identifier_type is 'name':
                try:
                    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{identifier_type}/CAS-{name}/' \
                          f'assaysummary/CSV'
                    activity_vectors.append(get_pubchem_activities(url, name))
                except (KeyError, ParserError):
                    skipped.append(name)
            else:
                skipped.append(name)
        except ParserError:
            skipped.append(name)
        except (HTTPError, URLError, TimeoutError) as err:
            log.error(err)
            skipped.append(name)
        sleep(0.5)

    return activity_vectors, skipped


def clean_bioprofile(bioprofile, min_actives=0):
    bioprofile.replace('Active', 1, True)
    bioprofile.replace('Inactive', -1, True)

    for col in bioprofile.columns:
        bioprofile[col] = pd.to_numeric(bioprofile[col], errors='coerce').fillna(0)

    bioprofile = bioprofile.transpose()

    return bioprofile.loc[:, ((bioprofile > 0).sum(axis=0) >= min_actives)]