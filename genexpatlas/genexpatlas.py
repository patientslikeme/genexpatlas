"""
Class for handling data retrieval from EMBL-EBI's Gene Expression Atlas.  Some methods translated from R package located
at http://www.bioconductor.org/packages/release/bioc/html/ExpressionAtlas.html
"""

import re
import requests
import json
import pandas as pd
import warnings
from operator import itemgetter
import urllib2
import xmltodict


def search_atlas_experiments(search=None, species='', summary=False):
    """
    Search against ArrayExpress API for Atlas datasets matching given terms
    :param search: Keywords describing search to be done. Defaults to None
    :param species: Desired species to search. Defaults to None
    :param summary: Return abbreviated list of experiment data, similar to R package. Defaults to True
    :type search: list
    :type species: str
    :type summary: bool
    :return: list dicts of results. If summary, sorted by species, experiment type, and accession
    """

    if search is None:
        search = []

    query = {'gxa': True}
    base_url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments"

    if search:
        if len(search) > 1:
            search = " OR ".join(search)
        query['keywords'] = search
    if species:
        query['species'] = species

    response = requests.get(base_url, params=query)
    response.raise_for_status()
    json_content = response.content

    # Parse JSON
    parsed = json.loads(json_content)

    # Get number of docs, return message if 0
    experiments_found = int(parsed['experiments']['total'])

    if experiments_found == 0:
        raise ValueError("No experiments found for search criteria {keywords} and species {species}".format(
            keywords=query['keywords'], species=query['species']))

    if summary is True:
        # Build list of dicts of the following for each experiment: accession, study name, organism, experimenttype
        result_set = [{'species': experiment['organism'], 'type': experiment['experimenttype'],
                       'accession': experiment['accession'], 'title': experiment['name']}
                      for experiment in parsed['experiments']['experiment']]

        # Sort by species, type, then accession
        return sorted(result_set, key=itemgetter('species', 'type', 'accession'))
    else:
        # Return full list of experiments
        return parsed['experiments']['experiment']


def get_atlas_experiments(experiments):
    """
    Download multiple experiments in one call and add to metadata
    :param experiments: a list of full experiments to pull full data from
    :type experiments: list
    :return: list of all successful results
    """

    valid_experiments = [x for x in experiments if __is_valid_experiment_accession(x['accession'])]

    if valid_experiments.__len__() == 0:
        raise ValueError("No valid accessions found")

    if valid_experiments.__len__() != experiments.__len__():
        warnings.warn("Invalid accessions found, will be removed", UserWarning)

    for experiment in valid_experiments:
        try:
            experiment['data'] = get_atlas_experiment(experiment)
        except Exception:
            warnings.warn("Experiment not in correct format, skipping")

    return valid_experiments


def get_atlas_experiment(experiment):
    """
    Downloads an individual experiment. We need a full experiment dict to construct a filename for differential studies
    :param experiment: an individual experiment dictionary
    :type experiment: dict
    :return: pandas DataFrame of experimental data
    """

    if not __is_valid_experiment_accession(experiment['accession']):
        raise ValueError("Invalid accession: {}".format(experiment['acession']))

    base_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
    file_name = experiment['accession']

    try:
        file_name = file_name + '/' + file_name + '_' + experiment['arraydesign'][0]['accession'] + '-analytics' + '.tsv'
        url = base_url + file_name
        loaded_data = pd.read_csv(url, sep='\t')
    except:
        try:
            file_name = file_name + '/' + file_name + '-analytics' + '.tsv'
            url = base_url + file_name
            loaded_data = pd.read_csv(url, sep='\t')
        except Exception as e:
            raise e

    # Get translations of contrast ids to comparison names
    configuration_url = base_url + experiment['accession'] + '/' + experiment['accession'] + '-configuration.xml'

    config_file = urllib2.urlopen(configuration_url)
    loaded_config = config_file.read()
    config_file.close()
    parsed_config = xmltodict.parse(loaded_config)
    compare_dict = {}
    for name, contrast in parsed_config['configuration']['analytics']['contrasts'].items():
        if type(contrast) is list:
            for ind_contrast in contrast:
                compare_dict[ind_contrast['@id']] = ind_contrast['name']
        else:
            compare_dict[contrast['@id']] = contrast['name']

    readable_data = __translate_data_headers(loaded_data, compare_dict)

    return readable_data


def get_atlas_experiment_summaries(accessions):
    """
    Download multiple R data summaries in one call, translated from R package
    :param accessions: a list of accessions to pull
    :type accessions: list
    :return: list of all successful results
    """

    result = []
    valid_accessions = [x for x in accessions if __is_valid_experiment_accession(x)]

    if len(valid_accessions) == 0:
        raise ValueError("No valid accessions found")

    if len(valid_accessions) != len(accessions):
        raise ValueError("No valid accessions found")

    for accession in valid_accessions:
        try:
            result.append(get_atlas_experiment_summary(accession))
        except ValueError as e:
            warnings.warn(e)

    return result


def get_atlas_experiment_summary(accession):
    """
    TODO: Download R data summary and load into Python-supported structure
    :param accession:
    :return:
    """
    pass


def __is_valid_experiment_accession(accession):
    """
    Helper method to determine if a given accession is valid
    :param accession: A single experiment accession
    :type accession: str
    :return: Boolean indicating if accession is valid
    """

    match = re.match(r"^E-\w{4}-\d+$", accession)
    if match is None:
        return False

    return True


def __translate_data_headers(experiment_data, translation_table):
    """
    Private method to translate dataframe headers to human-readable comparison names
    :param experiment_data: Loaded experiment data
    :type experiment_data: DataFrame
    :param translation_table: A table of translations of contrast IDs to contrast names
    :type translation_table: dict
    :return: pandas DataFrame with translated headers
    """

    columns = experiment_data.columns.values.tolist()
    trans_dict = {}
    for header in columns:
        if re.match(r"^g\d+_g\d+\..*", header):
            trans_dict[header] = translation_table[header.split('.')[0]] + '.' + header.split('.')[1]

    experiment_data = experiment_data.rename(columns=trans_dict)
    return experiment_data
