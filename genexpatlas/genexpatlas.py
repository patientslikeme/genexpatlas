"""
Class for handling data retrieval from EMBL-EBI's Gene Expression Atlas.  Some methods translated from R package located
at http://www.bioconductor.org/packages/release/bioc/html/ExpressionAtlas.html
"""

import re
import requests
import json
import pandas as pd
from operator import itemgetter
import urllib2
import xmltodict
import os
import errno

base_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"


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
    :return: generator
    """

    valid_experiments = (x for x in experiments if __is_valid_experiment_accession(x['accession']))
    for experiment in valid_experiments:
        try:
            exp_data = get_atlas_experiment(experiment)
            experiment['data'] = exp_data[0]
            experiment['contrasts'] = exp_data[1]
            yield experiment
        # We want to ignore some expected errors, otherwise iteration would stop
        except ValueError:
            pass
        except IOError:
            pass


def get_atlas_experiment(experiment):
    """
    Downloads an individual experiment. We need a full experiment dict to construct a filename for differential studies
    :param experiment: an individual experiment dictionary
    :type experiment: dict
    :return: pandas DataFrame of experimental data
    """

    if not __is_valid_experiment_accession(experiment['accession']):
        raise ValueError("Invalid accession: {}".format(experiment['acession']))

    file_name = experiment['accession']

    try:
        file_name = file_name + '/' + file_name + '_' + experiment['arraydesign'][0]['accession'] + '-analytics.tsv'
        url = base_url + file_name
        loaded_data = pd.read_csv(url, sep='\t')
    except:
        try:
            file_name = file_name + '/' + file_name + '-analytics.tsv'
            url = base_url + file_name
            loaded_data = pd.read_csv(url, sep='\t')
        except IOError as e:
            raise e

    # Get translations of contrast ids to comparison names
    configuration_url = base_url + experiment['accession'] + '/' + experiment['accession'] + '-configuration.xml'

    config_file = urllib2.urlopen(configuration_url)
    loaded_config = config_file.read()
    config_file.close()
    compare_dict = __get_comparison_translations(xmltodict.parse(loaded_config))
    readable_data = __translate_data_headers(loaded_data, compare_dict)

    return readable_data, compare_dict.values()


def get_atlas_experiment_summaries(accessions, directory):
    """
    Download multiple R data summaries in one call, translated from R package
    :param accessions: a list of accessions to pull
    :type accessions: list
    :param directory: directory to download files to
    :type directory: str
    :return: list of all successful results
    """

    result = []
    valid_accessions = [x for x in accessions if __is_valid_experiment_accession(x)]

    if len(valid_accessions) == 0:
        raise ValueError("No valid accessions found")

    if len(valid_accessions) != len(accessions):
        raise ValueError("No valid accessions found")

    # Set up download directory
    directory = __dir_setup(directory)

    for accession in valid_accessions:
        try:
            result.append(get_atlas_experiment_summary(accession, directory))
        except (ValueError, IOError, urllib2.URLError) as e:
            raise e

    return result


def get_atlas_experiment_summary(accession, directory):
    """
    Currently, downloads RData summary to supplied directory.
    TODO: load into Python-supported structure, likely with rpy2
    :param accession: a string of an experiment's accession
    :type accessions: list
    :param directory: directory to download files to
    :type directory: str
    :return: successful accession
    """

    filename = accession + "-atlasExperimentSummary.Rdata"
    full_url = base_url + "/" + accession + "/" + filename

    try:
        f = urllib2.urlopen(full_url)
    except urllib2.URLError as e:
        raise e

    try:
        with open(directory + filename, "wb") as d:
            d.write(f.read())
    except IOError as e:
        raise e

    return accession


def __get_comparison_translations(parsed_config):
    """
    Builds a translation dictionary of g1, etc. abbreviations to human-readable comparison names
    :param parsed_config: a parsed XML configuration file
    :return: translation dictionary
    """
    compare_dict = {}

    if type(parsed_config['configuration']['analytics']) is list:
        for item in parsed_config['configuration']['analytics']:
            for name, contrast in item['contrasts'].items():
                if type(contrast) is list:
                    for ind_contrast in contrast:
                        compare_dict[ind_contrast['@id']] = ind_contrast['name']
                else:
                    compare_dict[contrast['@id']] = contrast['name']
    else:
        for name, contrast in parsed_config['configuration']['analytics']['contrasts'].items():
            if type(contrast) is list:
                for ind_contrast in contrast:
                    compare_dict[ind_contrast['@id']] = ind_contrast['name']
            else:
                compare_dict[contrast['@id']] = contrast['name']

    return compare_dict


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


def __dir_setup(directory):
    """
    Sets up data directory
    :return: directory name
    """

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    return directory


def run_execution_tests():
    print("Running simple execution test that determines that code paths run without error")

    print("search_atlas_experiments with summary data")
    exp = search_atlas_experiments(species='homo sapiens', summary=True)

    print("get_atlas_experiment with first experiment from search results")
    data = get_atlas_experiment(exp[0])

    print("get_atlas_experiment_summaries with first 3 experiments from search results")
    try:
        summary = get_atlas_experiment_summaries([x['accession'] for x in exp[0:2]], directory=os.getcwd())
    except Exception as e:
        print(e.msg)

    print("search_atlas_experiments without summary data")
    exp = search_atlas_experiments(species='homo sapiens', summary=False)

    print("get_atlas_experiments with all results (produces generator)")
    data = get_atlas_experiments(exp)

    print("Methods executed successfully")


if __name__ == '__main__':
    run_execution_tests()
