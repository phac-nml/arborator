import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import json
import os
from datetime import datetime
import pandas as pd
import numpy as np
from statistics import mean, median
import shutil
from arborator.version import __version__
from arborator.classes.aggregator import summarizer
from profile_dists.utils import convert_profiles, calc_distances_hamming, process_profile
from arborator.classes.read_data import read_data
from arborator.classes.report import report
from arborator.classes.split_profiles import split_profiles
from genomic_address_service.classes.multi_level_clustering import multi_level_clustering
from genomic_address_service.utils import format_threshold_map
from genomic_address_service.mcluster import write_clusters
from genomic_address_service.constants import CLUSTER_METHODS
import fastparquet as pq
from multiprocessing import Pool, cpu_count

# ARGUMENTS
PROFILE_KEY = "profile"
PROFILE_LONG = "--" + PROFILE_KEY
PROFILE_SHORT = "-p"

METADATA_KEY = "metadata"
METADATA_LONG = "--" + METADATA_KEY
METADATA_SHORT = "-r"

CONFIG_KEY = "config"
CONFIG_LONG = "--" + CONFIG_KEY
CONFIG_SHORT = "-c"

OUTDIR_KEY = "outdir"
OUTDIR_LONG = "--" + OUTDIR_KEY
OUTDIR_SHORT = "-o"

PARTITION_COLUMN_KEY = "partition_col"
PARTITION_COLUMN_LONG = "--" + PARTITION_COLUMN_KEY
PARTITION_COLUMN_SHORT = "-a"

ID_COLUMN_KEY = "id_col"
ID_COLUMN_LONG = "--" + ID_COLUMN_KEY
ID_COLUMN_SHORT = "-i"

OUTLIER_THRESHOLD_KEY = "outlier_thresh"
OUTLIER_THRESHOLD_LONG = "--" + OUTLIER_THRESHOLD_KEY

MINIMUM_MEMBERS_KEY = "min_members"
MINIMUM_MEMBERS_LONG= "--" + MINIMUM_MEMBERS_KEY
MINIMUM_MEMBERS_SHORT = "-m"

COUNT_MISSING_KEY = "count_missing"
COUNT_MISSING_LONG = "--" + COUNT_MISSING_KEY
COUNT_MISSING_SHORT = "-n"

MISSING_THRESHOLD_KEY = "missing_thresh"
MISSING_THRESHOLD_LONG = "--" + MISSING_THRESHOLD_KEY

DISTANCE_METHOD_KEY = "distm"
DISTANCE_METHOD_LONG = "--" + DISTANCE_METHOD_KEY

SKIP_QC_KEY = "skip_qc"
SKIP_QC_LONG = "--" + SKIP_QC_KEY
SKIP_QC_SHORT = "-s"

THRESHOLDS_KEY = "thresholds"
THRESHOLDS_LONG = "--" + THRESHOLDS_KEY
THRESHOLDS_SHORT = "-t"

DELIMITER_KEY = "delimiter"
DELIMITER_LONG = "--" + DELIMITER_KEY
DELIMITER_SHORT = "-d"

CLUSTER_METHOD_KEY = "method"
CLUSTER_METHOD_LONG = "--" + CLUSTER_METHOD_KEY
CLUSTER_METHOD_SHORT = "-e"

TREE_DISTANCES_KEY = "tree_distances"
TREE_DISTANCES_LONG = "--" + TREE_DISTANCES_KEY

FORCE_KEY = "force"
FORCE_LONG = "--" + FORCE_KEY
FORCE_SHORT = "-f"

SORT_MATRIX_KEY = "sort_matrix"
SORT_MATRIX_LONG = "--" + SORT_MATRIX_KEY

THREADS_KEY = "n_threads"
THREADS_LONG = "--" + THREADS_KEY

VERSION_KEY = "version"
VERSION_LONG = "--" + VERSION_KEY
VERSION_SHORT = "-V"

ONLY_REPORT_LABELED_KEY = "only_report_labeled_columns"
ONLY_REPORT_LABELED_LONG = "--" + ONLY_REPORT_LABELED_KEY

GROUPED_METADATA_COLUMNS_KEY = "grouped_metadata_columns"
LINELIST_COLUMNS_KEY = "linelist_columns"
DISPLAY_KEY = "display"
LABEL_KEY = "label"
GAS_CLUSTER_ADDRESS_KEY = "gas_denovo_cluster_address"

METADATA_INCLUDED_FILEPATH_TSV = "metadata.included.tsv"
METADATA_INCLUDED_FILEPATH_EXCEL = "metadata.included.xlsx"
METADATA_INCLUDED_SHEET_NAME = "Included Metadata"

CLUSTER_SUMMARY_FILEPATH_TSV = "cluster_summary.tsv"
CLUSTER_SUMMARY_FILEPATH_EXCEL = "cluster_summary.xlsx"
CLUSTER_SUMMARY_SHEET_NAME = "Cluster Summary"

PARAMETER_KEYS = [PROFILE_KEY, METADATA_KEY, CONFIG_KEY, OUTDIR_KEY,
                  PARTITION_COLUMN_KEY, ID_COLUMN_KEY, OUTLIER_THRESHOLD_KEY,
                  MINIMUM_MEMBERS_KEY, COUNT_MISSING_KEY, MISSING_THRESHOLD_KEY,
                  DISTANCE_METHOD_KEY, SKIP_QC_KEY, THRESHOLDS_KEY,
                  DELIMITER_KEY, CLUSTER_METHOD_KEY, TREE_DISTANCES_KEY,
                  FORCE_KEY, SORT_MATRIX_KEY, THREADS_KEY, VERSION_KEY,
                  ONLY_REPORT_LABELED_KEY, GROUPED_METADATA_COLUMNS_KEY,
                  LINELIST_COLUMNS_KEY]

BOOLEAN_KEYS = [COUNT_MISSING_KEY, SKIP_QC_KEY, FORCE_KEY, SORT_MATRIX_KEY, ONLY_REPORT_LABELED_KEY]

# Expected to check lowercase:
TRUE_STRINGS = ["t", "true"]

# Expected to check lowercase:
FALSE_STRINGS = ["f", "false"]

def parse_args():
    """ Argument Parsing method.

        A function to parse the command line arguments passed at initialization of Clade-o-matic,
        format these arguments,  and return help prompts to the user shell when specified.

        Returns
        -------
        ArgumentParser object
            The arguments and their user specifications, the usage help prompts and the correct formatting
            for the incoming argument (str, int, etc.)
        """
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        """
                Class to instantiate the formatter classes required for the argument parser.
                Required for the correct formatting of the default parser values

                Parameters
                ----------
                ArgumentDefaultsHelpFormatter object
                    Instatiates the default values for the ArgumentParser for display on the command line.
                RawDescriptionHelpFormatter object
                    Ensures the correct display of the default values for the ArgumentParser
                """
        pass

    parser = ArgumentParser(
        description="Arborator, an aggregate tool for producing summary reports of genetic distances within groups v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument(PROFILE_LONG, PROFILE_SHORT, type=str, required=True, help='Allelic profiles')
    parser.add_argument(METADATA_LONG, METADATA_SHORT, type=str, required=True, help='Matched metadata for samples in the allele profile')
    parser.add_argument(CONFIG_LONG, CONFIG_SHORT, type=str, required=False,
                        help='Configuration json')
    parser.add_argument(OUTDIR_LONG, OUTDIR_SHORT, type=str, required=True, help='Result output files')
    parser.add_argument(PARTITION_COLUMN_LONG, PARTITION_COLUMN_SHORT, type=str, required=False, help='Metadata column name for aggregating samples' )
    parser.add_argument(ID_COLUMN_LONG, ID_COLUMN_SHORT, type=str, required=False, help='Sample identifier column' )
    parser.add_argument(OUTLIER_THRESHOLD_LONG, type=float, required=False, help='Threshold to flag outlier comparisons within a group',default=100)
    parser.add_argument(MINIMUM_MEMBERS_LONG, MINIMUM_MEMBERS_SHORT, type=int, required=False,
                        help='Minimum number of members to perform clustering',default=2)
    parser.add_argument(ONLY_REPORT_LABELED_LONG, required=False, help='Only report labeled columns',
                        action='store_true')

    #profile dists
    parser.add_argument(COUNT_MISSING_LONG, COUNT_MISSING_SHORT, required=False, help='UNUSED: Count missing as differences',
                        action='store_true')
    parser.add_argument(MISSING_THRESHOLD_LONG, type=float, required=False,
                        help='UNUSED: Maximum percentage of missing data allowed per locus (0 - 1)')
    parser.add_argument(DISTANCE_METHOD_LONG, type=str, required=False, help='UNUSED: Distance method raw hamming or scaled difference [hamming, scaled]')
    parser.add_argument(SKIP_QC_LONG, SKIP_QC_SHORT, required=False, help='UNUSED: Skip QA/QC steps',
                        action='store_true')
    #GAS
    parser.add_argument(THRESHOLDS_LONG, THRESHOLDS_SHORT, type=str, required=False, help='thresholds delimited by ,',default='100')
    parser.add_argument(DELIMITER_LONG, DELIMITER_SHORT, type=str, required=False, help='UNUSED: delimiter desired for nomenclature code')
    parser.add_argument(CLUSTER_METHOD_LONG, CLUSTER_METHOD_SHORT, type=str, required=False, help='cluster method [single, complete, average]',
                        default='average')
    parser.add_argument(TREE_DISTANCES_LONG, type=str, required=False, default='patristic', choices=multi_level_clustering.VALID_TREE_DISTANCES,
                        help=('Defines how distances in distance matrices are interpretted by GAS and represented in the output tree (Newick file). '
                             'Use "patristic" to interpret distances in the matrix as sum of branch lengths between clusters or leaves, '
                             'and "cophenetic" to interpret distances in the matrix as the minimum distance two clusters or leaves need '
                             'to be in order to be grouped into the same cluster.'))

    parser.add_argument(FORCE_LONG, FORCE_SHORT, required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument(SORT_MATRIX_LONG, required=False,
                        help=('Sorts the samples in the distance matrix generated by GAS. The order of sample rarely '
                             'has an effect on the assigned cluster labels and sorting them ensures the same inputs always generate the same outputs.'),
                        action='store_true')
    parser.add_argument(THREADS_LONG, type=int, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument(VERSION_LONG, VERSION_SHORT, action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def remove_columns(df,missing_value,max_missing_frac=1):
    if max_missing_frac != 1:
        columns = list(df.columns)
        columns_to_remove = []
        num_records = len(df)
        for col in columns:
            unique_values = dict(df[col].astype(str).value_counts())
            if missing_value in unique_values:
                n = unique_values[missing_value]
                frac = n / num_records
                if frac > max_missing_frac:
                    columns_to_remove.append(col)

        return df.drop(columns_to_remove, axis=1)
    else:
        columns_to_remove = []
        return df.drop(columns_to_remove, axis=1)

def get_pairwise_outliers(distance_matrix, thresh):
    # Upper triangle of matrix to avoid duplicates:
    pairwise_distances = distance_matrix.where(np.triu(distance_matrix, k=1).astype(bool)).stack()

    pairwise_outliers = pairwise_distances[pairwise_distances.abs().gt(thresh)]
    pairwise_outliers_list = []

    for outlier in pairwise_outliers.items():
        # Every item in the series will look like:
        # [[index1, index2], value]
        indices = list(outlier[0])
        value = [outlier[1]]

        # New format: [index1, index2, value]
        outlier = indices + value

        pairwise_outliers_list.append(outlier)

    return pairwise_outliers_list

def get_average_outliers(distance_matrix, thresh):
    num_samples = len(distance_matrix)
    # Dividing by num_samples - 1, because the distance matrix
    # includes the distance of each sample to itself (0):
    averages = distance_matrix.sum(axis=1) / (num_samples - 1)
    outliers = averages[averages.abs().gt(thresh)]
    average_outliers_list = list(outliers.index)

    return average_outliers_list

def get_outliers(file_path, thresh, delim="\t"):
    PROFILE_DISTS_ID_INDEX = "dists" # This is not exposed in profile_dists.
    distance_matrix = pd.read_csv(file_path, sep=delim)
    distance_matrix = distance_matrix.set_index(PROFILE_DISTS_ID_INDEX)

    average_outliers_list = get_average_outliers(distance_matrix, thresh)
    pairwise_outliers_list = get_pairwise_outliers(distance_matrix, thresh)

    return (average_outliers_list, pairwise_outliers_list)

def write_outliers(outliers,outfile):
    with open(outfile, 'w') as f:
        f.write("id1\tid2\tdist\n")
        for row in outliers:
            f.write("{}\n".format("\t".join([str(x) for x in row])))

def stage_data(groups, outdir, metadata_df, id_col, group_file_mapping, max_missing_frac=1):
    files = {}
    for group_id in groups:
        directory_name = group_file_mapping[group_id]
        directory_path = os.path.join(outdir,f"{directory_name}")

        if not os.path.isdir(directory_path):
            os.makedirs(directory_path, 0o755)

        files[group_id] = {
            "profile": os.path.join(directory_path, "profile.tsv"),
            "parquet_matrix": os.path.join(directory_path, "matrix.pq"),
            "matrix": os.path.join(directory_path, "matrix.tsv"),
            "clusters": os.path.join(directory_path, "clusters.tsv"),
            "metadata": os.path.join(directory_path, "metadata.tsv"),
            "tree": os.path.join(directory_path, "tree.nwk"),
            "summary": os.path.join(directory_path, "loci.summary.tsv"),
            "outliers": os.path.join(directory_path, "outliers.tsv"),

        }

        #remove existing files if they exist
        for fname in files[group_id]:
            if os.path.isfile(files[group_id][fname]):
                os.remove(files[group_id][fname])

        df = remove_columns(groups[group_id], '0', max_missing_frac=max_missing_frac)
        df.to_csv(files[group_id][PROFILE_KEY], sep="\t", header=True, index=False)
        metadata_df[metadata_df[id_col].isin(list(groups[group_id][id_col]))].to_csv(files[group_id]['metadata'], sep="\t", header=True, index=False)

    return files

def process_data(group_files, id_col, group_col, thresholds, outlier_thresh, method, min_members,
                 tree_distance_representation, sort_matrix, num_cpus=1):
    try:
        sys_num_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        sys_num_cpus = cpu_count()

    if num_cpus > sys_num_cpus:
        num_cpus = sys_num_cpus

    pool = Pool(processes=num_cpus)

    results = []
    for group_id in group_files:
        results.append(pool.apply_async(process_group, (group_id, group_files[group_id], id_col, group_col, thresholds,
                                                        outlier_thresh, method, tree_distance_representation, sort_matrix,
                                                        min_members)))

    pool.close()
    pool.join()

    r = []
    for x in results:
        if isinstance(x,dict):
            r.append(x)
        else:
            r.append(x.get())

    return r

def process_group(group_id, output_files, id_col, group_col, thresholds,
                  outlier_thresh, method, tree_distance_representation,
                  sort_matrix, min_members=2):
    (allele_map, df) = process_profile(output_files[PROFILE_KEY], column_mapping={})
    l, p = convert_profiles(df)
    min_dist = 0
    mean_dist = 0
    med_dist = 0
    max_dist = 0
    outliers = {}
    outlier_ids = []
    metadata_summary = report(read_data(output_files[METADATA_KEY]).df,[id_col,group_col]).get_data()

    if len(l) >= min_members:
        # compute distances
        calc_distances_hamming(p, l, p, l, output_files['parquet_matrix'], len(l))
        pq.ParquetFile(output_files['parquet_matrix']).to_pandas().to_csv(output_files['matrix'], index=False, header=True,
                                                                          sep="\t")

        # perform clustering
        mc = multi_level_clustering(output_files['matrix'], thresholds, method, sort_matrix, tree_distances=tree_distance_representation)
        memberships = mc.get_memberships()
        with open(output_files['tree'], 'w') as fh:
            fh.write(f"{mc.newick}\n")

        write_clusters(memberships, len(thresholds), output_files['clusters'], ".")
        labels, dists = mc.read_distance_matrix(output_files['matrix'])
        min_dist = min(dists)
        mean_dist = mean(dists)
        med_dist = median(dists)
        max_dist = max(dists)
        report(df, [id_col]).write_data(output_files['summary'])
        (outlier_ids, pairwise_outlier) = get_outliers(output_files['matrix'], outlier_thresh, delim="\t")
        write_outliers(pairwise_outlier, output_files['outliers'])

        if os.path.isfile(output_files['clusters']) and os.path.isfile(output_files["metadata"]):
            clust_df = pd.read_csv(output_files['clusters'], sep="\t", header=0, dtype=str)
            clust_df = clust_df[['id','address']]
            clust_df['address'] = str(group_id) + "|" + clust_df['address'].astype(str) # appends "{group_id}|" to the address
            clust_df = clust_df.rename(columns={'id': id_col,'address':GAS_CLUSTER_ADDRESS_KEY})
            clust_df.to_csv(output_files['clusters'],header=True,sep="\t",index=False)
            metadata_df = pd.read_csv(output_files[METADATA_KEY], sep="\t", header=0, dtype=str)
            pd.merge(metadata_df, clust_df, on=id_col).to_csv(output_files[METADATA_KEY],sep="\t",header=True,index=False)
            del(clust_df)
            del(metadata_df)

    return { group_id:{
        'count_members': len(l),
        'min_dist': min_dist,
        'mean_dist': mean_dist,
        'median_dist': med_dist,
        'max_dist': max_dist,
        'count_outliers': len(outlier_ids),
        'outlier_ids':",".join([str(x) for x in outlier_ids]),
        'metadata':metadata_summary
    }
}

def compile_group_data(group_metrics, field_data_types,id_col,field_name_key,field_name_value,header=[]):
    s = summarizer(header,group_metrics,field_data_types,field_name_key,field_name_value)
    data = s.get_data()
    header = s.header
    for id in data:
        record = data[id]
        record[id_col] = id
        for k in group_metrics[id]:
            if k == 'metadata':
                continue
            data[id][k] = str(group_metrics[id][k])
        for k in record:
            value = record[k]
            if isinstance(value, datetime):
                data[id][k] = value.strftime('%Y-%m-%d')
            else:
                data[id][k] = str(value)
    return data

def format_df(column_map,df):
    df_cols = list(df.columns)
    cols_to_remove = list(set(df_cols) - set(column_map.keys()))
    df = df.drop(cols_to_remove, axis=1)
    df = df[list(column_map.keys())]
    return df.rename(columns=column_map)

def prepare_column_map(column_info,columns):
    column_map = {}
    for f in column_info:
        column_map[f] = column_info[f][LABEL_KEY]
        if column_map[f] == "":
            column_map[f] = f

    for f in columns:
        column_map[f] = f

    return column_map

def prepare_linelist(column_info,df,columns):
    if len(columns) == 0:
        columns = df.columns.to_list()
    column_map = prepare_column_map(column_info,columns)
    return format_df(column_map, df)


def update_column_order(df,col_properties,restrict=False):
    cols = {}
    df_cols = list(df.columns)
    num_rows = len(df)
    for col in col_properties:
        if DISPLAY_KEY in col_properties[col] and col_properties[col][DISPLAY_KEY]:
            if LABEL_KEY in col_properties[col]:
                cols[col] = col_properties[col][LABEL_KEY]
            else:
                cols[col] = str(col)
            if col not in df_cols:
                df[col] = [col_properties[col]['default']] * num_rows

    order = list(cols.keys())
    if not restrict:
        to_add = set(df_cols ) - set(order)
        for c in df_cols:
            if c in to_add:
                order.append(c)
                cols[c] = c

    df = df.reindex(columns=order)
    df = df.rename(columns=cols)

    return df[list(cols.values())]

def convert_to_bool(input):

    if not isinstance(input, bool):
            if str(input).lower() in TRUE_STRINGS:
                result = True
            elif str(input).lower() in FALSE_STRINGS:
                result = False
            else:
                message = f'Expected a boolean-like string but found {input}'
                raise Exception(message)
    else:
        result = input

    return result

def validate_params(config):
    params = [PROFILE_KEY, METADATA_KEY, OUTDIR_KEY, ID_COLUMN_KEY, PARTITION_COLUMN_KEY, MINIMUM_MEMBERS_KEY]
    missing = []
    for p in params:
        if p not in config or config[p] == '' or config[p] == None:
            missing.append(p)
    if len(missing) > 0:
        message = f"Error, parameters not set for: {missing}"
        raise Exception(message)

    for key in config.keys():
        # Check for any unexpected config keys:
        if key not in PARAMETER_KEYS:
            print(f'WARNING: "{key}" parameter unrecognized')

        # Convert string booleans into actual booleans:
        if key in BOOLEAN_KEYS:
            config[key] = convert_to_bool(config[key])

    # Convert string booleans into actual booleans:
    if GROUPED_METADATA_COLUMNS_KEY in config:
        summaries = config[GROUPED_METADATA_COLUMNS_KEY]
        for summary in summaries:
            if DISPLAY_KEY in summaries[summary]:
                display = summaries[summary][DISPLAY_KEY]
                summaries[summary][DISPLAY_KEY] = convert_to_bool(display)

    # Convert string booleans into actual booleans:
    if LINELIST_COLUMNS_KEY in config:
        summaries = config[LINELIST_COLUMNS_KEY]
        for summary in summaries:
            if DISPLAY_KEY in summaries[summary]:
                display = summaries[summary][DISPLAY_KEY]
                summaries[summary][DISPLAY_KEY] = convert_to_bool(display)

def cluster_reporter(config):
    validate_params(config)
    profile_file = config[PROFILE_KEY]
    partition_file = config[METADATA_KEY]
    outdir = config[OUTDIR_KEY]
    outlier_thresh = config[OUTLIER_THRESHOLD_KEY]
    thresholds = config[THRESHOLDS_KEY]
    method = config[CLUSTER_METHOD_KEY]
    tree_distance_representation = config[TREE_DISTANCES_KEY]
    force = config[FORCE_KEY]
    sort_matrix = config[SORT_MATRIX_KEY]
    id_col = config[ID_COLUMN_KEY]
    partition_col = config[PARTITION_COLUMN_KEY]
    min_members = config[MINIMUM_MEMBERS_KEY]
    num_threads = config[THREADS_KEY]
    restrict_output = config[ONLY_REPORT_LABELED_KEY]

    # Unused parameters:
    skip_qc = config[SKIP_QC_KEY]
    missing_thresh = config[MISSING_THRESHOLD_KEY]
    distm = config[DISTANCE_METHOD_KEY]
    count_missing = config[COUNT_MISSING_KEY]
    delimiter = config[DELIMITER_KEY]

    # We're leaving the skip_qc for later, but want to warn.
    # Since it's in argparse as a flag, it will always be false
    # if not provided.
    if(skip_qc):
        print(f'WARNING: skip QC ({SKIP_QC_LONG}/{SKIP_QC_SHORT}) was provided, but this parameter is currently unused.')

    if(missing_thresh):
        print(f'WARNING: missing threshold ({MISSING_THRESHOLD_LONG}) was provided, but this parameter is currently unused.')

    if(distm):
        print(f'WARNING: distance method ({DISTANCE_METHOD_LONG}) was provided, but this parameter is currently unused.')

    # See above comment for skip_qc.
    if(count_missing):
        print(f'WARNING: count missing ({COUNT_MISSING_LONG}/{COUNT_MISSING_SHORT}) was provided, but this parameter is currently unused.')

    if(delimiter):
        print(f'WARNING: delimiter ({DELIMITER_LONG}/{DELIMITER_SHORT}) was provided, but this parameter is currently unused.')

    try:
        sys_num_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        sys_num_cpus = cpu_count()

    if num_threads < 1:
        message = f'{THREADS_KEY} ({num_threads}) needs to be at least 1.'
        raise Exception(message)
    elif num_threads > sys_num_cpus:
        print(f'WARNING: {THREADS_KEY} ({num_threads}) exceeds the number of CPUs available ({sys_num_cpus}). Setting {THREADS_KEY} to {sys_num_cpus}.')
        num_threads = sys_num_cpus

    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config

    linelist_cols_properties = {}
    line_list_columns = []
    if LINELIST_COLUMNS_KEY in config:
        linelist_cols_properties = config[LINELIST_COLUMNS_KEY]
        for f in linelist_cols_properties:
            if DISPLAY_KEY in linelist_cols_properties[f]:
                if linelist_cols_properties[f][DISPLAY_KEY]:
                    line_list_columns.append(f)

    cluster_summary_cols_properties = {}
    cluster_summary_header = []
    cluster_display_cols_to_remove = []
    if GROUPED_METADATA_COLUMNS_KEY in config:
        cluster_summary_cols_properties = config[GROUPED_METADATA_COLUMNS_KEY]
        cluster_summary_header = list(cluster_summary_cols_properties.keys())
        for f in cluster_summary_cols_properties:
            cluster_summary_cols_properties[f]['data_type'] = cluster_summary_cols_properties[f]['data_type'].lower()
            if DISPLAY_KEY in cluster_summary_cols_properties[f]:
                if not cluster_summary_cols_properties[f][DISPLAY_KEY]:
                    cluster_display_cols_to_remove.append(f)
    else:
        cluster_summary_cols_properties[partition_col] = { "data_type": "none","label":partition_col,"default":"","display":"True"}
        cluster_summary_cols_properties['min_dist'] = { "data_type": "none","label":'min_dist',"default":"","display":"True"}
        cluster_summary_cols_properties['median_dist'] = { "data_type": "none","label":'median_dist',"default":"","display":"True"}
        cluster_summary_cols_properties['mean_dist'] = { "data_type": "none","label":'mean_dist',"default":"","display":"True"}
        cluster_summary_cols_properties['max_dist'] = { "data_type": "none","label":'max_dist',"default":"","display":"True"}
        cluster_summary_cols_properties['count_outliers'] = { "data_type": "none","label":'count_outliers',"default":"","display":"True"}
        cluster_summary_cols_properties['outlier_ids'] = { "data_type": "none","label":'outlier_ids',"default":"","display":"True"}
        cluster_summary_header = list(cluster_summary_cols_properties.keys())

    if len(cluster_summary_header) == 0:
        cluster_summary_header = [partition_col]
        display_cluster_header = [partition_col]

    if not os.path.isfile(profile_file):
        message = f'Profile path {profile_file} does not exist, please check path and try again'
        raise Exception(message)

    if not os.path.isfile(partition_file):
        message = f'Metadata file {partition_file} does not exist, please check path and try again'
        raise Exception(message)

    if not isinstance(outlier_thresh,int) or not isinstance(outlier_thresh,float):
        try:
            outlier_thresh = float(outlier_thresh)
        except:
            message = f'Outlier threshold needs to be numeric: {outlier_thresh}'
            raise Exception(message)

    if not isinstance(thresholds,list):
        thresholds = thresholds.split(',')

    thresholds = process_thresholds(thresholds)

    if not method in CLUSTER_METHODS:
        message = f'Linkage method supplied is invalid: {method}, it needs to be one of average, single, complete'
        raise Exception(message)

    if not isinstance(min_members, int):
        try:
            min_members = int(min_members)
        except:
            message = f'Min members needs to be an integer {min_members}'
            raise Exception(message)

    if min_members < 2:
        message = f'{MINIMUM_MEMBERS_KEY} ({min_members}) needs to be at least 2.'
        raise Exception(message)

    if not force and os.path.isdir(outdir):
        message = f'folder {outdir} already exists, please choose new directory or use --force'
        raise Exception(message)

    # initialize analysis directory
    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    (allele_map, profile_df) = process_profile(profile_file, column_mapping={})
    profile_df.insert(0, id_col, profile_df.index.to_list())

    #write allele mapping file
    with open(os.path.join(outdir,"allele_map.json"),'w' ) as fh:
        fh.write(json.dumps(allele_map, indent=4))

    metadata = read_data(partition_file)
    metadata_df = metadata.df

    if len(metadata_df) == 0:
        message = f'No metadata rows were provided.'
        raise Exception(message)

    input_profile_samples = set(profile_df[id_col])
    input_metadata_samples = set(metadata_df[id_col])

    ovl_samples = input_profile_samples & input_metadata_samples
    missing_profile_samples = input_profile_samples - ovl_samples
    missing_metadata_samples = input_metadata_samples - ovl_samples

    run_data['count_missing_profile_samples'] = len(missing_profile_samples)
    run_data['missing_profile_samples'] = ",".join(sorted(list(missing_profile_samples)))
    run_data['count_missing_metadata_samples'] = len(missing_metadata_samples)
    run_data['missing_metadata_samples'] = ",".join(sorted(list(missing_metadata_samples)))
    ovl_samples = list(ovl_samples)

    metadata_df[metadata_df[id_col].isin(ovl_samples)].to_csv(os.path.join(outdir,"metadata.overlap.tsv"),sep="\t",header=True,index=False)
    split = split_profiles(profile_df[profile_df[id_col].isin(ovl_samples)],os.path.join(outdir,"metadata.overlap.tsv"),id_col,partition_col)
    groups = split.subsets
    group_file_mapping = split.group_file_mapping

    filtered_samples = pd.concat(list(groups.values()), ignore_index=True)[id_col].to_list()
    linelist_df = prepare_linelist({}, metadata_df[metadata_df[id_col].isin(filtered_samples)], columns=[])
    ll_cols = list(set(linelist_df.columns.to_list()))
    if restrict_output:
        t = []
        for c in line_list_columns:
            if c in ll_cols:
                t.append(c)
        line_list_columns = t
    else:
        t = []
        for c in line_list_columns:
            if c in ll_cols:
                t.append(c)
        for c in ll_cols:
            if not c in t:
                t.append(c)
        line_list_columns = t

    linelist_df = prepare_linelist({}, metadata_df[metadata_df[id_col].isin(list(set(metadata_df[id_col].to_list()) - set(filtered_samples)))], columns=[])
    linelist_df = linelist_df[line_list_columns]
    linelist_df.to_csv(os.path.join(outdir, "metadata.excluded.tsv"), sep="\t", header=True, index=False)
    del(linelist_df)

    run_data['threshold_map'] = format_threshold_map(thresholds)
    with open(os.path.join(outdir,"threshold_map.json"),'w' ) as fh:
        fh.write(json.dumps(run_data['threshold_map'], indent=4))

    group_files = stage_data(groups, outdir, metadata_df, id_col, group_file_mapping, max_missing_frac=1)
    results = process_data(group_files, id_col, partition_col, thresholds, outlier_thresh, method, min_members, tree_distance_representation, sort_matrix, num_cpus=num_threads)
    group_metrics = {}
    for r in results:
        for k in r:
            group_metrics[k] = r[k]

    #merge metadata files

    summary_file = os.path.join(outdir, CLUSTER_SUMMARY_FILEPATH_TSV)


    summary_data = compile_group_data(group_metrics=group_metrics, field_data_types=cluster_summary_cols_properties,
                                      id_col=partition_col, field_name_key='metadata', field_name_value='value_counts',
                                      header=cluster_summary_header)
    summary_df = pd.DataFrame.from_dict(summary_data, orient='index')
    cluster_display_cols_to_remove = list(set(cluster_display_cols_to_remove) & set(list(summary_df.columns)))
    summary_df = summary_df.drop(cluster_display_cols_to_remove, axis=1)
    summary_cols = sorted(list(summary_df.columns))
    display_columns = []
    for col in cluster_summary_cols_properties:
        prop = cluster_summary_cols_properties[col]
        if DISPLAY_KEY in prop:
            if prop[DISPLAY_KEY]:
                display_columns.append(col)
        else:
            display_columns.append(col)
    for col in summary_cols:
        if col in display_columns:
            continue
        display_columns.append(col)

    summary_df = summary_df[display_columns]
    for k in cluster_display_cols_to_remove:
        del(cluster_summary_cols_properties[k])
    summary_df = update_column_order(summary_df, cluster_summary_cols_properties, restrict=restrict_output)
    summary_df.to_csv(summary_file, sep="\t", index=False, header=True)
    summary_df.to_excel(os.path.join(outdir, CLUSTER_SUMMARY_FILEPATH_EXCEL), header=True, index=False, sheet_name=CLUSTER_SUMMARY_SHEET_NAME)
    
    if LINELIST_COLUMNS_KEY in config:
        line_list_columns = []
        linelist_cols_properties = config[LINELIST_COLUMNS_KEY]
        for f in linelist_cols_properties:
            if DISPLAY_KEY in linelist_cols_properties[f]:
                if linelist_cols_properties[f][DISPLAY_KEY]:
                    line_list_columns.append(f)

    if not restrict_output and GAS_CLUSTER_ADDRESS_KEY not in line_list_columns:
        line_list_columns.append(GAS_CLUSTER_ADDRESS_KEY)

    metadata_dfs = []
    for group_id in group_files:
        #remove parquet file
        num_members = 0
        f = group_files[group_id]['parquet_matrix']

        if os.path.isfile(f):
            os.remove(f)
        f = group_files[group_id]["metadata"]

        if os.path.isfile(f):
            obj = read_data(f)

            if obj.status:
                num_members = len(obj.df)
                metadata_dfs.append(obj.df)

        if num_members < min_members:
            directory_name = group_file_mapping[group_id]
            shutil.rmtree(os.path.join(outdir, directory_name))

    linelist_df = pd.concat(metadata_dfs, ignore_index=True, sort=False)

    # Only try to load metadata columns that actually exists:
    intersection = set(line_list_columns).intersection(set(linelist_df.columns))

    # Ensure clustering was successful and therefore the GAS_CLUSTER_ADDRESS_KEY
    # column exists in both the line list and dataframe:
    if GAS_CLUSTER_ADDRESS_KEY in intersection:

        # Warn about metadata columns specified in the line list that don't exist
        # in the metadata. This warning is inside the conditional, otherwise
        # it will report that GAS_CLUSTER_ADDRESS_KEY was specified in the
        # line list, but doesn't exist, which isn't true. It's how arborator handles
        # this data.
        difference = set(line_list_columns).difference(set(linelist_df.columns))

        for item in difference:
            print(f'WARNING: "{item}" specified in the line list, but does not exist in the metadata.')

        linelist_df = linelist_df[list(intersection)]
        linelist_df = update_column_order(linelist_df, linelist_cols_properties, restrict=restrict_output)

        linelist_df.to_csv(os.path.join(outdir, METADATA_INCLUDED_FILEPATH_TSV), sep="\t", header=True, index=False)
        linelist_df.to_excel(os.path.join(outdir, METADATA_INCLUDED_FILEPATH_EXCEL), header=True, index=False, sheet_name=METADATA_INCLUDED_SHEET_NAME)

    else:
        print(f'WARNING: Failed to generate any clusters! No "{METADATA_INCLUDED_FILEPATH_TSV}" will be generated.')

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    sys.stdout.flush()

    with open(os.path.join(outdir, "run.json"), 'w') as fh:
        fh.write(json.dumps(run_data, indent=4))

def process_thresholds(thresholds):

    try:
        processed = [float(x) for x in thresholds]
    except ValueError:
        message = f'thresholds {thresholds} must all be integers or floats'
        raise Exception(message)

    # Thresholds must be strictly decreasing:
    if not all(processed[i] > processed[i+1] for i in range(len(processed)-1)):
        message = f'thresholds {thresholds} must be in decreasing order'
        raise Exception(message)

    # Thresholds must be non-negative:
    if not all(processed[i] >= 0 for i in range(len(processed))):
        message = f'thresholds {thresholds} must be non-negative'
        raise Exception(message)

    return processed

def main():
    cmd_args = parse_args()
    config_file = cmd_args.config

    # Initialize based on argparse (command-line arguments):
    config = vars(cmd_args)

    # Overwrite with config file parameters:
    if config_file is not None:

        if not os.path.isfile(config_file):
            message = f'Config path {config_file} does not exist, please check path and try again'
            raise Exception(message)

        with open(config_file) as fh:
            c = json.loads(fh.read())
            for field in c:
                config[field] = c[field]

    if not OUTLIER_THRESHOLD_KEY in config or config[OUTLIER_THRESHOLD_KEY] == '':
        message = f'Error you must supply an outlier threshold as a cmd line parameter or in the config file'
        raise Exception(message)

    if not THRESHOLDS_KEY in config or config[THRESHOLDS_KEY] == '':
        message = f'Error you must supply a threshold as a cmd line parameter or in the config file'
        raise Exception(message)

    cluster_reporter(config)


# call main function
if __name__ == '__main__':
    main()
