import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import json
import os
from datetime import datetime
import pandas as pd
from statistics import mean, median
import shutil
from arborator.version import __version__
from arborator.classes.aggregator import summarizer
from profile_dists.utils import convert_profiles, calc_distances_hamming, process_profile
from arborator.classes.read_data import read_data
from arborator.classes.report import report
from arborator.classes.split_profiles import split_profiles
from genomic_address_service.classes.multi_level_clustering import multi_level_clustering
from genomic_address_service.utils import format_threshold_map, write_threshold_map
from genomic_address_service.mcluster import write_clusters
import fastparquet as pq
from multiprocessing import Pool, cpu_count

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
    parser.add_argument('--profile','-p', type=str, required=True, help='Allelic profiles')
    parser.add_argument('--metadata','-r', type=str, required=True, help='Matched metadata for samples in the allele profile')
    parser.add_argument('--config', '-c', type=str, required=False,
                        help='Configuration json')
    parser.add_argument('--outdir', '-o', type=str, required=True, help='Result output files')
    parser.add_argument('--partition_col', '-a', type=str, required=False, help='Metadata column name for aggregating samples' )
    parser.add_argument('--id_col', '-i', type=str, required=False, help='Sample identifier column' )
    parser.add_argument('--outlier_thresh', type=float, required=False, help='Threshold to flag outlier comparisons within a group',default=100)
    parser.add_argument('--min_members','-m', type=int, required=False,
                        help='Minimum number of members to perform clustering',default=2)

    #profile dists
    parser.add_argument('-n', '--count_missing', required=False, help='Count missing as differences',
                        action='store_true')
    parser.add_argument('--missing_thresh', type=float, required=False,
                        help='Maximum percentage of missing data allowed per locus (0 - 1)',default=1.0)
    parser.add_argument('--distm', type=str, required=False, help='Distance method raw hamming or scaled difference [hamming, scaled]',default='scaled')
    parser.add_argument('-s', '--skip_qc', required=False, help='Skip QA/QC steps',
                        action='store_true')
    #GAS
    parser.add_argument('-t','--thresholds', type=str, required=False, help='thresholds delimited by ,',default='100')
    parser.add_argument('-d', '--delimeter', type=str, required=False, help='delimeter desired for nomenclature code',
                        default=".")
    parser.add_argument('-e', '--method', type=str, required=False, help='cluster method [single, complete, average]',
                        default='average')

    parser.add_argument('--force','-f', required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument('--n_threads', type=int, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def remove_columns(df,missing_value,max_missing_frac=1):
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

def get_outliers_matrix(file_path, thresh,delim="\t"):
    outliers = []
    ave_dists = {}
    with open(file_path, 'r') as f:
        header = next(f).strip().split(delim)[1:]
        offset = 2
        for idx,value in enumerate(header[1:]):
            ave_dists[value] = 0
        for line in f:
            
            line_split = line.strip().split(delim)

            if len(line_split) < 1:
                continue
            v = [float(x) for x in line_split[offset:]]
            h = header[offset-1:]
            for idx,value in enumerate(v):
                ave_dists[h[idx]]+=value
                if value > thresh:
                    outliers.append([line_split[0],h[idx],value])
            offset += 1
    num_samples = len(ave_dists)
    outlier_ids = list()
    for id in ave_dists:
        ave = ave_dists[id] / num_samples
        if ave > thresh:
            outlier_ids.append(id)

    return (outlier_ids,outliers)

def write_outliers(outliers,outfile):
    with open(outfile, 'w') as f:
        f.write("id1\tid2\tdist\n")
        for row in outliers:
            f.write("{}\n".format("\t".join([str(x) for x in row])))

def stage_data(groups,outdir,metadata_df,id_col,max_missing_frac=1):
    files = {}
    for group_id in groups:
        d = os.path.join(outdir,f"{group_id}")
        if not os.path.isdir(d):
            os.makedirs(d, 0o755)

        files[group_id] = {
            "profile": os.path.join(d, "profile.tsv"),
            "parquet_matrix": os.path.join(d, "matrix.pq"),
            "matrix": os.path.join(d, "matrix.tsv"),
            "clusters": os.path.join(d, "clusters.tsv"),
            "metadata": os.path.join(d, "metadata.tsv"),
            "tree": os.path.join(d, "tree.nwk"),
            "summary": os.path.join(d, "loci.summary.tsv"),
            "outliers": os.path.join(d, "outliers.tsv"),

        }

        #remove existing files if they exist
        for fname in files[group_id]:
            if os.path.isfile(files[group_id][fname]):
                os.remove(files[group_id][fname])

        df = remove_columns(groups[group_id], '0', max_missing_frac=max_missing_frac)
        df.to_csv(files[group_id]['profile'], sep="\t", header=True, index=False)
        metadata_df[metadata_df[id_col].isin(list(groups[group_id][id_col]))].to_csv(files[group_id]['metadata'], sep="\t", header=True, index=False)

    return files

def process_data(group_files,id_col,group_col,thresholds,outlier_thresh,method,num_cpus=1):
    try:
        sys_num_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        sys_num_cpus = cpu_count()

    if num_cpus > sys_num_cpus:
        num_cpus = sys_num_cpus

    pool = Pool(processes=num_cpus)

    results = []
    for group_id in group_files:
        results.append(pool.apply_async(process_group, (group_id,group_files[group_id],id_col,group_col,thresholds,outlier_thresh,method)))

    pool.close()
    pool.join()

    r = []
    for x in results:
        if isinstance(x,dict):
            r.append(x)
        else:
            r.append(x.get())

    return r

def process_group(group_id,output_files,id_col,group_col,thresholds,outlier_thresh,method,min_members=2):
    (allele_map, df) = process_profile(output_files['profile'], column_mapping={})
    l, p = convert_profiles(df)
    min_dist = 0
    mean_dist = 0
    med_dist = 0
    max_dist = 0
    outliers = {}
    outlier_ids = []
    metadata_summary = report(read_data(output_files['metadata']).df,[id_col,group_col]).get_data()

    if len(l) >= min_members:
        # compute distances
        calc_distances_hamming(p, l, p, l, output_files['parquet_matrix'], len(l))
        pq.ParquetFile(output_files['parquet_matrix']).to_pandas().to_csv(output_files['matrix'], index=False, header=True,
                                                                          sep="\t")

        # perform clustering
        mc = multi_level_clustering(output_files['matrix'], thresholds, method)
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
        (outlier_ids, pairwise_outlier) = get_outliers_matrix(output_files['matrix'], outlier_thresh, delim="\t")
        write_outliers(pairwise_outlier, output_files['outliers'])

        if os.path.isfile(output_files['clusters']) and os.path.isfile(output_files["metadata"]):
            clust_df = pd.read_csv(output_files['clusters'],sep="\t",header=0)
            clust_df = clust_df[['id','address']]
            values = clust_df['address']
            for idx,value in enumerate(values):
                values[idx] = f'{group_id}|{value}'
            clust_df['address'] = values
            clust_df = clust_df.rename(columns={'id': id_col,'address':'gas_denovo_cluster_address'})
            clust_df.to_csv(output_files['clusters'],header=True,sep="\t",index=False)     
            metadata_df = pd.read_csv(output_files['metadata'],sep="\t",header=0)
            pd.merge(metadata_df, clust_df, on=id_col).to_csv(output_files['metadata'],sep="\t",header=True,index=False)
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
        column_map[f] = column_info[f]['label']
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
        cols[col] = col_properties[col]['label']
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

def validate_params(config):
    params = ['profile','metadata','outdir','id_col','partition_col','min_members']
    missing = []
    for p in params:
        if p not in config or config[p] == '' or config[p] == None:
            missing.append(p)
    if len(missing) > 0:
        print(f"Error, parameters not set for: {missing}")
        sys.exit()

    


def cluster_reporter(config):
    validate_params(config)
    profile_file = config['profile']
    partition_file = config['metadata']
    outdir = config['outdir']
    outlier_thresh = config['outlier_thresh']
    thresholds = config['thresholds']
    method = config['method']
    force = config['force']
    id_col = config['id_col']
    partition_col = config['partition_col']
    min_members = config['min_members']
    skip_qc = config['skip_qc']
    num_threads = config['n_threads']

    missing_thresh = config['missing_thresh']
    distm = config['distm']
    count_missing = config['count_missing']
    delimeter = config['delimeter']

    try:
        sys_num_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        sys_num_cpus = cpu_count()

    if num_threads > sys_num_cpus:
        num_threads = sys_num_cpus

    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config
    restrict_output = False
    if "only_report_labeled_columns" in config:
        restrict_output = config["only_report_labeled_columns"]
        if not isinstance(restrict_output, bool):
            if restrict_output.lower() in ['f', 'false']:
                restrict_output = False
            elif restrict_output.lower() in ['t', 'true']:
                restrict_output = True
            else:
                print(f'only_report_labeled_columns needs to be true or false : you supplied {restrict_output}')
                sys.exit()


    linelist_cols_properties = {}
    line_list_columns = []
    if "linelist_columns" in config:
        linelist_cols_properties = config["linelist_columns"]
        for f in linelist_cols_properties:
            if 'display' in linelist_cols_properties[f]:
                v = linelist_cols_properties[f]['display'].lower()
                if v in ['t','true']:
                    line_list_columns.append(f)


    cluster_summary_cols_properties = {}
    cluster_summary_header = []
    cluster_display_cols_to_remove = []
    if "grouped_metadata_columns" in config:
        cluster_summary_cols_properties = config["grouped_metadata_columns"]
        cluster_summary_header = list(cluster_summary_cols_properties.keys())
        for f in cluster_summary_cols_properties:
            cluster_summary_cols_properties[f]['data_type'] = cluster_summary_cols_properties[f]['data_type'].lower()
            if 'display' in cluster_summary_cols_properties[f]:
                v = cluster_summary_cols_properties[f]['display'].lower()
                if v in ['f','false']:
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
        print(f'Profile path {profile_file} does not exist, please check path and try again')
        sys.exit()

    if not os.path.isfile(partition_file):
        print(f'Metadata file {partition_file} does not exist, please check path and try again')
        sys.exit()

    if not isinstance(outlier_thresh,int) or not isinstance(outlier_thresh,float):
        try:
            outlier_thresh = float(outlier_thresh)
        except:
            print(f'Outlier threshold needs to be numeric: {outlier_thresh}')
            sys.exit()

    if not isinstance(thresholds,list):
        thresholds = thresholds.split(',')

    try:
        thresholds = [float(x) for x in thresholds]
    except:
        print(f'Thresholds needs to be numeric and delineated by comma {thresholds}')
        sys.exit()

    if not isinstance(force, bool):
        if force.lower() in ['f','false']:
            force = False
        elif force.lower() in ['t','true']:
            force = True
        else:
            print(f'Force needs to be true or false : you supplied {force}')
            sys.exit()


    if not method in ['average','complete','single']:
        print(f'Linkage method supplied is invalid: {method}, it needs to be one of average, single, complete')
        sys.exit()

    if not isinstance(min_members, int):
        try:
            min_members = int(min_members)
        except:
            print(f'Min members needs to be an integer {min_members}')
            sys.exit()

    if not force and os.path.isdir(outdir):
        print(f'folder {outdir} already exists, please choose new directory or use --force')
        sys.exit()

    # initialize analysis directory
    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    (allele_map, profile_df) = process_profile(profile_file, column_mapping={})
    profile_df = profile_df.copy()
    profile_df.insert(0, id_col, profile_df.index.to_list())

    #write allele mapping file
    with open(os.path.join(outdir,"allele_map.json"),'w' ) as fh:
        fh.write(json.dumps(allele_map, indent=4))

    metadata = read_data(partition_file)
    metadata_df = metadata.df

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
    groups = split_profiles(profile_df[profile_df[id_col].isin(ovl_samples)],os.path.join(outdir,"metadata.overlap.tsv"),id_col,partition_col).subsets
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

    group_files = stage_data(groups, outdir, metadata_df, id_col, max_missing_frac=1)
    results = process_data(group_files, id_col, partition_col, thresholds, outlier_thresh, method, num_threads)
    group_metrics = {}
    for r in results:
        for k in r:
            group_metrics[k] = r[k]

    #merge metadata files

    summary_file = os.path.join(outdir, "cluster_summary.tsv")


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
        if 'display' in prop:
            if prop['display'] in ['true','True',True]:
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


    
    if "linelist_columns" in config:
        line_list_columns = []
        linelist_cols_properties = config["linelist_columns"]
        for f in linelist_cols_properties:
            if 'display' in linelist_cols_properties[f]:
                v = linelist_cols_properties[f]['display'].lower()
                if v in ['t','true']:
                    line_list_columns.append(f)

    if not restrict_output and 'gas_denovo_cluster_address' not in line_list_columns:
        line_list_columns.append('gas_denovo_cluster_address')
    

    metadata_dfs = []

    for group_id in group_files:
        #remove parquet file
        num_members = 0
        f = group_files[group_id]['parquet_matrix']
        if os.path.isfile(f):
            os.remove(f)
        f = group_files[group_id][ "metadata"]
        if os.path.isfile(f):
            obj = read_data(f)
            if obj.status:
                num_members = len(obj.df)
                metadata_dfs.append(obj.df)
        if num_members < min_members:
            shutil.rmtree(os.path.join(outdir,group_id))

    
    linelist_df = pd.concat(metadata_dfs, ignore_index=True, sort=False)
    linelist_df = linelist_df[line_list_columns]
    linelist_df.to_csv(os.path.join(outdir,"metadata.included.tsv"),sep="\t",header=True,index=False)

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    sys.stdout.flush()

    with open(os.path.join(outdir, "run.json"), 'w') as fh:
        fh.write(json.dumps(run_data, indent=4))


def main():
    cmd_args = parse_args()
    config_file = cmd_args.config
    config = vars(cmd_args)

    if config_file is not None:
        with open(config_file) as fh:
            c = json.loads(fh.read())
            for field in c:
                config[field] = c[field]
    

    if not 'outlier_thresh' in config or config['outlier_thresh'] == '':
        print(f'Error you must supply an outlier threshold as a cmd line parameter or in the config file')
        sys.exit()

    if not 'thresholds' in config or config['thresholds'] == '':
        print(f'Error you ust supply a threshold as a cmd line parameter or in the config file')
        sys.exit()


    cluster_reporter(config)




# call main function
if __name__ == '__main__':
    main()


