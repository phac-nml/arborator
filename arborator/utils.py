import os.path
import psutil
import pandas as pd
import numpy as np
import fastparquet as fp
from numba import jit
from numba.typed import List
import pyarrow.parquet as pq
import re
import sys
import gc
from profile_dists.constants import MIN_FILE_SIZE, FILE_FORMATS, VALID_INT_TYPES


def guess_format(unique_values):
    '''
    Function accepts a list of values and determines if the values are consistent with a list of integers,
    or a set of fixed length hash values
    :param unique_values: list of values
    :return: string (hash, mix, int)
    '''
    length_equal = is_all_same_len(unique_values)
    has_integers = contains_integers(unique_values)
    has_alpha  = contains_alpha(unique_values)

    format = ''
    #columns contains hash codes
    if length_equal and has_integers and has_alpha:
        format = 'hash'

    # columns contains a mix of integers and other info
    elif has_integers and has_alpha:
        format = 'mix'

    #columns contain only integers
    elif has_integers:
        format = 'int'

    return format


def is_all_same_len(unique_values):
    '''
    Accepts a list of values and calculates the length of each one
    :param unique_values: list of numeric or string values
    :return: True if all values are the same length
    '''
    l = set()
    for idx,value in enumerate(unique_values):
        if value != '0':
            l.add(len(str(value)))
    if len(l) == 1:
        status = True
    else:
        status = False
    return status


def contains_integers(unique_values):
    '''
    Accepts a list of values and determines if any represent an integer
    :param unique_values: list of numeric or string values
    :return: True if any value is an integer
    '''
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or re.search('[0-9]+',str(value)):
            status = True
            break
    return status

def contains_alpha(unique_values):
    '''
    Accepts a list of values and determines if any contain [aA -zZ]
    :param unique_values:
    :return:
    '''
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or isinstance(value, float):
            continue
        if re.search('[a-zA-Z]+',str(value)):
            status = True
            break
    return status

def convert_allele_codes(unique_values,method):
    '''
    Accepts a list of values and either casts them to an integer for a valid allele or 0 for missing
    :param unique_values: list of numeric or string values
    :param method: [int, hash] to toggle behaviour
    :return: dict of original: integer allele code
    '''
    converted_values = {}
    counter = 1
    for idx,value in enumerate(unique_values):
        if method == 'int':
            if isinstance(value,int) or str(value).isnumeric():
                converted_values[unique_values[idx]] = int(value)
            else:
                converted_values[unique_values[idx]] = 0
        elif method == 'hash':
            if value == '0':
                converted_values[unique_values[idx]] = 0
            else:
                converted_values[unique_values[idx]] = counter
                counter+=1
        else:
            if re.search('[a-zA-Z]+',str(value)) or re.search('[^0-9a-zA-Z]+',str(value)):
                value = '0'
            converted_values[unique_values[idx]] = int(value)
    return converted_values


def update_column_map(c1,c2):
    '''
    Adds a reference of each column from c2 into c1
    :param c1: dict
    :param c2: dict
    :return: dict
    '''
    for k in c2:
        if not k in c1:
            c1[k] = c2[k]

def is_all_columns_int(column_dtypes):
    '''

    :param column_dtypes: List of Pandas column dtypes
    :return: True if all columns are of an integer type
    '''
    count_non_int = 0
    for col in column_dtypes:
        if col in VALID_INT_TYPES:
            continue
        count_non_int+=1
    if count_non_int > 0:
        return False
    return True

def count_missing_data(df):
    '''
    Scans the df for missing data
    :param df: pandas dataframe
    :return: dict of column: and count of missing values per column
    '''
    counts = {}
    columns = df.columns.values.tolist()
    for c in columns:
        counts[c] = 0
        v = df[c].value_counts()
        if 0 in v:
            counts[c] = v[0]
    return counts

@jit(nopython=True)
def identify_cols_to_remove(column_counts,threshold):
    '''
    Identifies columns with missing data which exceeds the defined threshold
    :param column_counts: dict count of missing data per column
    :param threshold: float or in of the value to filter
    :return: list of columns which fail the test
    '''
    cols_to_remove = []
    for c in column_counts:
        if column_counts[c] > threshold:
            cols_to_remove.append(c)
    return cols_to_remove


def filter_columns(df,columns_to_remove):
    '''
    Removes columns from a data frame
    :param df:
    :param columns_to_remove:
    :return: data frame with columns removed
    '''
    return df.drop(columns_to_remove, axis=1)


def process_profile(profile_path,format="text",column_mapping={}):
    '''
    Reads in a file in (text, parquet) formats and applies processing to standardize the data and prepare it to be used
    for distance calculations
    :param profile_path: path to file
    :param format: format of the file [text, parquet]
    :param column_mapping: Previous allele code mapping to apply to the current file
    :return: (dict, pd)
    '''

    if format=='text':
        df = pd.read_csv(profile_path,header=0,sep="\t",index_col=0,low_memory=False)
    elif format=='parquet':
        df = pd.read_parquet(
            profile_path,
            engine='auto',
            columns=None,
            storage_options=None,
        )

    columns = df.columns.values.tolist()
    column_dtypes = df.dtypes.tolist()
    is_correct_format = is_all_columns_int(column_dtypes)

    #If all columns are already integers then skip the extra processing steps
    if is_correct_format:
        return (column_mapping, df)

    df = df.fillna('0')
    df = df.replace('?', '0', regex=False)
    df = df.replace(' ', '0', regex=False)
    df = df.replace('-', '0', regex=False)
    df = df.replace('', '0', regex=False)

    for column in columns:
        unique_col_values = sorted(df[column].unique().tolist())
        method = guess_format(List(unique_col_values))
        if not column in column_mapping:
            column_mapping[column] = convert_allele_codes(unique_col_values, method)
        else:
            update_column_map(column_mapping[column], convert_allele_codes(unique_col_values, method))

        df[column] = df[column].map(column_mapping[column])
    return (column_mapping, df)


def convert_profiles(df):
    '''
    Convert the pd data frame into a numpy array
    :param df: pandas df of integer allele profiles
    :return: (list,list) labels, numpy.arrays
    '''
    labels = df.index.tolist()
    profiles = []
    for index,row in df.iterrows():
        profiles.append(np.array(row.values.tolist()))
    return labels, profiles

def write_profiles(df,out_file,format):
    '''
    Write df to parquet
    :param df:
    :param out_file:
    :param format:
    :return:
    '''
    if format == 'parquet':
        df.to_parquet(out_file,compression='gzip')
    else:
        df.to_csv(out_file,sep="\t",header=True)

@jit(nopython=True)
def count_missing(p):
    '''
    Counts number of 0 occurances in the data
    :param p:
    :return:
    '''
    count = 0
    for idx,value in enumerate(p):
        if value ==0:
            count+=1

    return count


@jit(nopython=True)
def get_distance_raw(p1, p2):
    '''
    Calculates hamming distance with missing data skipped
    :param p1: numpy profile
    :param p2: numpy profile
    :return: int differences
    '''
    count = 0
    for v1,v2 in zip(p1,p2):
        if v1 == 0 or v2 == 0:
            continue
        if v1 != v2:
            count+=1
    return count

@jit(nopython=True)
def get_distance_scaled(p1, p2):
    '''
    Calculates % distance with missing data skipped
    :param p1: numpy profile
    :param p2: numpy profile
    :return: float differences
    '''
    count_compared_sites = 0
    count_match = 0
    for v1,v2 in zip(p1,p2):
        if v1 == 0 or v2 == 0:
            continue
        count_compared_sites+=1
        if v1 == v2:
            count_match+=1

    if count_compared_sites:
        return 100.0 * (float(count_compared_sites) - float(count_match)) / float(count_compared_sites)
    else:
        return 100.0


@jit(nopython=True)
def get_distance_raw_missing(p1, p2):
    '''
    Calculates hamming distance with missing data counted as differences
    :param p1: numpy profile
    :param p2: numpy profile
    :return: int differences
    '''
    count = 0
    for v1, v2 in zip(p1, p2):
        if v1 != v2:
            count += 1
    return count


@jit(nopython=True)
def get_distance_scaled_missing(p1, p2):
    '''
    Calculates % distance with missing data counted as differences
    :param p1: numpy profile
    :param p2: numpy profile
    :return: float differences
    '''
    count_compared_sites = 0
    count_match = 0
    for v1, v2 in zip(p1, p2):
        count_compared_sites += 1
        if v1 == v2:
            count_match += 1

    if count_compared_sites:
        return 100.0 * (float(count_compared_sites) - float(count_match)) / float(count_compared_sites)
    else:
        return 100.0

def calc_batch_size(num_records,num_columns,byte_value_size,max_mem=None):
    '''
    Calculation for the number of simultaneous calcualtions of distance to occur concurrently
    :param num_records: number of records
    :param num_columns: number of columns in a profile
    :param byte_value_size: integer of memory required
    :param max_mem: value to specify the maximum memory to use in the program
    :return: int of the number of records to process concurently
    '''
    if max_mem == None:
        max_mem = psutil.virtual_memory()
        avail = max_mem.available
    else:
        avail = max_mem
    p = (byte_value_size * num_columns) + 56
    profile_mem = p * num_records * 4
    dist_mem = ((byte_value_size * num_records) + 56) ** 2
    estimated_mem_needed = profile_mem + dist_mem
    if estimated_mem_needed < avail:
        return num_records

    avail -= profile_mem
    num_batches = int(avail / ((byte_value_size * num_records) + 56))
    batch_size = int(num_records / num_batches)

    if batch_size <= 0:
        batch_size = 100
    elif batch_size > num_records:
        batch_size = num_records

    return batch_size

@jit(nopython=True)
def validate_file(f):
    '''
    Helper function to determine a file exists and is not empty
    :param f: string path to file
    :return:  True on success
    '''
    if not os.path.isfile(f):
        return False

    if os.path.getsize(f) < MIN_FILE_SIZE:
        return False

    return True

def compare_headers(file1,file2):
    '''
    Compares two files to see if their headers are the same
    :param file1: string path to file 1
    :param file2: string path to file 2
    :return: True on success
    '''
    h1 = []
    h2 = []
    with open(file1,'r') as f1:
        h1 = next(f1).rstrip().split("\t")
        with open(file2, 'r') as f2:
            h2 = next(f2).rstrip().split("\t")
    if len(h1) > 0 and len(h2) > 0 and len(h1) == len(h2):
        ovl = set(h1) & set(h2)
        if len(ovl) == len(h1):
            return True
    return False

def guess_profile_format(f):
    '''
    Helper function to determine what file type a file is
    :param f: string path to file
    :return: string of the format
    '''
    ext = FILE_FORMATS
    ftype = ''

    for format in ext:
        for e in ext[format]:
            if f.endswith(e):
                ftype = format
                break
        if ftype != '':
            break

    return ftype


def get_file_length(f):
    '''
    Counts the number of lines in a file
    :param f: string path to file
    :return: int
    '''
    return int(os.popen(f'wc -l {f}').read().split()[0])


def calc_distances_scaled(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):
    '''
    Calculates pairwise distances between ref and query with the distances scaled according to the number of alleles
    considered between two samples with missing data ignored
    :param query_profiles: list of integer numpy profiles
    :param query_labels:  list of data labels for query
    :param ref_profiles: list of integer numpy profiles
    :param ref_labels: list of data labels for ref
    :param parquet_file: string path
    :param batch_size: int number of records to process concurrently
    :return: None
    '''

    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            if i != k:
                d.append(get_distance_scaled(query_profiles[i], ref_profiles[k]))
            else:
                d.append(0)
        dists.append(d)
        count += 1

        if count == batch_size:
            sys.stderr.write(f"{i} batch\n")
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                sys.stderr.write(f"write\n")
                fp.write(parquet_file, df, compression='GZIP')
            else:
                sys.stderr.write(f"append\n")
                fp.write(parquet_file, df, append=True, compression='GZIP')
            del(df)
            dists = []
            count = 0
            gc.collect()

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')

def calc_distances_scaled_missing(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):
    '''
    Calculates pairwise distances between ref and query with the distances scaled according to the number of alleles
    considered between two samples with missing data counted as a difference
    :param query_profiles: list of integer numpy profiles
    :param query_labels:  list of data labels for query
    :param ref_profiles: list of integer numpy profiles
    :param ref_labels: list of data labels for ref
    :param parquet_file: string path
    :param batch_size: int number of records to process concurrently
    :return: None
    '''
    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            if i != k:
                d.append(get_distance_scaled_missing(query_profiles[i], ref_profiles[k]))
            else:
                d.append(0)
        dists.append(d)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression='GZIP')
            else:
                fp.write(parquet_file, df, append=True, compression='GZIP')
            dists = []
            count = 0

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')

def calc_distances_hamming(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):
    '''
    Calculates pairwise hamming distances between ref and query with  with missing data ignored
    :param query_profiles: list of integer numpy profiles
    :param query_labels:  list of data labels for query
    :param ref_profiles: list of integer numpy profiles
    :param ref_labels: list of data labels for ref
    :param parquet_file: string path
    :param batch_size: int number of records to process concurrently
    :return: None
    '''
    count = 0
    columns = ["dists"] + ref_labels
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            if i != k:
                d.append(get_distance_raw(query_profiles[i], ref_profiles[k]))
            else:
                d.append(0)
        dists.append(d)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression='GZIP')
            else:
                fp.write(parquet_file, df, append=True, compression='GZIP')
            dists = []
            count = 0

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')


def calc_distances_hamming_missing(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):
    '''
    Calculates pairwise hamming distances between ref and query with  with missing data counted as differences
    :param query_profiles: list of integer numpy profiles
    :param query_labels:  list of data labels for query
    :param ref_profiles: list of integer numpy profiles
    :param ref_labels: list of data labels for ref
    :param parquet_file: string path
    :param batch_size: int number of records to process concurrently
    :return: None
    '''
    count = 0
    columns = ["dists"] + ref_labels
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            if i != k:
                d.append(get_distance_raw_missing(query_profiles[i], ref_profiles[k]))
            else:
                d.append(0)
        dists.append(d)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression='GZIP')
            else:
                fp.write(parquet_file, df, append=True, compression='GZIP')
            dists = []
            count = 0

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')




def is_file_ok(f):
    '''
    Helper function to determine if a profile file exists, has a header and >= 1 row of data
    :param f:
    :return: True on success
    '''
    status = True
    if not os.path.isfile(f):
        status = False
    elif get_file_length(f) < 2:
        status = False
    elif os.path.getsize(f) < MIN_FILE_SIZE:
        status = False

    return status

@jit(nopython=True)
def filter_dists(labels,distances,threshold):
    '''
    Filter a list of distances to ones less than or equal to the threshold
    :param labels: list of data labels
    :param distances: list of float/int distances
    :param threshold: float/int
    :return: dict {label: distance} of distances <= threshold
    '''
    results = {}
    for id, value in zip(labels,distances):
        if value <= threshold:
            results[id] = value
    return results


def format_pairwise_dist(df,threshold=-1):
    '''
    Accepts a pd dataframe in matrix format and converts it into query_id,ref_id, distance format df for
    distances <= thresh
    :param df: input df
    :param threshold: float/int
    :return: pd df query_id,ref_id with distances filtered by the threshold
    '''
    dists = {}
    columns = df.columns.values.tolist()[1:]
    for index,row in df.iterrows():
        dists[row[0]] = {row[0]:0}
        if threshold != -1:
            dists[row[0]] = filter_dists(List(columns), List(row[1:]), threshold)
        else:
            dists[row[0]] = dict(zip(columns, row[1:]))
        dists[row[0]] = {k: v for k, v in sorted(dists[row[0]].items(), key=lambda item: item[1])}

    results = {
        'query_id':[],
        'ref_id':[],
        'dist':[]
    }

    for qid in dists:
        results['query_id'] += [qid] * len(dists[qid])
        results['ref_id'] += list(dists[qid].keys())
        results['dist'] += list(dists[qid].values())


    return pd.DataFrame(results)


def write_dist_results(mat,outfile,outtype,outfmt,batch_size=1,threshold=-1):
    '''
    Writes the distance results into a matrix or pairwise distance result file
    :param mat: string path to distance matrix file
    :param outfile: string path for results to be written to
    :param outtype: string (pairwise, matrix)
    :param outfmt: string (text, parquet) file type of output file
    :param batch_size: int number of records to process concurrently
    :param threshold: int/float to filter results (only applies to pairwise)
    :return: None
    '''

    #If the desired output is a matrix in parquet format simply rename the mat file
    if outtype == 'matrix' and outfmt == 'parquet':
        os.rename(mat,outfile)
        return
    init_file = True
    parquet_file = pq.ParquetFile(mat)
    for batch in parquet_file.iter_batches(batch_size):
        batch_df = batch.to_pandas()
        del (batch)
        if outtype == 'pairwise':
            batch_df = format_pairwise_dist(batch_df, threshold=threshold)
        if init_file:
            init_file = False
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(outfile,index = False, header = True, sep="\t")
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(outfile, index=False, header=True, sep="\t")
            else:
                if not os.path.isfile(outfile):
                    fp.write(outfile, batch_df, compression='GZIP')
        else:
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(outfile, mode ='a', index = False, header = False, sep="\t")
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(outfile, mode ='a', index = False, header = False, sep="\t")
            else:
                fp.write(parquet_file, batch_df, append=True, compression='GZIP')


def get_missing_loci_counts(profiles,labels,count_loci):
    '''
    Counts the % missing data by sample
    :param profiles: list numpy arrays
    :param labels: list labels
    :param count_loci: total number of loci in profile
    :return: dict {label: % missing}
    '''
    n = len(labels)
    counts = {}
    for i in range(0, n):
        counts[labels[i]] = count_missing(profiles[i])/count_loci
    return counts

def flag_samples(missing_counts,threshold):
    '''
    Produces a list of labels which exceed the threshold for missing data
    :param missing_counts: dict of {label: % missing}
    :param threshold: int/float
    :return: list of labels
    '''
    r = []
    for sample_id in missing_counts:
        if missing_counts[sample_id] > threshold:
            r.append(sample_id)
    return sorted(r)

def filter_samples(labels,profiles,labels_to_remove):
    '''
    Removes samples from profile which are flagged to be removed
    :param labels: list of labels
    :param profiles: list of numpy profiles
    :param labels_to_remove:  List of labels to remove
    :return: (list, list) labels, list numpy profiles
    '''
    l = []
    p = []
    for idx,label in enumerate(labels):
        if label in labels_to_remove:
            continue
        l.append(label)
        p.append(profiles[idx])
    return l, p

