from profile_dists.version import __version__
MIN_FILE_SIZE = 32


EXTENSIONS = {'text': ['txt','tsv','mat','text'],
    'hd5': ['hd','h5','hdf5'],
    'parquet': ['parq','parquet','pq']}


FILE_FORMATS = ['tsv','parquet','json']

VALID_INT_TYPES = ['int64','int32','int16','int8']


OUTPUT_FILES = [
    'run.json',
    'allele_map.json',
    'results.{format}',
]


RUN_DATA = {
    'profile_dists': f'version: {__version__}',
    'analysis_start_time':'',
    'analysis_end_time':'',
    'parameters':{},
    'query_profile_info':{
        'num_samples':0,
        'num_samples_pass':0,
        'failed_samples':[],
        'parsed_file_path':'',
    },
    'ref_profile_info':{
        'num_samples': 0,
        'num_samples_pass': 0,
        'failed_samples': [],
        'parsed_file_path':'',
    },
    'loci_removed': [],
    'result_file':''
}