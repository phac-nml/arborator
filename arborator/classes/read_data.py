import pandas as pd
import os

class read_data:

    def __init__(self,input_file):
        self.input_file = input_file
        self.status = self.is_file_ok(self.input_file)
        self.messages = []

        if  self.status:
            self.df = self.process_profile(input_file)
        else:
            self.df = pd.DataFrame()
            self.messages.append(f"Error unable to process {input_file}: is_file:{os.path.isfile(input_file)}")

    def is_file_ok(self,f):
        '''
        Helper function to determine if a profile file exists, has a header and >= 1 row of data
        :param f:
        :return: True on success
        '''
        status = True
        if not os.path.isfile(f):
            status = False
        elif self.get_file_length(f) < 2:
            status = False

        return status

    def get_file_length(self,f):
        '''
        Counts the number of lines in a file
        :param f: string path to file
        :return: int
        '''
        return int(os.popen(f'wc -l {f}').read().split()[0])

    def process_profile(self,file_path, format="text"):
        '''
        Reads in a file in (text, parquet) formats and produces a df
        :param profile_path: path to file
        :param format: format of the file [text, parquet]
        :return:  pd
        '''

        if format == 'text':
            df = pd.read_csv(file_path, header=0, sep="\t", low_memory=False, dtype=str)
        elif format == 'parquet':
            df = pd.read_parquet(
                file_path,
                engine='auto',
                columns=None,
                storage_options=None,
            )

        return df