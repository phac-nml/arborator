import pandas as pd
import os

class read_data:
    status = True
    messages = []
    def __init__(self,input_file,MIN_FILE_SIZE=32):
        self.MIN_FILE_SIZE = MIN_FILE_SIZE
        self.input_file = input_file
        status = self.is_file_ok(self.input_file)
        if status:
            self.df = self.process_profile(input_file)
        else:
            self.df = pd.DataFrame()
            self.messages.append(f"Error unable to process {input_file}: is_file:{os.path.isfile(input_file)}")
        self.status = status

    def is_file_ok(self,f):
        '''
        Helper function to determine MIN_FILE_SIZEif a profile file exists, has a header and >= 1 row of data
        :param f:
        :return: True on success
        '''
        status = True
        if not os.path.isfile(f):
            status = False
        elif os.path.getsize(f) < self.MIN_FILE_SIZE:
            status = False
        elif self.is_file_reasonable(f) < 2:
            status = False


        return status

    def is_file_reasonable(self,f):
        '''
        Read up to first two lines of a file and return true if >= 2
        :param f: string path to file
        :return: boolean
        '''
        with open("r", encoding="utf-8", errors="replace") as f:
            # Only read up to 2 lines for efficiency
            line_count = 0
            for _ in f:
                line_count += 1
                if line_count >= 2:
                    return True
        return False

    def process_profile(self,file_path, format="text"):
        '''
        Reads in a file in (text, parquet) formats and produces a df
        :param profile_path: path to file
        :param format: format of the file [text, parquet]
        :return:  pd
        '''

        if format == 'text':
            df = pd.read_csv(file_path, header=0, sep="\t", low_memory=False)
        elif format == 'parquet':
            df = pd.read_parquet(
                file_path,
                engine='auto',
                columns=None,
                storage_options=None,
            )
        else:
            df = pd.DataFrame()

        return df