from arborator.classes.read_data import read_data
import re

class split_profiles:

    def __init__(self, df, partition_file, id_col, partition_col):
        self.df = df
        self.partition_file = partition_file
        self.id_col = id_col
        self.partition_col = partition_col

        self.parse_partition_file()
        self.create_groups()
        self.subsets = self.subset_df(df,id_col)

        return

    def parse_partition_file(self):
        partition = read_data(self.partition_file)
        data_frame = partition.df

        if self.partition_col not in data_frame:
            message = f'the partition column {self.partition_col} does not exist in the data'
            raise Exception(message)

        data_frame = data_frame.dropna(subset=[self.partition_col])
        data_frame = data_frame[[self.id_col, self.partition_col]].set_index(self.id_col)
        self.partitions = data_frame.to_dict()[self.partition_col]

    def create_groups(self):
        SPECIAL_REGEX = "[^A-Za-z0-9_\-.]" # For file names.
        REPLACEMENT_CHARACTER = "_"
        self.groups = {}
        self.used_filepaths = []
        self.group_file_mapping = {}
            # A file path-safe version of the ID for file writing.
            # Example: {"1": "1", "United Kingdom": "United_Kingdom"}

        for sample_id in self.partitions:
            group_id = str(self.partitions[sample_id])

            if not group_id in self.groups:
                self.groups[group_id] = []

                filepath = re.sub(SPECIAL_REGEX, REPLACEMENT_CHARACTER, group_id)

                # Collision check:
                # Replace different characters to the same character
                # can cause collisions.
                while filepath in self.used_filepaths:
                    filepath += "-1"

                self.used_filepaths.append(filepath)
                self.group_file_mapping[group_id] = filepath

            self.groups[group_id].append(sample_id)

    def subset_df(self,df,id_col):
        subsets = {}
        for group_id in self.groups:
            subsets[group_id] = df[df[id_col].isin(self.groups[group_id])]
        return subsets


