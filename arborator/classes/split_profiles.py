from arborator.classes.read_data import read_data

class split_profiles:
    partition_file = None
    partitions = None
    groups = None
    subsets = {}
    df = None

    def __init__(self,df,partition_file,id_col,partition_col):
        self.df = df
        self.partition_file = partition_file
        self.id_col = id_col
        self.partition_col = partition_col
        self.parse_partition_file(id_col,partition_col)
        self.create_groups()
        self.subsets = self.subset_df(df,id_col)
        return

    def parse_partition_file(self,id_col,partition_col):
        partition = read_data(self.partition_file)
        df = partition.df.dropna(subset=[partition_col])
        self.partitions = dict(zip(df[id_col].values.tolist(), df[partition_col].values.tolist()))


    def create_groups(self):
        self.groups = {}
        for sample_id in self.partitions:
            group_id = self.partitions[sample_id]
            if not group_id in self.groups:
                self.groups[group_id] = []
            self.groups[group_id].append(sample_id)

    def subset_df(self,df,id_col):
        subsets = {}
        for group_id in self.groups:
            subsets[group_id] = df[df[id_col].isin(self.groups[group_id])]
        return subsets


