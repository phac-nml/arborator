import sys
from statistics import mean, median
import pandas as pd


class summarizer:
    valid_types = [
        'categorical',
        'min_max',
        'desc_stats',
        'none'
    ]

    def __init__(self,header,dictionary,field_data_types,field_name_key,field_name_value):
        self.header = sorted(header + self.get_fields(dictionary,field_data_types,field_name_key,field_name_value))

        self.data = self.populate_records(dictionary,field_data_types,field_name_key,field_name_value)


    def get_fields(self,dictionary,field_data_types,field_name_key,field_name_value):
        output_fields = []
        for group_id in dictionary:
            m = dictionary[group_id][field_name_key]
            for col in m:
                col_dtype = 'categorical'
                c = col.split(' ')
                if 'age' in c:
                    col_dtype = 'desc_stats'
                if 'date' in c:
                    col_dtype = 'min_max'
                c = col.split('_')
                if 'age' in c:
                    col_dtype = 'desc_stats'
                if 'date' in c:
                    col_dtype = 'min_max'
                if col in field_data_types:
                    col_dtype = field_data_types[col]["data_type"]
                if col_dtype not in self.valid_types:
                    col_dtype = 'categorical'

                if col_dtype == 'categorical':
                    values = list(m[col][field_name_value].keys())
                    for value in values:
                        k = f"count_{col}_{value}"
                        output_fields.append(k)

                elif col_dtype == 'min_max':
                    output_fields += [f'{col}_min_value',f'{col}_max_value']
                elif col_dtype == 'desc_stats':
                    output_fields += [f'{col}_min_value', f'{col}_mean_value', f'{col}_median_value',f'{col}_max_value']
                elif col_dtype == 'None':
                    output_fields.append(col)
        output_fields = sorted(list(set(output_fields)))
        return output_fields

    def create_record(self,header,field_data_types):
        record = {}
        for f in header:
            record[f] = ''
            if ('count_' in f or '_value') and not '_date' in f:
                record[f] = 0
            if f in field_data_types and "default" in field_data_types[f]:
                record[f] = field_data_types[f]["default"]
        return record


    def calc_desc_stats_numerical(self,values):
        s = {
            'min':'nan',
            'median':'nan',
            'mean':'nan',
            'max':'nan',
        }
        if len(values) > 0:
            s = {
                'min': min(values),
                'median': mean(values),
                'mean': median(values),
                'max':max(values),
            }
        return s

    def calc_desc_stats_dates(self,values):
        values = sorted(self.convert_date(values))
        s = {
            'min': 'nan',
            'median': 'nan',
            'mean': 'nan',
            'max': 'nan',
        }
        if len(values) > 0:
            s = {
                'min': values[0],
                'median': '',
                'mean': '',
                'max':values[-1],
            }
        return s

    def is_numbers(self,values):
        df = pd.DataFrame(data={'column1':values})
        try:
            df.astype(float)
        except:
            return False
        return True


    def is_date(self,values):
        d = self.convert_date(values)
        if len(d) > 0:
            return True
        return False

    def convert_date(self,values):
        return pd.to_datetime(values, errors='coerce',format="%Y-%m-%d").dropna()



    def populate_records(self,dictionary,field_data_types,field_name_key,field_name_value):
        records = {}
        for group_id in dictionary:
            records[group_id] = self.create_record(self.header,field_data_types)
            data = dictionary[group_id][field_name_key]
            for col in data:
                col_dtype = 'categorical'
                records[group_id][col] = ",".join(sorted([str(x) for x in list(data[col][field_name_value].keys())]))
                c = col.split(' ')
                if 'age' in c:
                    col_dtype = 'desc_stats'
                if 'date' in c:
                    col_dtype = 'min_max'
                c = col.split('_')
                if 'age' in c:
                    col_dtype = 'desc_stats'
                if 'date' in c:
                    col_dtype = 'min_max'
                
                if col in field_data_types:
                    if 'data_type' in field_data_types[col]:
                        col_dtype = field_data_types[col]['data_type']

                if col_dtype == 'categorical':
                    names = list(data[col][field_name_value].keys())
                    for n in names:
                        k = f"count_{col}_{n}"
                        records[group_id][k] = data[col][field_name_value][n]
                else:
                    if col_dtype == 'none':
                        records[group_id][col] = ",".join([str(x) for x in list(data[col][field_name_value].keys())])
                    else:
                        if col_dtype == 'desc_stats' or col_dtype == 'min_max':
                            values = list(data[col][field_name_value].keys())
                            if self.is_numbers(values):
                                values = []
                                for k in data[col][field_name_value]:

                                    values += [float(k)] * data[col][field_name_value][k]
                                dstats = self.calc_desc_stats_numerical(values)
                            else:
                                if self.is_date(list(data[col][field_name_value].keys())):
                                    dstats = self.calc_desc_stats_dates(values)
                                else:
                                    values = list(data[col][field_name_value].values())
                                    if self.is_date(values):
                                        dstats = self.calc_desc_stats_dates(values)


                            records[group_id][f'{col}_min_value'] = dstats['min']
                            records[group_id][f'{col}_max_value'] = dstats['max']

                            if col_dtype == 'desc_stats':
                                records[group_id][f'{col}_mean_value'] = dstats['mean']
                                records[group_id][f'{col}_median_value'] = dstats['median']
                        else:
                            records[group_id][col] = data[col][field_name_value]

        return records

    def get_data(self):
        return self.data

class merger:

    def __init__(self,left_dict,right_dict):
        self.fields = list(set(self.get_fields(left_dict) + self.get_fields(right_dict)))
        self.data = self.populate_records(self.data,left_dict)
        self.data = self.populate_records(self.data, right_dict)

    def get_data(self):
        return self.data

    def get_fields(self, dictionary):
        output_fields = []
        for id in dictionary:
            output_fields += list(dictionary[id].keys())
        return output_fields

    def create_record(self,header):
        record = {}
        for f in header:
            record[f] = ''


        return record

    def populate_records(self,merge_dict,dictionary):
        for id in dictionary:
            if not id in merge_dict:
                merge_dict[id] = self.create_record(self.fields)
            data = dictionary[id]
            for col in data:
                merge_dict[id][col] = data[col]


        return merge_dict














