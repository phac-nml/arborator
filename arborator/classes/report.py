import sys

from scipy.stats import entropy

class report:

    def __init__(self,df,columns_to_skip=[],scale=True,missing_data='0'):
        self.loci = self.get_col_counts(df,columns_to_skip,scale=True,missing_data='0')
        pass

    def get_col_counts(self,df,columns_to_skip=[],scale=True,missing_data='0'):
        loci = {}
        columns = list(df.columns)
        for col in columns:
            if col in columns_to_skip:
                continue
            locus = {
                'num_values': 0,
                'num_missing': 0,
                'value_counts': {},
                'shannon_entropy': -1,
            }
            unique_values = df[col].value_counts(dropna=True)
            unique_values.index = unique_values.index.astype(str)
            unique_values = dict(unique_values)

            if missing_data in unique_values:
                locus['num_missing'] = unique_values[missing_data]
                del (unique_values[missing_data])
            locus['num_values'] = len(unique_values)
            locus['value_counts'] = unique_values
            value_list = list(unique_values.values())
            if len(value_list) > 1:
                locus['shannon_entropy'] = self.calc_shanon_entropy(value_list,scale)
            elif len(value_list) == 1:
                locus['shannon_entropy'] = 0
            loci[col] = locus
        return loci

    def calc_shanon_entropy(self,value_list,scale=True):

        total = sum(value_list)
        values = []
        for v in value_list:
            values.append(v / total)

        if scale:
            return entropy(values) / entropy([1] * len(values))
        else:
            return entropy(values)


    def write_data(self,outfile):
        with open(outfile,'w') as oh:
            oh.write("locus\tnum_values\tnum_missing\tshannon_entropy\tvalue_counts\n")
            for l in self.loci:
                row = [
                    l,
                    self.loci[l]['num_values'],
                    self.loci[l]['num_missing'],
                    self.loci[l]['shannon_entropy'],
                    self.loci[l]['value_counts']
                ]
                oh.write("{}\n".format("\t".join([str(x) for x in row])))


    def get_data(self):
        return self.loci


