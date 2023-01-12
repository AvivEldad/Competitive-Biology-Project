import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import general as ge

PART_A_CSV_PATH = os.path.join(ge.DATA_PATH, "part_a.csv")
UNIPROT_PATH = os.path.join(ge.DATA_PATH, "uniprot_BS168.csv")


class UniProt:

    def __init__(self, genbank_path, uniport_path):
        self.genbank_file = pd.read_csv(genbank_path)
        self.uniprot_file = pd.read_csv(uniport_path)
        self.trans_df = None

    def prepre_data(self, names):
        prepared = []
        nan_count = 0
        for name in names:
            if str(name) == 'nan':
                nan_count += 1
            if str(name) != 'nan':
                name = str(name).replace('_', '').replace(';', ' ').replace('/', ' ')
                multi = str(name).split()
                for i in multi:
                    prepared.append(i)
        return prepared, nan_count

    def plot_pie_charts(self, total, nan_values, unique_genes, name):
        if nan_values == 0:
            y = np.array([total, unique_genes])
            my_labels = ["total", "not in {}".format('Uniprot' if name == 'GeneBank' else 'GeneBank')]
            exp = [0.05, 0.3]
        else:
            y = np.array([total, nan_values, unique_genes])
            my_labels = ["total", "non-values", "not in {}".format('Uniprot' if name == 'GeneBank' else 'GeneBank')]
            exp = [0.05, 0.01, 0.3]

        plt.pie(y, labels=my_labels, explode=exp, shadow=True, autopct='%1.2f%%')
        plt.title("{} genes".format(name))
        plt.show()

    def compare_uni_and_gb(self):
        uniprot_names = np.asarray(self.uniprot_file['Gene Names (ordered locus)'])
        gb = self.genbank_file[self.genbank_file['type'] == 'CDS']
        gb_names = np.asarray(gb['id'])

        # prepre the data
        uni_to_compare, uni_nan = self.prepre_data(uniprot_names)
        gb_to_compare, gb_nan = self.prepre_data(gb_names)

        in_uni = list(set(uni_to_compare).difference(gb_to_compare))
        in_gb = list(set(gb_to_compare).difference(uni_to_compare))

        print('\033[93m' + "UniProt and not in GenBank:" + '\033[0m' + "\n{}".format(in_uni))
        print("missing: {}, nan values: {}, total: {}".format(len(in_uni), uni_nan, len(uniprot_names)))
        print('\033[93m' + "\nGenBank and not in UniProt:" + '\033[0m' + "\n{}".format(in_gb))
        print("missing: {}, nan values: {}, total: {}".format(len(in_gb), gb_nan, len(gb_names)))

        self.plot_pie_charts(len(uniprot_names), uni_nan, len(in_uni), "Uniprot")
        self.plot_pie_charts(len(gb_names), gb_nan, len(in_gb), "GeneBank")

    def trans_operations(self):
        non_empty_df = self.uniprot_file.dropna(subset=['Transmembrane'])
        start_list, end_list, len_list, gene_list, seq_list = [], [], [], [], []
        for index, row in non_empty_df.iterrows():
            name = row['Gene Names (ordered locus)']
            seq = row['Sequence']
            trans = row['Transmembrane'].replace('TRANSMEM ', '').replace(' ', '').split(';')
            trans = [t for t in trans if '\"' not in t]
            for i in trans:
                splitter = i.split('..')
                start, end = int(splitter[0]), int(splitter[1])
                start_list.append(start)
                end_list.append(end)
                len_list.append(np.abs(end - start))
                gene_list.append(name)
                seq_list.append(seq[start:end])

        self.trans_df = pd.DataFrame(zip(gene_list, start_list, end_list, len_list, seq_list),
                                     columns=['Id', 'Start', 'End', 'Length', 'Sequence'])
        total_len = np.asarray(self.trans_df['Length'])
        ge.show_stat(total_len, 'Transmembrane lengths stats')
        ge.plot_hist("Transmembrane lengths", total_len, 'length', 'count', np.max(total_len) + 5, 8000)
        plt.show()

    def hydro_operations(self):
        sequences = self.trans_df['Sequence']
        seq_percent = []
        for seq in sequences:
            hydro = sum([1 for c in seq if c in ge.HYDRO_AMINO])
            hydro_percent = (hydro / len(seq)) * 100
            seq_percent.append(hydro_percent)

        ge.show_stat(seq_percent, "Hydrophobic Amino Percent in Transmembrane sequences")
        ge.plot_hist("Hydrofobic amino acids", seq_percent, 'percent', 'count', 200, 4000)
        plt.show()

    def plot_subs_hist(self, at_percent_gb, intersection_seq_percent, complement_seq_percent):
        x_title = "AT percent"
        y_title = "numbers"
        x_max = 100
        y_max = 1500
        plt.figure(figsize=(15, 6))

        plt.subplot(1, 4, 1)
        ge.plot_hist("GenBank (A)", at_percent_gb, x_title, y_title, x_max, y_max)

        plt.subplot(1, 4, 2)
        ge.plot_hist("GenBank with Trans (B)",
                     intersection_seq_percent, x_title, y_title, x_max, y_max)

        plt.subplot(1, 4, 3)
        ge.plot_hist("GenBank without Trans (A/B)",
                     complement_seq_percent, x_title, y_title, x_max, y_max)

        plt.subplot(1, 4, 4)
        ge.plot_hist("GenBank with/out trans", complement_seq_percent, x_title, y_title, x_max, y_max,
                     color='black')
        ge.plot_hist("GenBank with/out trans", intersection_seq_percent, x_title, y_title, x_max, y_max)

        plt.suptitle("AT Percent")
        plt.tight_layout()
        plt.show()

    def AT_distribution(self):
        trans_names = np.asarray(self.trans_df['Id'])
        gb = self.genbank_file[self.genbank_file['type'] == 'CDS']
        gb_names = np.asarray(gb['id'])

        un_to_compare, uni_nan = gb_and_uni.prepre_data(trans_names)
        gb_to_compare, gb_nan = gb_and_uni.prepre_data(gb_names)

        at_percent_gb = np.asarray(gb['AT percent'])  # A team

        mask = np.array([(id in un_to_compare) for id in gb_to_compare])

        intersection = gb_names[mask]
        complement = gb_names[~mask]
        intersection_seq_percent = []  # B team
        for name in intersection:
            seq = gb.loc[gb['id'] == name, 'sub sequence'].iloc[0]
            intersection_seq_percent.append(ge.calc_percentage_in_genes(seq, 'AT'))
        ge.plot_hist("AT percent in intersection", intersection_seq_percent, 'percent', 'count', 100, 500)
        plt.show()

        complement_seq_percent = []
        for name in complement:
            seq = gb.loc[gb['id'] == name, 'sub sequence'].iloc[0]
            complement_seq_percent.append(ge.calc_percentage_in_genes(seq, 'AT'))

        ge.show_stat(at_percent_gb, "GenBank AT percent")
        ge.show_stat(intersection_seq_percent, "GenBank with Transmembrane intersection AT percents")
        ge.show_stat(complement_seq_percent, "GenBank without Transmembrane intersection AT percents")

        self.plot_subs_hist(at_percent_gb, intersection_seq_percent, complement_seq_percent)


if __name__ == "__main__":
    gb_and_uni = UniProt(PART_A_CSV_PATH, UNIPROT_PATH)

    print('\033[92m' + '--- PartB | Q1 ---' + '\033[0m')
    gb_and_uni.compare_uni_and_gb()

    print('\033[92m' + '\n--- PartB | Q2 ---' + '\033[0m')
    gb_and_uni.trans_operations()  # a
    gb_and_uni.hydro_operations()  # b

    print('\033[92m' + '\n--- PartB | Q3 ---' + '\033[0m')
    gb_and_uni.AT_distribution()
