from Bio import SeqIO
from os import path
from general import *

BS_PATH = os.path.join(DATA_PATH, 'BS168.gb')
PART_A_CSV_PATH = os.path.join(DATA_PATH, "part_a.csv")
EXCEPTION_PATH = os.path.join(DATA_PATH, "gene_exceptions.csv")


class GenBank:
    def __init__(self, file_path, id_header="locus_tag"):
        df, seq, g_names = create_dataframe_from_gb(file_path, id_header)
        self.df = df
        self.seq = seq
        self.g_names = g_names
        self.cds = None
        self.other = None  # not cds

    # Q1
    def show_elements_dict(self):
        gene_dict = {'gene': len(self.g_names)}
        types_dict = self.df.type.value_counts().to_dict()
        elements_dict = dict(gene_dict)
        elements_dict.update(types_dict)
        print(f"Types dictionary: {elements_dict}")

    # Q2.a
    def add_len_col(self):
        self.df['len'] = self.df.apply(lambda row: np.abs(row['start'] - row['end']), axis=1)

    # Q2.b
    def set_cds_and_not_cds(self):
        cds, other = cds_other(self.df)
        self.cds = cds
        self.other = other

    # Q2.c,d
    def plot_len_statistics(self):
        total_len = np.asarray(self.df['len'])
        cds_len = np.asarray(self.cds['len'])
        other_len = np.asarray(self.other['len'])

        show_stat(cds_len, 'Proteins lengths stats')
        show_stat(other_len, 'Non proteins lengths stats')

        x_max = max([np.max(cds_len), np.max(other_len)])
        y_max = len(self.g_names)
        x_label = 'length'
        y_label = 'count'
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 3, 1)
        plot_hist("All Genes", total_len, x_label, y_label, x_max, y_max)

        plt.subplot(1, 3, 2)
        plot_hist("Proteins Genes", cds_len, x_label, y_label, x_max, y_max)

        plt.subplot(1, 3, 3)
        plot_hist("Non Proteins Genes", other_len, x_label, y_label, x_max, y_max)

        plt.suptitle("Lengths Statistics")
        plt.tight_layout()
        plt.show()

    # Q3.a
    def show_AT_percentage(self):
        at_percent = calc_percentage_in_genes(self.seq, ['A', 'T'])
        print(f"Genome AT percent: {at_percent:.2f}%")

    # Q3.b
    def add_sub_sequence_col(self):
        self.df['sub sequence'] = self.df.apply(lambda row: self.seq[row['start']:row['end']], axis=1)

    def add_at_percent_col(self):
        self.df['AT percent'] = self.df.apply(lambda row: calc_percentage_in_genes(row['sub sequence'], ['A', 'T']),
                                              axis=1)
        self.set_cds_and_not_cds()

    # Q3.c
    def show_cds_AT_average_percent(self):
        gb.add_sub_sequence_col()
        gb.add_at_percent_col()
        cds_at_percent = self.cds['AT percent']
        at_average = cds_at_percent.mean()
        print(f'Proteins average AT percent: {at_average:.2f}%')

    # Q3.d
    def plot_cds_AT_percent_stats(self):
        cds_at_percent = self.cds['AT percent']
        at_max = max(cds_at_percent)
        x_max = at_max + 0.1 * at_max
        y_max = 1500
        plot_hist('Proteins AT Percent Histogram', cds_at_percent, 'AT percent', 'count', x_max, y_max)
        plt.show()

    # Q3.e
    def show_extreme_AT_percents_genes(self, n=5):
        high = self.df.nlargest(n, 'AT percent')
        low = self.df.nsmallest(n, 'AT percent')
        # high.to_csv(os.path.join(DATA_PATH, 'high_AT_percent.csv'))
        # low.to_csv(os.path.join(DATA_PATH, 'low_AT_percent.csv'))
        print('\033[93m' + "\nExtreme AT percents genes:" + '\033[0m')
        print(f"Top {n} AT percents genes details:\n{high.to_string()}")
        print(f"\nBottom {n} AT percents genes details:\n{low.to_string()}")

    # Q4
    def show_cell_wall_stat(self):
        # a
        cell_wall = self.df[self.df['cell wall'] == 'Yes']
        print(f'The number of genes with the word "cell wall": {len(cell_wall)}')

        # b
        cell_wall_len = np.asarray(cell_wall['len'])
        show_stat(cell_wall_len, 'cell wall stats(length)')
        len_max = np.max(cell_wall_len)
        x_max = len_max + 0.1 * len_max
        plot_hist('cell wall length', cell_wall_len, 'length', 'count', x_max, 25)
        plt.show()

        # c
        cell_wall_at_percent = np.asarray(cell_wall['AT percent'])
        show_stat(cell_wall_at_percent, 'cell wall stats(%AT)')
        at_max = np.max(cell_wall_at_percent)
        x_max = at_max + 0.1 * at_max
        plot_hist('cell wall %AT', cell_wall_at_percent, '%AT', 'count', x_max, 10)
        plt.show()

    # Q5
    def report_conflicts(self, results_path):
        self.df['check'] = self.df.apply(
            lambda row: check_translation(row['type'], row['translation'], row['sub sequence'], row['table'],
                                          row['strand'], row['codon_start']), axis=1)
        check_df = self.df.dropna(subset=['check'])
        err_df = check_df[check_df['check'] != 'OK']

        print(f'Genes with conflict in translations:\n{err_df.to_string()}')
        err_df.to_csv(results_path)


if __name__ == "__main__":
    gb = GenBank(BS_PATH)

    print('\033[92m' + '--- PartA | Q1 ---' + '\033[0m')
    gb.show_elements_dict()

    print('\n\033[92m' + '--- PartA | Q2 ---' + '\033[0m')
    gb.add_len_col()
    gb.set_cds_and_not_cds()
    gb.plot_len_statistics()

    print('\n\033[92m' + '--- PartA | Q3 ---' + '\033[0m')
    gb.show_AT_percentage()
    gb.show_cds_AT_average_percent()
    gb.plot_cds_AT_percent_stats()
    gb.show_extreme_AT_percents_genes()

    print('\n\033[92m' + '--- PartA | Q4 ---' + '\033[0m')
    gb.show_cell_wall_stat()

    print('\n\033[92m' + '--- PartA | Q5 ---' + '\033[0m')
    gb.report_conflicts(EXCEPTION_PATH)

    sorted_df = gb.df.sort_values(by='id')
    sorted_df.to_csv(PART_A_CSV_PATH, index=False)
