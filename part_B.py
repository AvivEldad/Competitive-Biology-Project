import pandas as pd
from Bio import SeqIO
from os import path
import numpy as np
import matplotlib.pyplot as plt
import general as ge


def prepre_data(names):
    prepred = []
    nan_count = 0
    for name in names:
        if (str(name) == 'nan'):
            nan_count += 1
        if (str(name) != 'nan'):
            name = str(name).replace('_', '').replace(';', ' ').replace('/', ' ')
            multi = str(name).split()
            for i in multi:
                prepred.append(i)
    return prepred, nan_count


def plot_pie_charts(total, nan_values, uniqe_genes, name):
    y = np.array([total, nan_values, uniqe_genes])
    mylabels = ["total", "non-values", "not in gb"]

    plt.pie(y, labels=mylabels, explode=[0.05, 0.05, 0.2], shadow=True, autopct='%1.3f%%')
    plt.title("{} genes".format(name))
    plt.show()


def Q1():
    uniprot_names = np.asarray(pd.read_csv('uniprot_BS168.csv')['Gene Names (ordered locus)'])
    df = pd.read_csv('part_a.csv')
    gb = df[df['type'] == 'CDS']
    gb_names = np.asarray(gb['id'])

    # prepre the data
    un_to_compre, uni_nan = prepre_data(uniprot_names)
    gb_to_compre, gb_nan = prepre_data(gb_names)

    in_uni = list(set(un_to_compre).difference(gb_to_compre))
    in_gb = list(set(gb_to_compre).difference(un_to_compre))

    print("uniprot and not in gb:\n{}".format(in_uni))
    print("\nmissing: {}, nan values: {}, total: {}".format(len(in_uni), uni_nan, len(uniprot_names)))
    print("gb and not in uniprot:\n{}".format(in_gb))
    print("\nmissing: {}, nan values: {}, total: {}".format(len(in_gb), gb_nan, len(gb_names)))

    plot_pie_charts(len(uniprot_names), uni_nan, len(in_uni), "Uniprot")
    plot_pie_charts(len(gb_names), gb_nan, len(in_gb), "GeneBank")


if __name__ == "__main__":

    # Q1
    # Q1()

    # Q2a

    non_empty_df = pd.read_csv('uniprot_BS168.csv').dropna(subset=['Transmembrane'])
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

    final_df = pd.DataFrame(zip(gene_list, start_list, end_list, len_list, seq_list),
                            columns=['Id', 'Start', 'End', 'Length', 'Sequence'])
    total_len = np.asarray(final_df['Length'])
    ge.show_stat(total_len, 'Transmembrane lengths stats')
    ge.plot_hist("Transmembrane lengths", total_len, 'length', 'count', np.max(total_len) + 5, 8000)
    plt.show()

    # Q2b

    HYDROPHOBIC_AMINO = ['A', 'F', 'L', 'I', 'V', 'M', 'P', 'W']
    sequences = final_df['Sequence']
    seq_prectange = []
    for seq in sequences:
        hydro = sum([1 for c in seq if c in HYDROPHOBIC_AMINO])
        hydro_percent = (hydro / len(seq)) * 100
        seq_prectange.append(hydro_percent)

    ge.show_stat(seq_prectange, "Hydrophobic Amino Percent in Transmembrane sequences")
