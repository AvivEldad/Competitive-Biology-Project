import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os import path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Data')

DNA = ['A', 'C', 'T', 'G']
HYDRO_AMINO = ['A', 'F', 'L', 'I', 'V', 'M', 'P', 'W']

TRANS_TABLE_1 = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


def check_translation(gene_type, org_trans, seq, trans_table, strand, codon_start):
    if gene_type != "CDS":
        return None

    codon_start = int(codon_start)
    org_seq = Seq(seq[(codon_start - 1):])
    if strand != 1:
        org_seq = org_seq.reverse_complement()

    try:
        trans = str(org_seq.translate(table=trans_table, cds=True))
        check = (trans == org_trans)
        if check:
            return "OK"

    except TranslationError as err:
        return str(err)


def parse_genbank(gb_path):
    assert (path.exists(gb_path))
    with open(gb_path, 'r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)

    return record_gb


def create_dataframe_from_gb(gb_path, id_header='locus_tag'):
    record_gb = parse_genbank(gb_path)
    sequence = record_gb.seq.upper()  # full genome
    genes_id, start, end, feat_types, strand, cell_wall = [], [], [], [], [], []
    tables, translations, codon_starts = [], [], []  # only for cds
    gene_names = []
    features = record_gb.features[1:]

    for feature in features:
        cw = 'No'
        if id_header in feature.qualifiers.keys():
            g_name = feature.qualifiers[id_header][0]
        else:
            g_name = None

        if feature.type == 'gene':
            gene_names.append(g_name)
            continue

        if feature.type == 'CDS':
            translation = feature.qualifiers['translation'][0]
            codon_start = feature.qualifiers['codon_start'][0]

            if 'transl_table' in feature.qualifiers.keys():
                table = feature.qualifiers['transl_table'][0]
            else:
                table = 1

            if 'product' in feature.qualifiers.keys():
                if 'cell wall' in feature.qualifiers['product'][0]:
                    cw = 'Yes'

        else:
            table = None
            translation = None
            codon_start = None

        genes_id.append(g_name)
        feat_types.append(feature.type)
        start.append(feature.location.start)
        end.append(feature.location.end)
        strand.append(feature.location.strand)
        cell_wall.append(cw)
        tables.append(table)
        translations.append(translation)
        codon_starts.append(codon_start)

    df = pd.DataFrame(zip(genes_id, start, end, strand, feat_types, tables, translations, codon_starts, cell_wall),
                      columns=['id', 'start', 'end', 'strand', 'type', 'table', 'translation', 'codon_start',
                               'cell wall'])
    df['sub sequence'] = df.apply(lambda row: sequence[row['start']:row['end']], axis=1)
    df['check'] = df.apply(
        lambda row: check_translation(row['type'], row['translation'], row['sub sequence'], row['table'], row['strand'],
                                      row['codon_start']), axis=1)

    return df, sequence, gene_names


def cds_other(df):
    cds = df[df['type'] == 'CDS']
    other = df[df['type'] != 'CDS']
    return cds, other


def show_stat(arr, title):
    if type(arr) != np.ndarray:
        arr = np.asarray(arr)
    minimum = np.min(arr)
    maximum = np.max(arr)
    avg = np.average(arr)
    std = np.std(arr)

    print('\033[93m' + f"{title}:\n" + '\033[0m'
                                       f"Average: {avg:.2f}\n"
                                       f"Minimum: {minimum:.2f}\n"
                                       f"Maximum: {maximum:.2f}\n"
                                       f"Standard deviation: {std:.2f}\n")


def plot_hist(title, arr, x_label, y_label, x_max, y_max, color='purple'):
    plt.title(title)
    plt.hist(arr, color=color)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim([0, x_max])
    plt.ylim([0, y_max])


def calc_percentage_in_genes(seq, letters):
    counts = {}
    for letter in letters:
        count = seq.count(letter)
        counts[letter] = count

    return (sum(counts.values()) / len(seq)) * 100
