from Bio import Align
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Align import substitution_matrices
import os
import general as ge
from tabulate import tabulate
import re


class VirusAnalyzer:

    def __init__(self, file1, file2):
        assert (os.path.exists(file1))  # making sure that the path is valid
        assert (os.path.exists(file2))  # making sure that the path is valid

        self.df_1, sequences, self.gene_names_1 = ge.create_dataframe_from_gb(file1, 'gene')
        self.cds_1, self.other_1 = ge.cds_other(self.df_1)
        self.file_id_1 = os.path.basename(file1)[:-3]
        self.df_2, sequences, self.gene_names_2 = ge.create_dataframe_from_gb(file2, 'gene')
        self.cds_2, self.other_2 = ge.cds_other(self.df_2)
        self.file_id_2 = os.path.basename(file2)[:-3]
        self.common_genes_list = None

    # Q1
    def calc_synonymous(self):
        # Create a dictionary of codons and their synonymous positions
        synonymous_positions = {}

        # Iterate over all codons
        for codon in ge.TRANS_TABLE_1:
            syn = 0
            nonsyn = 0
            # Count the number of positions at which the codon is synonymous with others
            trans = 'ATGC'
            # for each site
            for i in range(3):
                # for each transaction
                for t in trans:
                    trans_codon = codon[:i] + t + codon[i + 1:]

                    if trans_codon == codon or ge.TRANS_TABLE_1[trans_codon] == '_':
                        continue
                    if ge.TRANS_TABLE_1[trans_codon] == ge.TRANS_TABLE_1[codon]:
                        syn += 1
                    else:
                        nonsyn += 1
            synonymous_positions[codon] = float("%.3f" % (3 * syn / (syn + nonsyn)))

        print(synonymous_positions)

    # Q2_a
    def compare_gb_files(self):
        all_gene_len_1 = len(self.other_1) + len(self.gene_names_1)
        all_gene_len_2 = len(self.other_2) + len(self.gene_names_2)

        format_massage = '{}\n=============\nGenes Number: {}\nProteins Number: {}\n'
        print(format_massage.format(self.file_id_1, all_gene_len_1, len(self.cds_1)))
        print(format_massage.format(self.file_id_2, all_gene_len_2, len(self.cds_2)))

    # Q2_b
    def compare_genes(self):
        self.common_genes_list = list(set(self.df_1['id']) & set(self.df_2['id']))
        def1 = list(set(self.df_1['id']) - set(self.df_2['id']))
        def2 = list(set(self.df_1['id']) - set(self.df_2['id']))

        print('There are {} common values between the files.\nthe values are: {}\n'.format(
            len(self.common_genes_list), self.common_genes_list))

        formate_massage = 'the values that only in {} and not in {} are - {}'
        print(formate_massage.format(self.file_id_1, self.file_id_2, def1))
        print(formate_massage.format(self.file_id_2, self.file_id_1, def2))

    # Q3_c
    def common_genes_analysis(self):
        gene_lst = []
        table = []

        for gene in self.common_genes_list:
            if gene not in self.cds_1['id'].values:
                continue

            if self.cds_1.loc[self.cds_1['id'] == gene, 'check'].iloc[0] != 'OK' or \
                    self.cds_2.loc[self.cds_2['id'] == gene, 'check'].iloc[0] != 'OK':
                continue

            sub_seq1 = self.cds_1.loc[self.cds_1['id'] == gene, 'sub sequence'].iloc[0]
            sub_seq2 = self.cds_2.loc[self.cds_2['id'] == gene, 'sub sequence'].iloc[0]

            if self.validate_Seq(sub_seq1) is None or self.validate_Seq(sub_seq2) is None:
                continue

            gene_lst.append(gene)

            t1 = self.cds_1.loc[self.cds_1['id'] == gene, 'table'].iloc[0]
            t2 = self.cds_2.loc[self.cds_2['id'] == gene, 'table'].iloc[0]

            pro_id = self.cds_1.loc[self.cds_1['id'] == gene, 'protein id'].iloc[0]
            role = self.cds_1.loc[self.cds_1['id'] == gene, 'product'].iloc[0]

            aligned_seq1, aligned_seq2 = self.trans_and_align(sub_seq1, sub_seq2, t1, t2)

            dN, dS = cal_dn_ds(CodonSeq(aligned_seq1), CodonSeq(aligned_seq2))

            if dN == 0 and dS == 0:  # no mutations at all
                dN_dS = 0
                selection = "negative"

            elif dS == 0:  # there are no synonymous mutations
                dN_dS = None
                selection = "positive"

            else:
                dN_dS = dN / dS
                dN_dS = float("%0.3f" % dN_dS)  # round
                selection = ""
                if 0.95 <= dN_dS <= 1.05:
                    selection = "neutral"
                elif dN_dS > 1.05:
                    selection = "positive"
                else:
                    selection = "negative"

            # build table:
            table.append([gene, pro_id, role, dN, dS, dN_dS, selection])

            # choose only 5
            if len(gene_lst) == 5:
                break

        # define header names
        col_names = ["Gene name", "Id", "Product", "dN", "dS", "dN/dS", "Selection Type"]
        # display table
        print(tabulate(table, headers=col_names))
        print("\n", gene_lst)

    def trans_and_align(self, seq1, seq2, t1, t2):
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        messenger_rna1 = Seq(seq1).translate(t1)
        messenger_rna2 = Seq(seq2).translate(t2)

        alignments = aligner.align(messenger_rna1, messenger_rna2)
        return self.calc_by_codons(seq1, alignments[0][0][:-1]), self.calc_by_codons(seq2, alignments[0][1][:-1])

    def calc_by_codons(self, seq, trans_seq):
        ind = 0
        new_seq = ""
        for pro in range(0, len(trans_seq)):
            if trans_seq[pro] == '-':
                new_seq += '---'
            else:
                new_seq += seq[ind] + seq[ind + 1] + seq[ind + 2]
                ind += 3
        return new_seq

    def validate_Seq(self, seq):
        pattern = r'[^ACGTacgt]'  # pattern of non-nucleotide
        test_string = str(seq)
        result = re.search(pattern, test_string)
        if result is not None:
            return None
        return "true"


if __name__ == "__main__":
    covid_2021 = r'Data\MZ383039.1.gb'  # corona virus from 2021
    covid_2022 = r'Data\OQ065689.1.gb'  # corona virus from 2022

    va = VirusAnalyzer(covid_2021, covid_2022)

    print('\033[92m' + '--- PartC | Q1 ---' + '\033[0m')
    va.calc_synonymous()

    print('\n\033[92m' + '--- PartC | Q2 ---' + '\033[0m')
    print('\033[93m' + '2.a' + '\033[0m')
    va.compare_gb_files()

    print('\n\033[93m' + '2.b' + '\033[0m')
    va.compare_genes()

    print('\n\033[93m' + '2.c' + '\033[0m')
    va.common_genes_analysis()
