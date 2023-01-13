from Bio import Align
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Align import substitution_matrices
import os
import general as ge
from tabulate import tabulate
import re

class Virus_Analyzer:

    def __init__(self, file1, file2):

        assert (os.path.exists(file1))  # making sure that the path is valid
        assert (os.path.exists(file2))  # making sure that the path is valid

        self.df_1, sequences, self.gene_names_1 = ge.create_dataframe_from_gb(
            file1, 'gene')
        self.cds_1, self.other_1 = ge.cds_other(self.df_1)
        self.file_id_1 = os.path.basename(file1)[:-3]

        self.df_2, sequences, self.gene_names_2 = ge.create_dataframe_from_gb(
            file2, 'gene')
        self.cds_2, self.other_2 = ge.cds_other(self.df_2)
        self.file_id_2 = os.path.basename(file2)[:-3]

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
                    trans_codon = codon[:i] + t + codon[i+1:]

                    if trans_codon == codon or ge.TRANS_TABLE_1[trans_codon] == '_':
                        continue
                    if ge.TRANS_TABLE_1[trans_codon] == ge.TRANS_TABLE_1[codon]:
                        syn += 1
                    else:
                        nonsyn += 1
            synonymous_positions[codon] = float("%.3f" % (3*syn/(syn+nonsyn)))

        print(synonymous_positions)

    # Q2_a
    def compare_gb_files(self):

        all_gene_len_1 = len(self.cds_1) + \
            len(self.other_1) + len(self.gene_names_1)
        all_gene_len_2 = len(self.cds_2) + \
            len(self.other_2) + len(self.gene_names_2)

        format_massage = '{}\n=============\nGenes Number: {}\nProteins Number: {}\n'
        print(format_massage.format(self.file_id_1,
              all_gene_len_1, len(self.cds_1)))
        print(format_massage.format(self.file_id_2,
              all_gene_len_2, len(self.cds_2)))

    # Q2_b

    def compare_genes(self):

        formate_massage = 'the values that only in {} and not in {} are - {}\n'

        self.common_genes_list = list(
            set(self.cds_1['id']) & set(self.cds_2['id']))
        def1 = list(set(self.cds_1['id'])-set(self.cds_2['id']))
        def2 = list(set(self.cds_2['id'])-set(self.cds_1['id']))

        print('There are {} common values between the files.\n the values are: {}\n'.format(
            len(self.common_genes_list), self.common_genes_list))

        formate_massage = 'the values that only in {} and not in {} are - {}\n'
        print(formate_massage.format(self.file_id_1, self.file_id_2, def1))
        print(formate_massage.format(self.file_id_2, self.file_id_1, def2))

    # Q3_c
    def common_genes_analysis(self):

        gene_lst = []
        table = []

        for gene in self.common_genes_list:

            sub_seq1 = self.cds_1.loc[self.cds_1['id']
                                      == gene, 'gene_seq'].iloc[0]
            sub_seq2 = self.cds_2.loc[self.cds_2['id']
                                      == gene, 'gene_seq'].iloc[0]

            if self.validate_Seq(sub_seq1) == None or self.validate_Seq(sub_seq2) == None:
                continue

            gene_lst.append(gene)

            t1 = self.cds_1.loc[self.cds_1['id']
                                == gene, 'table'].iloc[0]
            t2 = self.cds_2.loc[self.cds_2['id']
                                == gene, 'table'].iloc[0]

            pro_id = self.cds_1.loc[self.cds_1['id']
                                    == gene, 'protein id'].iloc[0]
            role = self.cds_1.loc[self.cds_1['id']
                                  == gene, 'product'].iloc[0]

            aligned_seq1, aligned_seq2 = self.trans_and_align(
                sub_seq1, sub_seq2, t1, t2)

            dN, dS = cal_dn_ds(CodonSeq(aligned_seq1), CodonSeq(aligned_seq2))

            if dN == 0 and dS == 0:  # no mutations at all
                dN_dS = 0
                selection = "negative"

            elif dS == 0:  # there are no synonimus mutations
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
        col_names = ["Gene name", "Id", "Product",
                     "dN", "dS", "dN/dS", "Selection Type"]
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
                new_seq += seq[ind]+seq[ind+1]+seq[ind+2]
                ind += 3
        return new_seq

    def validate_Seq(self, seq):

        pattern = r'[^ACGTacgt]'  # pattern of non-nucleotide
        test_string = str(seq)
        result = re.search(pattern, test_string)
        if result != None:
            return None
        return "true"


if __name__ == "__main__":

    # "Data\BS168.gb"
    # Q2:
    covid_2021 = 'Data\MZ383039.1.gb'  # corona virus from 2021
    covid_2022 = 'Data\OQ065689.1.gb'  # corona virus from 2022
    va = Virus_Analyzer(covid_2021, covid_2022)

    # Q1
    print("\nQ1:\n")
    va.calc_synonymous()
    # Q2_a
    print("\nQ2 - A:\n")
    va.compare_gb_files()
    # Q2_b
    print("Q2 - B:\n")
    va.compare_genes()
    # Q2_c
    print("Q2 - C:\n")
    va.common_genes_analysis()
