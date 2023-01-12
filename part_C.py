from Bio import SeqIO
from Bio.Seq import Seq
import os
import general as ge
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio import Entrez, SeqIO, Align
from Bio.Align import substitution_matrices
import numpy as np
from tabulate import tabulate
import re


class Virus_Analyzer:

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
    def compare_gb_files(self, files):
        assert (os.path.exists(files[0]))  # making sure that the path is valid
        assert (os.path.exists(files[1]))  # making sure that the path is valid
        self.gene_lists = {}

        for file in files:
            self.cur_file = file
            file_name = os.path.basename(file)
            gene_lst = self.count_genes(file)
            num_proteins = self.count_proteins(file)

            print('{}\n=============\nGenes Number: {}\nProteins Number: {}\n'.format(
                file_name, len(gene_lst), num_proteins))
            self.gene_lists[file_name] = gene_lst

    def count_genes(self, file):
        gene_list = []
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                gene = feature.qualifiers.get("gene", None)
                if gene:
                    gene_list.append(gene[0])

        return gene_list

    def count_proteins(self, file):
        num_cds = 0
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                protein_id = feature.qualifiers.get("protein_id", None)
                if protein_id:
                    num_cds = num_cds + 1

        return num_cds

    # Q2_b
    def compare_genes(self):

        formate_massage = 'the values that only in {} and not in {} are - {}\n'
        for i in range(0, len(self.gene_lists)-1):

            list1 = list(self.gene_lists.values())[i]
            list2 = list(self.gene_lists.values())[i+1]
            file_name1 = list(self.gene_lists)[i]
            file_name2 = list(self.gene_lists)[i+1]

            common_genes_lst = list(set(list1) & (set(list2)))
            res = list(set(list1) - set(list2))
            if len(res) != 0:
                print(formate_massage.format(file_name1, file_name2, res))

            res = list(set(list2) - set(list1))
            if len(res) != 0:
                print(formate_massage.format(file_name2, file_name1, res))

        print('There are {} common values between the files.\n the values are: {}\n'.format(
            len(common_genes_lst), common_genes_lst))
        self.common_genes_lst = common_genes_lst

    # Q3_c

    def common_genes_analysis(self, file1, file2):

        gene_lst = []
        table = []
        for gene in self.common_genes_lst:

            seq1, translation_table1, role1, id1 = self.get_sequence(gene, file1)
            seq2, translation_table2, role2, id2 = self.get_sequence(gene, file2)

            if seq1 != None and seq2 != None:
                # the gene must be coded to protein, if not- the seq will be none
                gene_lst.append(gene)
                aligned_seq1, aligned_seq2 = self.trans_and_align(seq1, seq2,
                                                                  translation_table1, translation_table2)


                dN, dS = cal_dn_ds(CodonSeq(aligned_seq1), CodonSeq(aligned_seq2))

                if dS == 0: # there are no synonimus codons
                    dN_dS = 0
                else:
                    dN_dS = dN / dS
                
                dN_dS = float("%0.3f" % dN_dS) #round
                selection = ""
                if dN_dS == 1 :  
                    selection = "neutral"
                elif dN_dS < 1:
                    selection = "negative"
                else:
                    selection = "positive"

                #build table:
                table.append([gene, id1, role1, dN_dS, selection])
          

            # choose only 5
            if len(gene_lst) == 5:
                break
           
        #define header names
        col_names = ["Gene name", "Id", "Product", "dN/dS", "Selection Type"] 
        #display table
        print(tabulate(table, headers=col_names))
        print("\n",gene_lst)

    def get_sequence(self, gene, file):
        # record_gb = ge.parse_genbank(file)
        # features_list = record_gb.features

        # for i in range(len(features_list)) :
        #     f = features_list[i]
        #     if f.type == 'CDS':
        #         gene_name = f.qualifiers["gene"][0]
        #         if gene_name != gene:
        #             continue

        #         if 'transl_table' in f.qualifiers.keys():
        #             transl_table = f.qualifiers['transl_table'][0]
        #         else:
        #             transl_table = "Standard"
        #         if 'product' in f.qualifiers.keys():
        #             gene_role = f.qualifiers['product'][0]
        #         gene_Id = f.id
        #         start_s = f.location.start.position
        #         end_s = f.location.end.position
        #         gene_seq = record_gb.seq[start_s:end_s]
        #         pattern = r'[^ACGTacgt]' # pattern of non-nucleotide
        #         test_string = str(gene_seq)
        #         result = re.search(pattern, test_string)
        #         if result !=None:
        #             continue
        #         return Seq(gene_seq), transl_table, gene_role, gene_Id
        # return None
        

                


        # Open the GenBank file
        with open(file, "r") as handle:
            # Iterate through the records in the file
            for record in SeqIO.parse(handle, "genbank"):
                # Check if the gene name of the feature matches the desired gene
                for feature in record.features:
                    if feature.type == "CDS" and "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name == gene:
                            # Extract the sequence of the gene
                            transl_table = feature.qualifiers.get(
                                "transl_table", "Standard")
                            gene_role = feature.qualifiers.get("product")[0]
                            gene_Id = feature.qualifiers.get("protein_id")[0]
                            gene_seq = feature.extract(record.seq)
                            pattern = r'[^ACGTacgt]' # pattern of non-nucleotide
                            test_string = str(gene_seq)
                            result = re.search(pattern, test_string)
                            if result !=None:
                                continue
                            return Seq(gene_seq), transl_table, gene_role, gene_Id

        return None,None,None,None

    def trans_and_align(self, seq1, seq2, t1, t2):

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        messenger_rna1 = seq1.translate(t1)
        messenger_rna2 = seq2.translate(t2)

        alignments = aligner.align(messenger_rna1, messenger_rna2)
        print("seq: ",seq1, "\n")
        print("trans : ",messenger_rna1,"\n")
        print("ali: ", alignments[0][0][:-1],"\n")


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


if __name__ == "__main__":

    va = Virus_Analyzer()

    # Q1:

    # Q2:
    covid_2021 = 'Data\MZ383039.1.gb'  # corona virus from 2021
    covid_2022 = 'Data\OQ065689.1.gb'  # corona virus from 2022

    # Q1
    print("\nQ1:\n")
    va.calc_synonymous()
    # Q2_a
    print("\nQ2 - A:\n")
    va.compare_gb_files([covid_2021, covid_2022])
    # Q2_b
    print("Q2 - B:\n")
    va.compare_genes()
    # Q2_c
    print("Q2 - C:\n")
    va.common_genes_analysis(covid_2021, covid_2022)
