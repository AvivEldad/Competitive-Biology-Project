from Bio import Entrez, SeqIO
import os

class Virus_Analyzer:

    # def __init__(self):

    #     pass
    #     # self.CV21_record_gb = self.read_gb_file(covid_2021)
    #     # self.CV22_record_gb = self.read_gb_file(covid_2022)

    # Q2_a_b
    def compare_gb_files(self, files):
        assert (os.path.exists(files[0]))  # making sure that the path is valid
        assert (os.path.exists(files[1]))  # making sure that the path is valid
        self.gene_lists = {}

        for file in files:
            self.cur_file = file
            file_name = os.path.basename(file)
            record_gb = self.read_gb_file()
            self.cur_features_list = record_gb.features
            num_genes = self.count_genes()
            num_proteins = self.count_proteins()

            print('{}\n=============\nGenes Number: {}\nProteins Number: {}\n'.format(
                file_name, num_genes, num_proteins))
            self.gene_lists[file_name] = self.build_genes_lst()

    def compare_genes(self):

        formate_massage = 'the values that only in {} and not in {} are - {}\n'
        for i in range(0, len(self.gene_lists)-1):

            list1 = list(self.gene_lists.values())[i]
            list2 = list(self.gene_lists.values())[i+1]
            file_name1 = list(self.gene_lists)[i]
            file_name2 = list(self.gene_lists)[i+1]

            common_genes_lst = list(set(list1).intersection(list2))

            res = list(set(list1) - set(list2))
            if len(res) != 0:
                print(formate_massage.format(file_name1, file_name2, res))

            res = list(set(list2) - set(list1))
            if len(res) != 0:
                print(formate_massage.format(file_name2, file_name1, res))

        print('There are {} common values between the files.\n the values are: {}\n'.format(
            len(common_genes_lst), common_genes_lst))
        self.common_genes_lst = common_genes_lst

    def build_genes_lst(self):
        lst = []
        for f in self.cur_features_list:
            if 'gene' in f.qualifiers.keys():
                lst.append(f.qualifiers['gene'][0])
        return lst

    def read_gb_file(self):

        # read file:
        with open(self.cur_file, "r") as input_handle:
            gen = SeqIO.parse(input_handle, "genbank")
            record_gb = next(gen)  # content of 1st record
        return record_gb

    def count_genes(self):

        num_genes = 0
        for i in range(len(self.cur_features_list)):
            f = self.cur_features_list[i]
            if f.type == 'gene':
                num_genes += 1

        return num_genes

    def count_proteins(self):
        num_cds = 0
        for i in range(len(self.cur_features_list)):
            f = self.cur_features_list[i]
            if f.type == 'CDS':
                num_cds += 1

        return num_cds

    def calc_dnds(self, gene):
        pass

    def common_genes_analysis(self):
        gene_lst = self.common_genes_lst[:5]
        print(gene_lst)
        for gene in gene_lst:
            self.calc_dnds(gene)


if __name__ == "__main__":

    va = Virus_Analyzer()
    
    # Q2:
    covid_2021 = 'Data\MZ383039.1.gb'  # corona virus from 2021
    covid_2022 = 'Data\OQ065689.1.gb'  # corona virus from 2022

    # Q2_a
    print("\nQ2 - A:\n")
    va.compare_gb_files([covid_2021, covid_2022])
    # Q2_b
    print("Q2 - B:\n")
    va.compare_genes()
    # Q2_c
    print("Q2 - C:\n")
    va.common_genes_analysis()

