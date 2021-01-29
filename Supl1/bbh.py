import json
import re
import ast

def retrieve_best_hist(input_filename, output_filename):
    """
    function that takes from an input file in format of BLASTP search output best score from each query
    :param input_filename: output from BLASTP search
    :param output_filename: homologs genes with best score
    :return: homologs genes with best score
    """
    best_hits = {}
    with open(input_filename, 'r') as file:
        for line in file:
            if line.startswith("Query="):
                gene = re.search("(?<=\|)(.*?)(?=\|)", line).group()
                best_hit_only = 1
            if line.startswith("tr|") or line.startswith("sp|"):
                if best_hit_only == 1:
                    hit = re.search("(?<=\|)(.*?)(?=\|)", line).group()
                    score = re.search("(?<=\s\s)(.*)(?=\s...)", line).group()
                    if gene != hit:
                        best_hit_only = 0
                        best_hits[gene] = [hit, score.strip()]
    with open(output_filename, 'w') as file:
        file.write(json.dumps(best_hits))
    print("Done writing file: ", output_filename)
    return best_hits


def bbh_2species(species1_file, species2_file, itself_species1_file, itself_species2_file):
    """
    finds best bi-directional hits between species1 and species two taking in consideration in-paralogs
    :param species1_file: txt file with best hits from species1 into species2
    :param species2_file: txt file with best hits from species2 into species1
    :param itself_species1_file: txt file with best hits from species1 into species1 used for in-paralogs determination
    :param itself_species2_file: txt file with best hits from species2 into species2 used for in-paralogs determination
    :return: a dictionary of form {SpeciesC.elegans_ortholog1: [SpeciesD.pulex_ortholog1, [SpeciesC.elegans_ortholog1,
                                SpeciesC.elegans_paralog1,..],[SpeciesD.pulex_ortholog1,SpeciesD.pulex_paralog1,..] ] }
    """
    file = open(species1_file, "r")
    contents = file.read()
    species1 = ast.literal_eval(contents)
    file.close()

    file = open(species2_file, "r")
    contents = file.read()
    species2 = ast.literal_eval(contents)
    file.close()

    ###inparalogs:
    file = open(itself_species1_file, "r")
    contents = file.read()
    itself_species1 = ast.literal_eval(contents)
    file.close()

    file = open(itself_species2_file, "r")
    contents = file.read()
    itself_species2 = ast.literal_eval(contents)
    file.close()

    bbh = {}
    i=0
    j=0
    for sp1_gene in species1.keys():
        sp_bh1 = species1[sp1_gene][0]
        if sp_bh1 in species2.keys():  #species1 value is in species2 key
            if sp1_gene == species2[sp_bh1][0]:  #species 1 key is equal with species 2 value
                in_paralog_sp1 = get_inparalogs(sp1_gene, species1[sp1_gene][1], itself_species1)
                in_paralog_sp2 = get_inparalogs(sp_bh1, species2[sp_bh1][1], itself_species2)
                if len(in_paralog_sp1) > 1:
                    i = i+1
                if len(in_paralog_sp2) > 1:
                    j = j+1
                bbh[sp1_gene] = [sp_bh1, in_paralog_sp1, in_paralog_sp2]

    with open("Ce_Dp_orthologs.txt", 'w') as file:
        file.write("Number of protein coding genes found in C. elegans : " + str(len(itself_species1)) + "\n")
        file.write("Number of protein coding genes found in P.duplex : " + str(len(itself_species2)) + "\n")
        file.write("Number of orthologs: " + str(len(bbh)) + "\n")
        file.write("Number of C.elegans genes with paralogs: " + str(i) + "\n")
        file.write("Number of D.pulex genes with paralogs: " + str(j) + "\n\n")
        file.write(json.dumps(bbh))
    return bbh


def get_inparalogs(gene, score, search_space):
    """
    function to check the existence of inparalogs genes for a given gene and a search space
    :param gene: the gene for which we search paralogs
    :param score: the score of the gene
    :param search_space: the dictionary in which the search is performed; has the form: {gene:[best_hit_gene,score]}
    :return: returns the paralog gene if the score in search space is higher than the score from ortholog search
    """
    inparalog_genes = [gene]
    if gene in search_space.keys():
        itself_score = search_space[gene][1]
        if float(score) < float(itself_score):
            inparalog_genes = [gene, search_space[gene][0]]
    return inparalog_genes
