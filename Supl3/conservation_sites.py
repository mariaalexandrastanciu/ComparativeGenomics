from Bio import AlignIO
from collections import Counter
import json

def Simpson_Diversity_score(distinct_residues, total_nr_rows):
    """
    function that calculates the Simpson diversity given the distinct residues, how many times they occur and total
    number of residues in the column; treats gaps as another residue, but marks it with a flag
    :param distinct_residues: dictionary with residue and how many times occurred
    :param total_nr_rows: total number of rows in the column
    :return: Simpson diversity, a flag that states if there are gaps in the column or not(contains_gaps ="N" = no gaps)
    and an array of residues(residues)
    """
    total_sum = 0
    contains_gaps = "N"
    residues = []
    for residue, counts in distinct_residues.items():
        total_sum += counts*(counts-1)
        if residue != "-":
            residues.append(residue)
        if residue == "-" and contains_gaps == "N":
            contains_gaps = "Y"
    D_score = 1 - (total_sum / (total_nr_rows*(total_nr_rows-1)))
    return D_score, contains_gaps, residues

def conserved_sites(filename):
    """
    reads a multiple sequence alignment created by clustal and calculates the Simpson diversity score of each column
    in the alignment; create a file (Conserved_Sites.txt) containing a dictionary with the column number ,score and
     residues array for Simpson Diversity less or equal than 0.2 and without gaps
    without gaps
    :param filename: filename given from clustal multiple alignment sequences
    :return: a dictionary with the column number ,score and residues array for Simpson Diversity less or equal than 0.2
    and without gaps
    """
    alignment = AlignIO.read(filename, "clustal")
    print("Number of rows: %i" % len(alignment))

    all_columns = {}
    conserved_sites_dict = {}
    for i in range(len(alignment[1, :])):
        distinct_residues = Counter(alignment[:, i])
        # for each cloumn number distinct resdues are calculated and Simpson Diversity is calculated
        D_score, contains_gaps, residues = Simpson_Diversity_score(distinct_residues, len(alignment))
        all_columns[i] = [D_score, contains_gaps, residues]
    for column, score in all_columns.items():
        if score[0] < 0.2 and score[1] == "N":
            conserved_sites_dict[column] = [score[0], score[2]]
    with open("Conserved_Sites.txt", 'w') as file:
        file.write("Conserved sites in multiple sequence alignment  " + "\n\n")
        file.write(json.dumps(conserved_sites_dict))
    return conserved_sites_dict

