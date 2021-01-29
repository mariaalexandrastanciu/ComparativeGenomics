import conservation_sites as cs

conserved_sites = cs.conserved_sites("blast_result.aln")
print(conserved_sites)
print(len(conserved_sites))
