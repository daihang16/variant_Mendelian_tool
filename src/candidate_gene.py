#!/usr/bin/env python
# Copyright (c) 2015, Hang Dai <hang.dai@bcm.edu>

import os, sys

def multihit(file_name):
	"""(str) -> NoneType
	read the data in the input_file, write a file called candidate_genes.tsv
	sort genes according to the number of occurences
	"""

	gene_family_dict={} #dictionary
	with open(file_name,'r') as input_file:
		for line in input_file:
			info=line.split()
			gene=info[0]
			family_id=info[1]
			if gene in gene_family_dict:
				gene_family_dict[gene].append(family_id)
			else:
				gene_family_dict[gene]=[family_id]
	
	with open('candidate_genes_unordered.tsv','w') as output_file:
		for key in gene_family_dict:
			new_line=key+'\t'+', '.join(list(set(gene_family_dict[key])))+'\t'+str(len(list(set(gene_family_dict[key]))))+'\n'
			output_file.write(new_line)
	
	os.system("""sort -g -r -k3,3 -t '\t' candidate_genes_unordered.tsv >candidate_genes_ordered.tsv""")
	os.system("""awk 'BEGIN {print "gene\tfamily_id\toccurence_number"};{print}' candidate_genes_ordered.tsv >candidate_genes.tsv""")
	os.system('rm -f candidate_genes_unordered.tsv candidate_genes_ordered.tsv')
	
if __name__=='__main__':
	file_name=sys.argv[1]
	multihit(file_name)
