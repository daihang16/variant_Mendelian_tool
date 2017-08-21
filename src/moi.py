#!/usr/bin/env python
# Copyright (c) 2015, Hang Dai <hang.dai@bcm.edu>

import sys, os, contextlib, logging, multiprocessing
import Command

#number of processors
nproc=multiprocessing.cpu_count()

#set logging:
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

##Utils
@contextlib.contextmanager
def stdout_redirect(where):
    sys.stdout = where
    try:
        yield where
    finally:
        sys.stdout = sys.__stdout__

##function to find variants in autosomal dominant assumption:
def AD(family_id, affected_individuals, table_name):
	"""(str, list, str) -> NoneType
	get the family_id, a list of affected individuals, a table_name you want to analyze
	then generate a table with variants
	"""
	genotypes=[]
	for item in affected_individuals:
		genotypes.append('genotype("{}")=1'.format(item))  #genotype("")=1, checking for heterozygous variants.
		
	logger.info('selecting variants for family {} from table {} according to autosomal dominant MOI'.format(family_id, table_name))
	#select shared heterozygous variants in affected individuals.
	Command.Select(From_table=table_name, Condition=genotypes, To_table=["{}_{}_AD".format(family_id,table_name)]).Run()
	
##function to find variants in autosomal dominant with phenocopy assumption, that is, some affected individuals do not share variant:
def AD_phenocopy(family_id, affected_individuals, table_name, proportion_of_affected_individuals_carrying_variant):
	"""(str, list, str, float) -> NoneType
	get the family_id, a list of affected individuals, a table_name you want to analyze
	then generate a table with variants
	"""
	#select heterozygous variants in affected individuals, even if the variant exists in only one affected individual.
	Command.Select(From_table=table_name, Samples=['family_id="{}" and affection_status=2'.format(family_id)], To_table=["{}_affected".format(family_id)]).Run()
	#count the number of heterozygotes for that small number of variants, for affected individuals in that family only
	Command.Update(Table="{}_affected".format(family_id), Jobs=nproc, From_stat=['het_count_temp=#(het)'], Samples=['family_id="{}" and affection_status=2'.format(family_id)]).Run()
	#if a variant occurs several times among affected individuals, more than designated proportion, select the variant
	designated_number_of_affected_individuals=len(affected_individuals)*proportion_of_affected_individuals_carrying_variant
	logger.info('selecting variants for family {} from table {} according to autosomal dominant MOI with phenocopy'.format(family_id, table_name))
	Command.Select(From_table="{}_affected".format(family_id), Condition=['het_count_temp>={}'.format(str(designated_number_of_affected_individuals))], To_table=["{}_{}_AD_phenocopy_{}".format(family_id,table_name,str(proportion_of_affected_individuals_carrying_variant))]).Run()
	#delete the temp field and temp table
	Command.Remove(Type='tables', Items=["{}_affected".format(family_id)]).Run()
	Command.Remove(Type='fields', Items=['het_count_temp']).Run()
	
##function to find variants in De Novo assumption:
def De_Novo(family_id, affected_individuals, table_name):
	"""(str, list, str) -> NoneType
	get the family_id, a list of affected individuals, a table_name you want to analyze
	then generate a table with variants
	"""
	genotypes=[]
	for item in affected_individuals:
		genotypes.append('genotype("{}")=1'.format(item))  #genotype("")=1, checking for heterozygous variants.
		
	logger.info('selecting variants for family {} from table {} according to De Novo MOI'.format(family_id, table_name))
	#first we select the heterozygous variants shared by affected individuals or just owned by one affected individuals
	Command.Select(From_table=table_name, Condition=genotypes, To_table=["{}_affected_shared_het".format(family_id)]).Run()
	#Then we select all the variants in unaffected individuals
	Command.Select(From_table=table_name, Condition=genotypes, Samples=['family_id="{}" and affection_status=1'.format(family_id)], To_table=["{}_unaffected".format(family_id)]).Run()
	#Then we shall exclude any variant in unaffected individuals from heterozygous variants in affected individual
	Command.Compare(Expression=['{0}_{1}_De_Novo={0}_affected_shared_het-{0}_unaffected'.format(family_id, table_name)]).Run()
	Command.Remove(Type='tables', Items=['{}_affected_shared_het'.format(family_id), '{}_unaffected'.format(family_id)]).Run()

##function to find variants in autosomal recessive assumption:
def AR(family_id, affected_individuals, table_name):
	"""(str, list, str) -> NoneType
	get the family_id, a list of affected individuals, a table_name you want to analyze
	then generate a table with variants
	"""
	genotypes_cph=[]  #this list is used for searching for compound het variants
	genotypes_hom=[]  #this list is used for searching for homozygous variants
	for item in affected_individuals:
		genotypes_cph.append('genotype("{}")=1'.format(item))  #genotype("")=1, checking for compound heterozygous variants.
		genotypes_hom.append('genotype("{}")=2'.format(item))  #genotype("")=2, checking for homozygous variants.
		
	logger.info('selecting variants for family {} from table {} according to autosomal recessive MOI'.format(family_id, table_name))
	#Find the shared genes and heterozygous variants by all affected individuals
	logger.info('checking for possible compound heterozygous variants')
	with open('{}_cph_genes'.format(family_id),'w') as output_file:
		with stdout_redirect(output_file):
			Command.Select(From_table=table_name, Condition=genotypes_cph, Output=['refGene.name2', 'count()'], Group_by=['refGene.name2']).Run()
	family_cph_genes=[]  #list to store compound het genes
	with open('{}_cph_genes'.format(family_id),'r') as input_file:
		for line in input_file:
			gene=line.split()[0]
			count=line.split()[1]
			try:
				if gene!='.' and int(count)>1:
					family_cph_genes.append('refGene.name2="{}"'.format(gene))
			except ValueError: #skip the header: refgene_name2   count__
				continue
	#delete the temp file
	os.remove('{}_cph_genes'.format(family_id))
	
	if len(family_cph_genes)!=0: #there are compound het variants
		logger.info('there are possible compound heterozygous variants')
		#select compound het variants for the family
		logger.info('selecting possible compound heterozygous variants')
		Command.Select(From_table=table_name, Condition=[' or '.join(family_cph_genes), ' and '.join(genotypes_cph)], To_table=['{}_cph'.format(family_id)]).Run()
		#nowthat we got compound het genes and variants, let's get homozygous variants then
		logger.info('selecting homozygous variants')
		Command.Select(From_table=table_name, Condition=[' and '.join(genotypes_hom)], To_table=['{}_hom'.format(family_id)]).Run()
		#let's combine the family_id_cph and family_id_hom tables to form a new table: family_id_AR
		Command.Compare(Expression=['{0}_{1}_AR={0}_cph|{0}_hom'.format(family_id, table_name)]).Run()
		Command.Remove(Type='tables', Items=['{}_cph'.format(family_id), '{}_hom'.format(family_id)]).Run()
	else: #there is no compound het variants, select hom variants directly
		logger.info('there are no possible compound heterozygous variants')
		logger.info('selecting homozygous variants')
		Command.Select(From_table=table_name, Condition=[' and '.join(genotypes_hom)], To_table=['{0}_{1}_AR'.format(family_id, table_name)]).Run()
		
def AR_check_segregation(family_id, affected_individuals, unaffected_individuals, table_name):
	"""(str, list, list, str) -> NoneType
	get the family_id, a list of affected individuals, a list of unaffected individuals (2 unaffected parents), a table_name you want to analyze
	then generate a table with variants
	"""
	affected_genotypes_cph=[]  #this list is used for searching for compound het variants in affected individuals
	affected_genotypes_hom=[]  #this list is used for searching for homozygous variants in affected individuals
	for item in affected_individuals:
		affected_genotypes_cph.append('genotype("{}")=1'.format(item))  #genotype("")=1, checking for compound heterozygous variants.
		affected_genotypes_hom.append('genotype("{}")=2'.format(item))  #genotype("")=2, checking for homozygous variants.
		
	unaffected_genotypes_het=[]  #this list is used for searching for heterozygous variants in unaffected individuals (parents)
	for item in unaffected_individuals:
		unaffected_genotypes_het.append('genotype("{}")=1'.format(item))  #genotype("")=1, checking for heterozygous variants.
		
	logger.info('selecting variants for family {} from table {} according to autosomal recessive MOI'.format(family_id, table_name))
	#Find the shared genes and heterozygous variants by all affected individuals
	logger.info('checking for possible compound heterozygous variants')
	with open('possible_{}_cph_genes'.format(family_id),'w') as output_file:
		with stdout_redirect(output_file):
			Command.Select(From_table=table_name, Condition=affected_genotypes_cph, Output=['refGene.name2', 'count()'], Group_by=['refGene.name2']).Run()
	possible_family_cph_genes=[]  #list to store possible compound het genes
	with open('possible_{}_cph_genes'.format(family_id),'r') as input_file:
		for line in input_file:
			gene=line.split()[0]
			count=line.split()[1]
			try:
				if gene!='.' and int(count)>1:
					possible_family_cph_genes.append('refGene.name2="{}"'.format(gene))
			except ValueError: #skip the header: refgene_name2   count__
				continue
	#delete the temp file
	os.remove('possible_{}_cph_genes'.format(family_id))
	
	#now check segregation
	seg_cph_variants=[] #list to store segregated compound het variants
	for gene in possible_family_cph_genes:  #for each gene
		possible_parents_cph_variants=[]  #list containing 2 lists, each for a parent
		for unaffected in unaffected_genotypes_het:  #for each unaffected individual (parent)
			parent=[]  #list
			with open('temp','w') as output_file:
				with stdout_redirect(output_file):
					Command.Select(From_table=table_name, Condition=[gene, unaffected, ' and '.join(affected_genotypes_cph)], Output=['variant_id']).Run()
			with open('temp','r') as input_file:
				for line in input_file:
					if line.strip()!='variant_id':  #skip header
						parent.append(line.strip())
			possible_parents_cph_variants.append(parent)
		if len(possible_parents_cph_variants[0])>0 and len(possible_parents_cph_variants[1])>0:  #both len>0
			gene_seg_cph_variants=list(set(possible_parents_cph_variants[0]+possible_parents_cph_variants[1]))
			if len(gene_seg_cph_variants)>1:  #set len>1
				for item in gene_seg_cph_variants:
					seg_cph_variants.append("variant_id={}".format(item))
	#finish checking

	if len(seg_cph_variants)!=0: #there are (possibly) segregated compound het variants
		logger.info('there are possible compound heterozygous variants')
		#select compound het variants for the family
		logger.info('selecting possible compound heterozygous variants')
		Command.Select(From_table=table_name, Condition=[' or '.join(seg_cph_variants)], To_table=['{}_cph'.format(family_id)]).Run()
		#nowthat we got compound het genes and variants, let's get homozygous variants then
		logger.info('selecting homozygous variants')
		Command.Select(From_table=table_name, Condition=[' and '.join(affected_genotypes_hom), ' and '.join(unaffected_genotypes_het)], To_table=['{}_hom'.format(family_id)]).Run()
		#let's combine the family_id_cph and family_id_hom tables to form a new table: family_id_AR
		Command.Compare(Expression=['{0}_{1}_AR={0}_cph|{0}_hom'.format(family_id, table_name)]).Run()
		Command.Remove(Type='tables', Items=['{}_cph'.format(family_id), '{}_hom'.format(family_id)]).Run()
	else: #there is no compound het variants, select hom variants directly
		logger.info('there are no possible compound heterozygous variants')
		logger.info('selecting homozygous variants')
		Command.Select(From_table=table_name, Condition=[' and '.join(affected_genotypes_hom), ' and '.join(unaffected_genotypes_het)], To_table=['{0}_{1}_AR'.format(family_id, table_name)]).Run()
		
if __name__=='__main__':
	table_name=sys.argv[1] #str
	affected_individuals=[] #list
	for line in sys.stdin:
		info=line.split()
		affected_individuals.append(info[0])
		family_id=info[1] #str
	if sys.argv[2].upper()=='AD':
		AD(family_id, affected_individuals, table_name)
	elif sys.argv[2].upper()=='AR':
		AR(family_id, affected_individuals, table_name)
	elif sys.argv[2].upper()=='DE_NOVO':
		De_Novo(family_id, affected_individuals, table_name)
	else:
		logger.error('Please choose the correct MOI: AD, AR, De_Novo')
