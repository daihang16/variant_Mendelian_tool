#!/usr/bin/env python
# Copyright (c) 2015, Hang Dai <hang.dai@bcm.edu>

import sys, os, sqlite3, logging, multiprocessing, contextlib, pkg_resources
import Command
import choose_damaging_variants
import moi
import candidate_gene

#number of processors
nproc=multiprocessing.cpu_count()
#nproc=psutil.cpu_count

#set logging
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

	
##init project
def init(args):
	logger.info('initiating project {}'.format(args.proj))
	Command.Init(Project=args.proj, Force=True).Run()


##import vcf and phenotype, index genotype, export INDELs
def import_data(args):
	#import vcf file:
	logger.info('importing vcf file')
	if args.format: #if user specifies a format file, use that file
		format_file=args.format
	else:  #if user does not specify format file, use the default
		format_file=pkg_resources.resource_filename(__name__, '../../../../format_files/vcf.fmt')
	Command.Import(Input_files=args.vcf, Format=format_file, Build=args.build, Jobs=nproc).Run()
	#import phenotype file:
	logger.info('importing family information file')
	Command.Phenotype(From_file=[args.sample_info], Jobs=nproc).Run()
	#index genotype:
	logger.info('indexing genotype')
	logger.warning('This may take a while, do not close the Bash session')
	Command.Export(Table='variant',  Output='tobedeleted.vcf', Samples=['1'], Format='vcf', Jobs=nproc).Run()
	os.remove('tobedeleted.vcf')
	#export INDELs for future use:
	Command.Select(From_table='variant', Condition=["ref='-' or alt='-'"], To_table=['INDEL']).Run()
	Command.Export(Table='INDEL', Output='indel.vcf', Jobs=nproc).Run()
	Command.Remove(Type='tables', Items=['INDEL']).Run()
	#execute ANNOVAR:
	logger.info('executing ANNOVAR')
	logger.info('Please cite the following paper: ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Wang K, Li M, Hakonarson H')
	logger.warning('This may take a while, do not close the Bash session')
	Command.Execute(Specfile='ANNOVAR', Jobs=nproc).Run()
	#use annotations:
	logger.info('using refGene')
	Command.Use(Source='refGene').Run()
	logger.info('using dbNSFP')
	logger.warning('This may take a while if this is the first time you use VMT, do not close the Bash session while downloading the dbNSFP database')
	Command.Use(Source='dbNSFP').Run()
	logger.info('using ExAC')
	logger.warning('This may take a while if this is the first time you use VMT, do not close the Bash session while downloading the ExAC database')
	Command.Use(Source='../testvtools/ExAC').Run()  #fix me 
	
	
##use annotations, output SNV, calculate "damaging" SNVs and INDELs, combine SNVs and INDELs annotation, update
def annotate(args):
	if args.database: #use the databases specified by user
		Command.Use(Source=args.database, Linked_by=args.linked_by).Run()
	elif args.source:  #build the database from source specified by user. Note that --source and --database are mutually exclusive
		if args.format:  #if user provide an annotation format file
			Command.Use(Source=args.format, Files=args.source, Linked_by=args.linked_by, Rebuild=True, Jobs=nproc).Run()  
		else:
			logger.error('Must provide an annotation format file with name *.ann')
			raise SystemExit(0)
	if args.vmt_annotation:  #then do vmt_annotation
		#output SNVs:
		logger.info('annotating SNVs with VMT_annotation')
		logger.warning('This may take a while, do not close the Bash session')
		with open('snv.txt','w') as output_file:
			with stdout_redirect(output_file):
				Command.Select(From_table='variant', Condition=['dbNSFP.chr is not null'], Output=['chr', 'pos', 'ref', 'alt', 'sift_pred', 'lrt_pred', 'fathmm_pred', 'mutationtaster_pred', 'mutationassessor_pred', 'polyphen2_hdiv_pred', 'polyphen2_hvar_pred', 'provean_pred', 'MetaLR_pred', 'MetaSVM_pred']).Run()
		#annotate SNVs:
		choose_damaging_variants.SNV('snv.txt')
		try:
			#annotate INDELs:
			logger.info('annotating INDELs with VMT_annotation')
			choose_damaging_variants.INDEL('indel.tsv')
		except:
			#fail to annotate INDELs:
			logger.warning('fail to annotate INDELs')
		#combine the parsed files together:
		os.system('cat *.parsed >snv_indel.parsed')
		#update variants with VMT_annotation
		logger.info('updating variants with VMT_annotation')
		logger.warning('This may take a while, do not close the Bash session')
		format_file=pkg_resources.resource_filename(__name__, '../../../../format_files/VMT_annotation.fmt')
		Command.Update(From_file=['snv_indel.parsed'], Format=format_file, Jobs=nproc).Run()
		os.remove('snv.txt')
		os.remove('indel.tsv.parsed')
		os.remove('snv_indel.parsed')
		os.remove('snv.txt.parsed')

##filter
def filtering(args):
	#Filter variants according to annotation fields
	annotation_condition=[item for item in [args.pop_maf, args.pos, args.damaging, args.other] if item is not None]  #if a user does not use an argument, drop it
	if args.filter_pass:
		annotation_condition.append('variant.filter="PASS"')
	#Filter variants according to sample genotypes
	sample_condition=[item for item in [args.samples] if item is not None]  #if a user does not use an argument, drop it
	
	Command.Select(From_table=args.from_table, Condition=annotation_condition, Samples=sample_condition, To_table=[args.to_table]).Run()

##choose
def choose(args):
	families=[]  #list, store family_id in the file provided by user
	with open(args.family_id,'r') as input_file:
		for line in input_file:
			families.append(line.strip())
	for family_id in families:
		with open('{}_samples'.format(family_id),'w') as output_file:
			with stdout_redirect(output_file):
				try:
					Command.Phenotype(Samples=['family_id="{}"'.format(family_id)], Output=['sample_name', 'affection_status', 'unaffected_parental_genotypes']).Run()
				except:
					Command.Phenotype(Samples=['family_id="{}"'.format(family_id)], Output=['sample_name', 'affection_status']).Run()
		affected_individuals=[]#list
		unaffected_individuals=[]#list
		unaffected_parental_genotypes=''  #str
		with open('{}_samples'.format(family_id),'r') as input_file:
			for line in input_file:
				info=line.split()
				try:
					unaffected_parental_genotypes=info[2]
				except:
					unaffected_parental_genotypes=''
				if info[1]=='2':
					affected_individuals.append(info[0])
				else:
					unaffected_individuals.append(info[0])
		os.remove('{}_samples'.format(family_id))
		
		if args.moi=='AD':
			moi.AD(family_id, affected_individuals, args.from_table)
		elif 'AD_phenocopy_' in args.moi:
			proportion_of_affected_individuals_carrying_variant=float(args.moi.split('_')[2])
			moi.AD_phenocopy(family_id, affected_individuals, args.from_table, proportion_of_affected_individuals_carrying_variant)
		elif args.moi=='De_Novo':
			if unaffected_parental_genotypes.upper()=='YES':
				moi.De_Novo(family_id, affected_individuals, args.from_table)
			else: #otherwise no point to check for De Novo variants, raise an error and continue
				logger.error('unaffected_parental_genotypes not available for family {0}\ncannot check for De Novo variants'.format(family_id))
				continue
		elif args.moi=='AR':
			if unaffected_parental_genotypes.upper()=='YES':  #then will check for segregation
				moi.AR_check_segregation(family_id, affected_individuals, unaffected_individuals, args.from_table)
			else:
				moi.AR(family_id, affected_individuals, args.from_table)
		#output the variants for user
		logger.info('putting out variants in table {0}_{1}_{2}'.format(family_id, args.from_table, args.moi))
		with open('{0}_{1}_{2}.tsv'.format(family_id, args.from_table, args.moi),'w') as output_file:  #replace any old file
			with stdout_redirect(output_file):
				Command.Select(From_table='{0}_{1}_{2}'.format(family_id, args.from_table, args.moi), Output=['variant_id', 'chr', 'raw_pos', 'raw_ref', 'raw_alt', "samples('geno_filter=GT=1')", "samples('geno_filter=GT=2')", 'refgene.name2', 'variant.filter', 'region_type', 'region_name', 'mut_type', 'function', 'ExAC.allaltcount', 'ExAC.all_maf', 'ExAC.afraltcount', 'ExAC.afr_maf', 'ExAC.amraltcount', 'ExAC.amr_maf', 'ExAC.easaltcount', 'ExAC.eas_maf', 'ExAC.finaltcount', 'ExAC.fin_maf', 'ExAC.nfealtcount', 'ExAC.nfe_maf', 'ExAC.sasaltcount', 'ExAC.sas_maf', 'dbNSFP.sift_pred', 'dbNSFP.lrt_pred', 'dbNSFP.fathmm_pred', 'VMT_MutationTaster_pred', 'dbNSFP.mutationassessor_pred', 'dbNSFP.polyphen2_hdiv_pred', 'dbNSFP.polyphen2_hvar_pred', 'dbNSFP.provean_pred', 'dbNSFP.MetaLR_pred', 'dbNSFP.MetaSVM_pred', 'dbNSFP.cadd_phred'], Header=['variant_id', 'chr', 'pos', 'ref', 'alt', 'het_samples', 'homo_samples', 'gene_name', 'filter', 'region_type', 'region_name', 'mut_type', 'function', 'ExAC_all_alt_count', 'ExAC_all_maf', 'ExAC_afr_alt_count', 'ExAC_afr_maf', 'ExAC_amr_alt_count', 'ExAC_amr_maf', 'ExAC_eas_alt_count', 'ExAC_eas_maf', 'ExAC_fin_alt_count', 'ExAC_fin_maf', 'ExAC_nfe_alt_count', 'ExAC_nfe_maf', 'ExAC_sas_alt_count', 'ExAC_sas_maf', 'SIFT', 'LRT', 'FATHMM', 'MutationTaster', 'MutationAssessor', 'Polyphen2_HDIV', 'Polyphen2_HVAR', 'PROVEAN', 'MetaLR', 'MetaSVM', 'CADD'], Delimiter='\t').Run()

		#output gene name and family_id into a file for future use if user wants multihit
		with open('{0}_{1}.temp.genename'.format(args.from_table, args.moi),'w') as output_file:
			with stdout_redirect(output_file):
				Command.Select(From_table='{0}_{1}_{2}'.format(family_id, args.from_table, args.moi), Condition=['refGene.name2 is not null'], Output=['refGene.name2']).Run()
		with open('{0}_{1}.genename'.format(args.from_table, args.moi),'a') as output_file:
			with open('{0}_{1}.temp.genename'.format(args.from_table, args.moi),'r') as input_file:
				for line in input_file:
					if line.strip()!='refGene_name2':  #skip header
						output_file.write(line.strip()+'\t'+family_id+'\n')
	#delete the temp file
	os.remove('{0}_{1}.temp.genename'.format(args.from_table, args.moi))
					
	if args.multihit: #if user wants to know the replicated genes
		candidate_gene.multihit('{0}_{1}.genename'.format(args.from_table, args.moi))
		#then remove the file to avoid confusion and future possible appending lines to the file
		os.remove('{0}_{1}.genename'.format(args.from_table, args.moi))
	else: #remove the file to avoid confusion and future possible appending lines to the file
		os.remove('{0}_{1}.genename'.format(args.from_table, args.moi))
	

##output
def output(args):
	Command.Output(Table=args.from_table, Fields=args.fields, Header=args.header, Group_by=args.group_by, Delimiter=args.delimiter, Order_by=args.order_by).Run()
	
##show
def show(args):
	Command.Show(Type=args.type, Items=args.item).Run()
	
##remove
def remove(args):
	Command.Remove(Type=args.type, Items=args.item).Run()
	
##compare
def compare(args):
	Command.Compare(Tables=args.tables, Union=args.union, Intersection=args.intersection, Difference=args.difference, Expression=args.expression).Run()
	

