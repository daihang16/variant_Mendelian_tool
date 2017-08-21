#!/usr/bin/env python
# Copyright (c) 2015, Hang Dai <hang.dai@bcm.edu>

import sys, logging

#choose damaging SNVs:
def SNV(file_name):
	"""(str) -> NoneType
	read the data in the input_file, write a file called input_file.parsed
	annotate each SNV in the input_file
	"""
	with open(file_name+'.parsed','w') as output_file:
		with open(file_name,'r') as read_file:
			for line in read_file:
				if line[0:3].upper()!="CHR": #check header, if there is header, skip it
					available_results=0
					info=line.split()  #list
					for i in range(4,9):  #accumulate the fields of SIFT, LRT, FATHMM, MutationTaster, MutationAssessor
						if info[i]!='.':
							available_results=available_results+1
					if info[9]!='.' or info[10]!='.':  #accumulate the fields of PolyPhen2_HDIV and PolyPhen2_HVAR, in total we have 6 prediction tools
						available_results=available_results+1
					for i in range(11,14):  #accumulate the fields of PROVEAN, MetaLR, MetaSVM
						if info[i]!='.':
							available_results=available_results+1
					#calculate the number of damaging annotations
					damaging_number=0
					if info[4]=='D':  #SIFT_pred=D
						damaging_number=damaging_number+1
					if info[5]=='D':  #LRT_pred=D
						damaging_number=damaging_number+1
					if info[6]=='D':  #FATHMM_pred=D
						damaging_number=damaging_number+1
					if info[7]=='D' or info[7]=='A':  #MutationTaster_pred=D or A
						damaging_number=damaging_number+1
					if info[8]=='H' or info[8]=='M':  #MutationAssessor_pred=H or M
						damaging_number=damaging_number+1
					if info[9]=='D' or info[9]=='P' or info[10]=='D' or info[10]=='P':  #any D or P in PolyPhen2_HDIV_pred and PolyPhen2_HVAR_pred
						damaging_number=damaging_number+1
					if info[11]=='D':  #PROVEAN_pred=D
						damaging_number=damaging_number+1
					if info[12]=='D':  #MetaLR_pred=D
						damaging_number=damaging_number+1
					if info[13]=='D':  #MetaSVM_pred=D
						damaging_number=damaging_number+1
					#write new line
					variant=info[0:4] #list [chr,pos,ref,alt]
					#check how much percentage of available tools predicting the variant as "damaging"
					if available_results!=0 and float(damaging_number)/available_results>=0.2:
						variant.append('damaging_SNV')  #VMT_annotation
						variant.append(info[7])  #VMT_MutationTaster_pred
					else:
						variant.append('.')  #VMT_annotation
						variant.append(info[7])  #VMT_MutationTaster_pred
					new_line='\t'.join(variant)+'\n'
					output_file.write(new_line)
	
#choose damaging INDELs:
def INDEL(file_name):
	"""(str) -> NoneType
	read the data in the input_file, write a file called input_file.parsed
	annotate each INDEL in the input_file
	"""
	variant_pred_dict={}  #dictionary, keys are variants, values are lists with predictions
	
	with open(file_name+'.parsed','w') as output_file:
		with open(file_name,'r') as read_file:
			read_file.readline()  #read header
			for line in read_file:
				try:
					info=line.split('\t')  #list
					if info[0].strip()=='23':
						chromosome='X'
					elif info[0].strip()=='24':
						chromosome='Y'
					else:
						chromosome=info[0].strip()
					position=info[1].strip()
					prediction=info[3].strip()
					if info[9].strip()=='':
						ref='-'
					else:
						ref=info[9].strip()
					if info[10].strip()=='':
						alt='-'
					else:
						alt=info[10].strip()
					variant=chromosome+'_'+position+'_'+ref+'_'+alt
					if variant in variant_pred_dict:
						variant_pred_dict[variant].append(prediction)
					else:
						variant_pred_dict[variant]=[prediction]
				except IndexError:
					continue
					
		for item in variant_pred_dict:
			variant=item.split('_') #list
			if 'disease_causing_automatic' in variant_pred_dict[item]:  #known to be deleterious
				variant.append('damaging_by_MT')  #VMT_annotation
				variant.append('A')  #VMT_MutationTaster_pred
			elif 'disease_causing' in variant_pred_dict[item]: #probably deleterious
				variant.append('damaging_by_MT')  #VMT_annotation
				variant.append('D')  #VMT_MutationTaster_pred
			elif 'polymorphism_automatic' in variant_pred_dict[item]: #known to be harmless
				variant.append('nondamaging_by_MT')  #VMT_annotation
				variant.append('P')  #VMT_MutationTaster_pred
			elif 'polymorphism' in variant_pred_dict[item]: #probably harmless
				variant.append('nondamaging_by_MT')  #VMT_annotation
				variant.append('N')  #VMT_MutationTaster_pred
			else:  #in case there is sth else
				variant.append('nondamaging_by_MT')  #VMT_annotation
				variant.append('.')  #VMT_MutationTaster_pred
			new_line='\t'.join(variant)+'\n'
			output_file.write(new_line)
	
if __name__=='__main__':
	#set logging:
	logger = logging.getLogger()
	handler = logging.StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	handler.setFormatter(formatter)
	logger.addHandler(handler)
	logger.setLevel(logging.INFO)
	#file name and variant type:
	file_name=sys.argv[1]
	variant_type=sys.argv[2]
	if variant_type.upper()=='DBNSFP':
		logger.info('annotating SNVs with VMT_annotation')
		SNV(file_name)
	elif variant_type.upper()=='MUTATIONTASTER':
		logger.info('annotating variants with MutationTaster annotation')
		INDEL(file_name)
