#!/usr/bin/python2.7
import os
import json
import re
from sys import argv

samples=argv[2]
output ={}
output["sample_captures"] ={}
output["diagnosis"] = {}
output["sample_type"] = {}
output["DNASeq"] = {}
output["library"] = {}
output["matched_normal"] = {}
output["matched_rnaseq"] = {}
output["RNASeq"] = {}
output["source"] = {}
#patientIndex=0
#TypeIndex=1
#diagnosisIndex=2
#captureIndex=4
#FCIDIndex=5
#libraryIndex=6
#normRefIndex=7
#rnaRefIndex=8
#source=??
### Column headers in current txt input file
#custom ID
#Type
#diagnosis
#Type of sequencing
#Enrichment step
#Source
#FCID
#Library ID
#Matched normal
#Matched RNA-seq lib
#Case Name
f = open(argv[1], 'r')
for line in f:
	#line = line.rstrip()	
	column = line.split("\t")
	if re.search("custom ID", line):
		patientIndex	=column.index('##custom ID')
		TypeIndex	=column.index('Type')
		diagnosisIndex	=column.index('Diagnosis')
		captureIndex	=column.index('Enrichment step')
		sourceIndex	=column.index('Source')
		FCIDIndex	=column.index('FCID')
		libraryIndex	=column.index('Library ID')
		normRefIndex	=column.index('Matched normal')
		rnaRefIndex	=column.index('Matched RNA-seq lib')
	column[captureIndex] =column[captureIndex].lower()
	# Change it to work with current pipeline
	if column[TypeIndex]  == 'tumor RNA':
		column[TypeIndex] = 'RNASeq'
	elif column[TypeIndex]  == 'tumor DNA':
		column[TypeIndex] = 'Tumor'
	elif column[TypeIndex]  == 'normal DNA' or column[TypeIndex]  == 'blood DNA':
		column[TypeIndex] = 'Normal'

	for sample in samples.split(','):
		if column[patientIndex] == sample:
			output["sample_captures"][column[libraryIndex]]=column[captureIndex]
			output["diagnosis"][column[patientIndex]]=column[diagnosisIndex]
			output["source"][column[libraryIndex]]=column[sourceIndex]
			output["sample_type"][column[libraryIndex]]=column[TypeIndex]
			if column[normRefIndex]:
				output["matched_normal"][column[libraryIndex]]=[column[normRefIndex]]
			if column[rnaRefIndex]:
				output["matched_rnaseq"][column[libraryIndex]] =[column[rnaRefIndex]]
			if not column[FCIDIndex]:
				output["library"][column[libraryIndex]] = [column[libraryIndex]]
			else:
				output["library"][column[libraryIndex]] = [column[libraryIndex]+"_"+column[FCIDIndex]]
			if 'Normal' in column[TypeIndex] or 'Tumor' in column[TypeIndex]:
				if column[patientIndex] not in output["DNASeq"].keys():	
					output["DNASeq"][column[patientIndex]]=[column[libraryIndex]]
				else:
					output["DNASeq"][column[patientIndex]].append(column[libraryIndex])
			else:
				if column[patientIndex] not in output["RNASeq"].keys():
					output["RNASeq"][column[patientIndex]]=[column[libraryIndex]]
				else:
					output["RNASeq"][column[patientIndex]].append(column[libraryIndex])
		

print(json.dumps(output, sort_keys=True, indent=4))
