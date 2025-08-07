library(TCGAbiolinks)
library(DESeq2)
library(dplyr)

source('self_defined_functions.R')

#PDX
#=====================================================================
load('../data/mutational_landscape_normalizedCount.RData') #loaded as: data.final
load('../data/meta_in_pub.RData')

patients_in_meta=unique(meta_in_pub$patient)
RNASeq_patients=sub('^([^~]+)~.+','\\1',colnames(data.final)[!grepl('organoid|PDC',colnames(data.final))])
patients_in_RNASeq=unique(RNASeq_patients)
patients_in_both=intersect(patients_in_RNASeq,patients_in_meta)
meta_in_pub_RNASeq=meta_in_pub[meta_in_pub$patient %in% patients_in_both,]

expression_no_organoid_no_PDC=data.final[,!grepl('organoid|PDC',colnames(data.final))]



#mutational_landscape_oncotree: 'BLCA','COADREAD','HNSC','MEL','NSCLC','PAAD','SARC'
#Not doing sarcoma.Too heterogenous.
#LUAD, LUSC are the oncotree code for PDX's NSCLC (check sup tab 1)
#SKCM is the oncotree code for PDX's MEL

BLCA_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('BLCA')]]
COAD_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('COAD')]]
READ_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('READ')]]
COADREAD_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('COADREAD')]]
HNSC_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$histology_in_figures %in% c('HNSC')]]
MEL_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('MEL')]]
LUAD_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('LUAD')]]
LUSC_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('LUSC')]]
PAAD_ind=colnames(expression_no_organoid_no_PDC)[RNASeq_patients %in% meta_in_pub_RNASeq$patient[meta_in_pub_RNASeq$oncotree %in% c('PAAD')]]

originator_ind=colnames(expression_no_organoid_no_PDC)[grepl('(?i)ORIGINATOR',colnames(expression_no_organoid_no_PDC))]


BLCA_originator_ind=intersect(originator_ind,BLCA_ind)
BLCA_sample_ind=setdiff(BLCA_ind,originator_ind)
COAD_originator_ind=intersect(originator_ind,COAD_ind)
COAD_sample_ind=setdiff(COAD_ind,originator_ind)
READ_originator_ind=intersect(originator_ind,READ_ind)
READ_sample_ind=setdiff(READ_ind,originator_ind)
COADREAD_originator_ind=intersect(originator_ind,COADREAD_ind)
COADREAD_sample_ind=setdiff(COADREAD_ind,originator_ind)
HNSC_originator_ind=intersect(originator_ind,HNSC_ind)
HNSC_sample_ind=setdiff(HNSC_ind,originator_ind)
MEL_originator_ind=intersect(originator_ind,MEL_ind)
MEL_sample_ind=setdiff(MEL_ind,originator_ind)
LUAD_originator_ind=intersect(originator_ind,LUAD_ind)
LUAD_sample_ind=setdiff(LUAD_ind,originator_ind)
LUSC_originator_ind=intersect(originator_ind,LUSC_ind)
LUSC_sample_ind=setdiff(LUSC_ind,originator_ind)
PAAD_originator_ind=intersect(originator_ind,PAAD_ind)
PAAD_sample_ind=setdiff(PAAD_ind,originator_ind)


expression_no_organoid_no_PDC=expression_no_organoid_no_PDC[,c(BLCA_originator_ind,BLCA_sample_ind,
                                                               COAD_originator_ind,COAD_sample_ind,
                                                               READ_originator_ind,READ_sample_ind,
                                                               COADREAD_originator_ind,COADREAD_sample_ind,
                                                               HNSC_originator_ind,HNSC_sample_ind,
                                                               MEL_originator_ind,MEL_sample_ind,
                                                               LUAD_originator_ind,LUAD_sample_ind,
                                                               LUSC_originator_ind,LUSC_sample_ind,
                                                               PAAD_originator_ind,PAAD_sample_ind)]

save(expression_no_organoid_no_PDC,file='../data/expression_no_organoid_no_PDC.RData')
#=====================================================================




#TCGA
#=====================================================================
query <- GDCquery(project = c("TCGA-BLCA","TCGA-COAD","TCGA-READ","TCGA-HNSC","TCGA-LUAD","TCGA-LUSC","TCGA-PAAD","TCGA-SKCM"), 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)


#stranded_counts=SummarizedExperiment::assay(data, "stranded_second")
unstranded_counts=SummarizedExperiment::assay(data, "unstranded")
rownames(unstranded_counts)=rowData(data)$gene_name #change ensemble id to hugo symbols
unstranded_counts=unstranded_counts[!is.na(rownames(unstranded_counts)) & !duplicated(rownames(unstranded_counts)), ] # Handle duplicates and NAs (common in TCGA data)


counts=unstranded_counts


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData(data),  # Sample metadata
  design = ~1               # Intercept-only design (no groups specified)
)
keep <- rowSums(counts(dds)) >= 10 #this is the same cutoff of our PDX pipeline
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
TCGA_DESeq2_normalized_counts <- counts(dds, normalized = TRUE)

#I only need the primary tumor, not metastatic or solid tissue normal
projects=colData(data)$project_id 
sample_type=colData(data)$sample_type
TCGA_DESeq2_normalized_counts=TCGA_DESeq2_normalized_counts[,sample_type=='Primary Tumor'] 
projects=projects[sample_type=='Primary Tumor']

TCGA_inds=sub('TCGA-','',projects)

TCGA_BLCA_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='BLCA']
TCGA_COAD_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='COAD']
TCGA_READ_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='READ']
TCGA_HNSC_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='HNSC']
TCGA_LUAD_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='LUAD']
TCGA_LUSC_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='LUSC']
TCGA_PAAD_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='PAAD']
TCGA_SKCM_ind=colnames(TCGA_DESeq2_normalized_counts)[TCGA_inds=='SKCM']



save(TCGA_DESeq2_normalized_counts,file='../data/TCGA_DESeq2_normalized_counts.RData')
#=====================================================================



#CCLE
#=====================================================================
#All CCLE files downloaded here are from: https://depmap.org/portal
CCLE_expression=read.csv('./CCLE_data/OmicsExpressionGenesExpectedCountProfile.csv') #This is unstranded data + RSEM. description of this data is: https://depmap.org/portal/data_page/?tab=allData

rownames(CCLE_expression)=CCLE_expression[,1]
CCLE_expression=CCLE_expression[,2:ncol(CCLE_expression)]
gene_symbols=sub('\\.\\.ENSG.+','',colnames(CCLE_expression))
duplicated_gene_symbols=check_replicates(gene_symbols) %>% names
CCLE_expression=CCLE_expression[,!(gene_symbols %in% duplicated_gene_symbols)]
colnames(CCLE_expression)=gene_symbols[!(gene_symbols %in% duplicated_gene_symbols)]
CCLE_expression=t(CCLE_expression)
CCLE_expression=convert_val_type(CCLE_expression,as.integer) #some values are not integers

Model=read.csv('../data/CCLE_Model.csv')
ModelID_OncotreeCode=unique(Model[,c('ModelID','OncotreeCode')])
OmicsProfiles=read.csv('./CCLE_data/OmicsProfiles.csv')
ProfileID_ModelID=unique(OmicsProfiles[,c('ProfileID','ModelID')])
ProfileID_ModelID_OncotreeCode=left_join(ProfileID_ModelID,ModelID_OncotreeCode,by='ModelID')


CCLE_meta=left_join(data.frame(ProfileID=colnames(CCLE_expression)),
                    ProfileID_ModelID_OncotreeCode,
                    by='ProfileID')

used_oncotree=c("BLCA","COAD","READ","HNSC","LUAD","LUSC","PAAD","SKCM")
CCLE_meta=CCLE_meta[CCLE_meta$OncotreeCode %in% used_oncotree,]
CCLE_expression=CCLE_expression[,colnames(CCLE_expression) %in% CCLE_meta$ProfileID]



dds <- DESeqDataSetFromMatrix(
  countData = CCLE_expression,
  colData=CCLE_meta[,'OncotreeCode',drop=F],
  design = ~1               # Intercept-only design (no groups specified)
)
keep <- rowSums(counts(dds)) >= 10 #this is the same cutoff of our PDX pipeline
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
CCLE_DESeq2_normalized_counts <- counts(dds, normalized = TRUE)


CCLE_BLCA_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='BLCA']
CCLE_COAD_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='COAD']
CCLE_READ_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='READ']
CCLE_HNSC_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='HNSC']
CCLE_LUAD_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='LUAD']
CCLE_LUSC_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='LUSC']
CCLE_PAAD_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='PAAD']
CCLE_SKCM_ind=CCLE_meta$ProfileID[CCLE_meta$OncotreeCode=='SKCM']

save(CCLE_DESeq2_normalized_counts,file='../data/CCLE_DESeq2_normalized_counts.RData')
#=====================================================================






#save inds
#=====================================================================

PDX_TCGA_CCLE_inds=list(PDX=list(BLCA_originator_ind=BLCA_originator_ind,
                                 BLCA_sample_ind=BLCA_sample_ind,
                                 COAD_originator_ind=COAD_originator_ind,
                                 COAD_sample_ind=COAD_sample_ind,
                                 READ_originator_ind=READ_originator_ind,
                                 READ_sample_ind=READ_sample_ind,
                                 COADREAD_originator_ind=COADREAD_originator_ind,
                                 COADREAD_sample_ind=COADREAD_sample_ind,
                                 HNSC_originator_ind=HNSC_originator_ind,
                                 HNSC_sample_ind=HNSC_sample_ind,
                                 MEL_originator_ind=MEL_originator_ind,
                                 MEL_sample_ind=MEL_sample_ind,
                                 LUAD_originator_ind=LUAD_originator_ind,
                                 LUAD_sample_ind=LUAD_sample_ind,
                                 LUSC_originator_ind=LUSC_originator_ind,
                                 LUSC_sample_ind=LUSC_sample_ind,
                                 PAAD_originator_ind=PAAD_originator_ind,
                                 PAAD_sample_ind=PAAD_sample_ind
                                 ),
                        TCGA=list(TCGA_BLCA_ind=TCGA_BLCA_ind,
                                  TCGA_COAD_ind=TCGA_COAD_ind,
                                  TCGA_READ_ind=TCGA_READ_ind,
                                  TCGA_HNSC_ind=TCGA_HNSC_ind,
                                  TCGA_LUAD_ind=TCGA_LUAD_ind,
                                  TCGA_LUSC_ind=TCGA_LUSC_ind,
                                  TCGA_PAAD_ind=TCGA_PAAD_ind,
                                  TCGA_SKCM_ind=TCGA_SKCM_ind
                                  ),
                        CCLE=list(CCLE_BLCA_ind=CCLE_BLCA_ind,
                                  CCLE_COAD_ind=CCLE_COAD_ind,
                                  CCLE_READ_ind=CCLE_READ_ind,
                                  CCLE_HNSC_ind=CCLE_HNSC_ind,
                                  CCLE_LUAD_ind=CCLE_LUAD_ind,
                                  CCLE_LUSC_ind=CCLE_LUSC_ind,
                                  CCLE_PAAD_ind=CCLE_PAAD_ind,
                                  CCLE_SKCM_ind=CCLE_SKCM_ind
                                  ))

save(PDX_TCGA_CCLE_inds,file='../data/PDX_TCGA_CCLE_inds.RData')
#=====================================================================