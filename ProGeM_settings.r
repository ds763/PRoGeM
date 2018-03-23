
######################################## USER-DEFINED SETTINGS ########################################
## 1. DIRECTORIES AND FILES

# Working directory:
base_dir <- ""

# Output directory: 
output_dir <- ""    # this can be the same as
                                                                        # base_dir

# Directory for GTEx eQTL data set:
eQTLdata_dir <- ""

# File containing a curated list of known metabolic phenotype-related genes:
filename_metabolic_genes <- "annotated_metabolic_genes.RData"		# default provided.

# File containing reference genes:
gene_model_filename <- "GRCh37_genes.RData"							# default provided.

# Files containing GWAS sentinel and proxy variants:
sentinel_filename <- "sentinel_information.txt"						# user provided (in 'base_dir').
proxy_filename <- "proxy_information.txt"							# user provided (in 'base_dir').

# File containing VEP annotation:
VEP_filename <- "VEP_output.txt"									# user provided (in 'base_dir').

#------------------------------------------------------------------------------------------------------
## 2. PARAMETERS FOR TOP-DOWN APPROACH

# Set filter for proxy variants based on r2 values:
filtering_required <- TRUE											# TRUE or FALSE.

# If TRUE, set threshold for r2 values:
r2_thresh <- 0.8													# default is 0.8.

# Genomic interval (in kb) either side of the sentinel SNP:
interval_kb <- 500													# default is 500kb.

#------------------------------------------------------------------------------------------------------
## 3. PARAMETERS FOR BOTTOM-UP APPROACH

# Genomic interval (in kb) of the overhang of the left-most and right-most proxy/sentinel variant at 
# each locus:
LD_region_overhang_kb <- 5											# default is 5kb.

# Number of nearest genes that reside nearest to the sentinel variant:
number_of_nearest <- 3												# default is 3.

# Biotype(s) of candidate genes:
# further information: http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
biotype_of_interest <- "protein_coding"								# default is "protein_coding".

# Column indices in the VEP output that contain the following information:
sentinel_rsID_column <- 
proxy_rsID_column <- 
ensembl_gene_id_column <- 
IMPACT_column <- 
r2_column <- NULL    # OPTIONAL: only required if your VEP output needs to be filtered based on r2_thresh

# Tissue(s) of interest for GTEx eQTL data set:
GTEx_tissues <- dir(eQTLdata_dir)[grep("nominal", dir(eQTLdata_dir))]     # variable contains filenames
                                                                          # for all GTEx tissues (n=44) 
                                                                          # as per below.

# [1] "Adipose_Subcutaneous_Analysis.nominal.filtered.txt"                 
# [2] "Adipose_Visceral_Omentum_Analysis.nominal.filtered.txt"             
# [3] "Adrenal_Gland_Analysis.nominal.filtered.txt"                        
# [4] "Artery_Aorta_Analysis.nominal.filtered.txt"                         
# [5] "Artery_Coronary_Analysis.nominal.filtered.txt"                      
# [6] "Artery_Tibial_Analysis.nominal.filtered.txt"                        
# [7] "Brain_Anterior_cingulate_cortex_BA24_Analysis.nominal.filtered.txt" 
# [8] "Brain_Caudate_basal_ganglia_Analysis.nominal.filtered.txt"          
# [9] "Brain_Cerebellar_Hemisphere_Analysis.nominal.filtered.txt"          
# [10] "Brain_Cerebellum_Analysis.nominal.filtered.txt"                     
# [11] "Brain_Cortex_Analysis.nominal.filtered.txt"                         
# [12] "Brain_Frontal_Cortex_BA9_Analysis.nominal.filtered.txt"             
# [13] "Brain_Hippocampus_Analysis.nominal.filtered.txt"                    
# [14] "Brain_Hypothalamus_Analysis.nominal.filtered.txt"                   
# [15] "Brain_Nucleus_accumbens_basal_ganglia_Analysis.nominal.filtered.txt"
# [16] "Brain_Putamen_basal_ganglia_Analysis.nominal.filtered.txt"          
# [17] "Breast_Mammary_Tissue_Analysis.nominal.filtered.txt"                
# [18] "Cells_EBV-transformed_lymphocytes_Analysis.nominal.filtered.txt"    
# [19] "Cells_Transformed_fibroblasts_Analysis.nominal.filtered.txt"        
# [20] "Colon_Sigmoid_Analysis.nominal.filtered.txt"                        
# [21] "Colon_Transverse_Analysis.nominal.filtered.txt"                     
# [22] "Esophagus_Gastroesophageal_Junction_Analysis.nominal.filtered.txt"  
# [23] "Esophagus_Mucosa_Analysis.nominal.filtered.txt"                     
# [24] "Esophagus_Muscularis_Analysis.nominal.filtered.txt"                 
# [25] "Heart_Atrial_Appendage_Analysis.nominal.filtered.txt"               
# [26] "Heart_Left_Ventricle_Analysis.nominal.filtered.txt"                 
# [27] "Liver_Analysis.nominal.filtered.txt"                                
# [28] "Lung_Analysis.nominal.filtered.txt"                                 
# [29] "Muscle_Skeletal_Analysis.nominal.filtered.txt"                      
# [30] "Nerve_Tibial_Analysis.nominal.filtered.txt"                         
# [31] "Ovary_Analysis.nominal.filtered.txt"                                
# [32] "Pancreas_Analysis.nominal.filtered.txt"                             
# [33] "Pituitary_Analysis.nominal.filtered.txt"                            
# [34] "Prostate_Analysis.nominal.filtered.txt"                             
# [35] "Skin_Not_Sun_Exposed_Suprapubic_Analysis.nominal.filtered.txt"      
# [36] "Skin_Sun_Exposed_Lower_leg_Analysis.nominal.filtered.txt"           
# [37] "Small_Intestine_Terminal_Ileum_Analysis.nominal.filtered.txt"       
# [38] "Spleen_Analysis.nominal.filtered.txt"                               
# [39] "Stomach_Analysis.nominal.filtered.txt"                              
# [40] "Testis_Analysis.nominal.filtered.txt"                               
# [41] "Thyroid_Analysis.nominal.filtered.txt"                              
# [42] "Uterus_Analysis.nominal.filtered.txt"                               
# [43] "Vagina_Analysis.nominal.filtered.txt"                               
# [44] "Whole_Blood_Analysis.nominal.filtered.txt"  

tissues_of_interest <- GTEx_tissues[]     # default is all tissues.
                                          # alternatively the user can select specific tissues by
                                          # providing the appropriate indices in the square brackets.

#------------------------------------------------------------------------------------------------------
## 4. EXECUTE ANNOTATION

source(file = file.path(base_dir, "ProGeM_functions.R"))			# provided.
source(file = file.path(base_dir, "ProGeM_commands.R"))	    # provided.

#######################################################################################################
