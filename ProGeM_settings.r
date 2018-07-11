
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
                                                                          # for all GTEx v7 tissues (n=44) 
                                                                          # as per below.

# [1] "Adipose_Subcutaneous.signifpairs.txt"                 
# [2] "Adipose_Visceral_Omentum.signifpairs.txt"             
# [3] "Adrenal_Gland.signifpairs.txt"                        
# [4] "Artery_Aorta.signifpairs.txt"                         
# [5] "Artery_Coronary.signifpairs.txt"                      
# [6] "Artery_Tibial.signifpairs.txt"                        
# [7] "Brain_Amygdala.signifpairs.txt"                       
# [8] "Brain_Anterior_cingulate_cortex_BA24.signifpairs.txt" 
# [9] "Brain_Caudate_basal_ganglia.signifpairs.txt"          
# [10] "Brain_Cerebellar_Hemisphere.signifpairs.txt"          
# [11] "Brain_Cerebellum.signifpairs.txt"                     
# [12] "Brain_Cortex.signifpairs.txt"                         
# [13] "Brain_Frontal_Cortex_BA9.signifpairs.txt"             
# [14] "Brain_Hippocampus.signifpairs.txt"                    
# [15] "Brain_Hypothalamus.signifpairs.txt"                   
# [16] "Brain_Nucleus_accumbens_basal_ganglia.signifpairs.txt"
# [17] "Brain_Putamen_basal_ganglia.signifpairs.txt"          
# [18] "Brain_Spinal_cord_cervical_c-1.signifpairs.txt"       
# [19] "Brain_Substantia_nigra.signifpairs.txt"               
# [20] "Breast_Mammary_Tissue.signifpairs.txt"                
# [21] "Cells_EBV-transformed_lymphocytes.signifpairs.txt"    
# [22] "Cells_Transformed_fibroblasts.signifpairs.txt"        
# [23] "Colon_Sigmoid.signifpairs.txt"                        
# [24] "Colon_Transverse.signifpairs.txt"                     
# [25] "Esophagus_Gastroesophageal_Junction.signifpairs.txt"  
# [26] "Esophagus_Mucosa.signifpairs.txt"                     
# [27] "Esophagus_Muscularis.signifpairs.txt"                 
# [28] "Heart_Atrial_Appendage.signifpairs.txt"               
# [29] "Heart_Left_Ventricle.signifpairs.txt"                 
# [30] "Liver.signifpairs.txt"                                
# [31] "Lung.signifpairs.txt"                                 
# [32] "Minor_Salivary_Gland.signifpairs.txt"                 
# [33] "Muscle_Skeletal.signifpairs.txt"                      
# [34] "Nerve_Tibial.signifpairs.txt"                         
# [35] "Ovary.signifpairs.txt"                                
# [36] "Pancreas.signifpairs.txt"                             
# [37] "Pituitary.signifpairs.txt"                            
# [38] "Prostate.signifpairs.txt"                             
# [39] "Skin_Not_Sun_Exposed_Suprapubic.signifpairs.txt"      
# [40] "Skin_Sun_Exposed_Lower_leg.signifpairs.txt"           
# [41] "Small_Intestine_Terminal_Ileum.signifpairs.txt"       
# [42] "Spleen.signifpairs.txt"                               
# [43] "Stomach.signifpairs.txt"                              
# [44] "Testis.signifpairs.txt"                               
# [45] "Thyroid.signifpairs.txt"                              
# [46] "Uterus.signifpairs.txt"                               
# [47] "Vagina.signifpairs.txt"                               
# [48] "Whole_Blood.signifpairs.txt" 

tissues_of_interest <- GTEx_tissues[]     # default is all tissues.
                                          # alternatively the user can select specific tissues by
                                          # providing the appropriate indices in the square brackets.

#------------------------------------------------------------------------------------------------------
## 4. EXECUTE ANNOTATION

source(file = file.path(base_dir, "ProGeM_functions.R"))			# provided.
source(file = file.path(base_dir, "ProGeM_commands.R"))	    # provided.

#######################################################################################################
