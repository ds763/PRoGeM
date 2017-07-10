---
title: "ProGeM"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```
  
ProGeM (**Pr**ioritisation **o**f candidate causal **Ge**nes at **M**olecular QTLs) is a pipeline for the prioritisation of candidate causal genes - or mediating genes - for molecular quantitative trait loci (QTL) in an automatable fashion.

#### 1) How ProGeM works (in brief)

ProGeM adopts a two-pronged strategy consisting of a *"bottom-up"* and a *"top-down"* component, based on the underlying assumptions that in order for a given gene to mediate a molecular QTL that gene must be: (i) impacted functionally by either the sentinel or a corresponding proxy variant (bottom-up) and, (ii) that gene must then in turn impact the functioning of a molecular mechanism or biological process involved in the phenotype in question (top-down).

**1.1. Bottom-up**. ProGeM will identify bottom-up candidate causal genes in two broad ways:

**Proximal candidates**: 
    (i) ProGeM will identify the nearest genes to each sentinel variant.  The user can specify how many nearest genes they would like to pull out (default is 3), which can also be filtered on the basis of biotype of interest (default is "protein-coding").
    (ii) ProGeM will also identify all genes that overlap the so-called "LD region" surrounding each sentinel variant.  This is defined as the range between the left-most and right-most sentinel/proxy variant plus an overhang.  The user can specify the size of the overhang (default is 5kb). 

**Distal candidates**:
    (iii) ProGeM will search GTEx v6p data to identify cis-eQTL target genes of sentinel and proxy variants.  The user can specify the tissues of interest (default is all GTEx tissues, n=44).
    
ProGeM will also annotate bottom-up candidates with information from a user-supplied Variant Effect Predictor (VEP) output.  The pipeline annotates using the IMPACT column, which categorises variants according to whether their expected impact on gene function is "HIGH", "MODERATE", "LOW", or "MODIFIER".  ProGeM will highlight any bottom-up candidate mediating genes that contain either a sentinel or a proxy variant expected to exert a HIGH or a MODERATE IMPACT.  According the VEP, the IMPACT categories encompass the following variant CONSEQUENCES (http://www.ensembl.org/info/genome/variation/predicted_data.html):

``` r
HIGH_impact <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
                 "stop_gained", "frameshift_variant", "stop_lost", "start_lost", 
                 "transcript_amplification")

MODERATE_impact <- c("inframe_insertion", "inframe_deletion", "missense_variant", 
                     "protein_altering_variant", "regulatory_region_ablation")

LOW_impact <- c("splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant",                    "synonymous_variant")

MODIFIER_impact <- c("coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", 
                     "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", 
                     "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
                     "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", 
                     "TF_binding_site_variant", "regulatory_region_amplification", 
                     "feature_elongation", "regulatory_region_variant", "feature_truncation", 
                     "intergenic_variant")
```

**1.2. Top-down**. ProGeM has been designed to identify top-down candidates by pulling out genes that (i) reside "locally" to a sentinel variant AND (ii) have previously been implicated in the broad molecular phenotype in question.  For the former, the user can specify how "locally" is defined (the default is +/-500kb), and for the latter, we provide an .RData file containing a list of genes involved in metabolic-related molecular phenotypes.  The user can of course prepare a similar .RData file containing genes relating to other molecular phenotypes (i.e., lipidomics) for input to ProGeM.  OUr .RData file containing metabolic-related genes was prepared as follows:

- **(i)** Gene Ontology (GO): all genes associated with a GO term containing the keyword "metabolic"

- **(ii)** Kyoto Encyclopedia of Genes and Genomes (KEGG): all genes annotated to the top-level metabolism pathway (hsa01100)

- **(iii)** Orphanet: all genes annotated to the "rare inborn errors of metabolism" (ORPHA:68367) category

- **(iv)** Reactome: all genes annotated to the top-level metabolism pathway (R-HSA-1430728)

- **(v)** Mouse Genome Informatics (MGI): all human orthologs of mouse genes annotated to the Mammalian Phenotype term "abnormal metabolism" (MP:0005266)

**1.3. Cooccurring candidates**. ProGeM will also prioritise genes that appear as both bottom-up AND top-down candidates for a given sentinel variant.  These genes are so-called "concurrent" genes and should be considered primary causal candidates based on the fact that they satisfy both of the above-mentioned underlying assumptions of a causal gene.

#### 2) Before you begin

**2.1.** Download the necessary files to run ProGeM and place them ALL locally into a directory of your choosing.   Your directry should contain the following files:

- **"ProGeM_settings.R"** (a series of ProGeM settings that require user input)
- **"ProGeM_functions.R"** (a series of required functions)
- **"ProGeM_commands.R"** (a series of execution commands)
- **"GRCh37_genes.RData"** (a GRCh37 gene model based on a GTF file ["Homo_sapiens.GRCh37.82.gtf"] downloaded from ensembl)
- **"annotated_metabolic_genes.RData"** (OPTIONAL: a list of GRanges objects containing so-called metabolic-related genes according to five open source databases - see further details below)

**2.2.** Download a zipped folder containing cis-eQTL data (v6p) from the Genotype-Tissue Expression (GTEx) project: https://www.gtexportal.org/home/.  At this site click on the drop-down menu "Datasets" and then click on "Download" - you'll need to register for an account with GTEx first.  The name of the folder you need to download is "GTEx_Analysis_v6p_eQTL.tar" (728M), which contains eGene and significant SNP-gene associations based on permutations.  Unzip this folder and place ALL files locally in a directory of your choosing (preferably not the one containing the ProGeMM_*.R and .RData files from above).

**2.3.** Finally, ensure you have downloaded the following R packages and dependencies:

- **GenomicAlignments** (https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
- **GenomicFeatures** (http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
- **GenomicRanges** (https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- **biomaRt** (https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- **data.table** (https://cran.r-project.org/web/packages/data.table/data.table.pdf)
    
``` r
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments", "GenomicFeatures", "GenomicRanges", "biomaRt")
install.packages("data.table")
```


#### 3) ProGeM input files

The user must prepare and save three files to the same directory already containing the ProGeM_*.R and .RData files.  These three files must contain and be formatted as follows:

**3.1.** A tab-separated .txt file containing rsIDs, chromosomes, and GRCh37 coordinates (both start and end) of your sentinel variants of interest across four columns with the below column names.  In cases where there is no rsID for a sentinel variant then the notation "chr:start" should be used (i.e., 1:11856378).
    
rsID | CHR | START | END 
--- | --- | --- | ---
rs1801133 | 1 | 11856378 | 11856378
rs532545 | 1 | 20915172 | 20915172
rs1697421 | 1 | 21823292 | 21823292

**3.2.** A tab-separated .txt file containing rsIDs, chromosomes, and GRCh37 coordinates (both start and end) of proxies for your sentinel variants; along with the corresponding sentinel rsIDs and r2 value across six columns in total with the below column names.  Your proxies can be derived either directly from your sample or, for example, from 1000 genomes data.  If proxies have not been pre-filtered prior to input, then ProGeMM will filter based on a user-defined r2 threshold (default r2 is 0.8) if desired.

PROXY_rsID | PROXY_CHR | PROXY_START | PROXY_END | LEAD_rsID | r2
--- | --- | --- | --- | --- | ---
rs602950 | 1 | 20915531 | 20915531 |  rs532545 | 0.991
rs1253904 | 1 | 20913519 | 20913519 | rs532545 | 0.963
rs589942 | 1 | 20916080 | 20916080 | rs532545 | 0.931

**3.3.** The output (as a .txt file) after running the Variant Effect Predictor (VEP: http://grch37.ensembl.org/Homo_sapiens/Tools/VEP) using both sentinels and proxies as input.  The user is free to decide upon the parameters used to run the VEP, though the output must contain a column with (i) sentinel rsIDs as well as a column with (ii) proxy rsIDs - in the "ProGeM_settings.R" files the user is asked to indicate the column indices of their VEP output that contains this information.  For VEP records/rows corresponding directly to a sentinel variant, then the sentinel rsID should be entered in both the sentinel rsID and the proxy rsID column.  Further, in the VEP output there should be a separate record/row for each sentinel-proxy combination, meaning that if a variant is a proxy for more than one sentinel variant in your dataset, then there should be more than one VEP record/row for that proxy variant.  


#### 4) Running ProGeM

After downloading, preparing, and making local copies of the necessary files (see above), the user should specify their preferred parameters in the "ProGeM_settings.R" file where indicated.  This file can then be saved and sourced, and ProGeMM will then run in its entirety - simple as that.  In our tests, ProGeM will fully process 213 sentinel variants in ~12-15 minutes on a standard desktop computer (4GB RAM, dual core i3 processor at 3.40GHz) running Window's 7.


#### 5) ProGeM output files

After successful completion, ProGeM will write multiple data tables to file within a user-specified output directory.  The three main summary files are:

**5.1. "OUTPUT_bottom_up_summary.txt"** 

LEAD_rsID	| ensembl_id	| hgnc_symbol |	gene_biotype |	nearest |	second_nearest |	third_nearest |	LD_overlapping |	lead_eQTL |	proxy_eQTL |	lead_IMPACT |	proxy_IMPACT
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
rs1801133 |	ENSG00000177000 |	MTHFR |	protein_coding |	1 |	0 |	0 |	1 |	1 |	0 |	1 |	0
rs1801133 |	ENSG00000215910 |	C1orf167 |	protein_coding |	0 |	1 |	0 |	0 |	0 |	0 |	0 |	0
rs1801133 |	ENSG00000011021 |	CLCN6 |	protein_coding |	0 |	0 |	1 |	0 |	0 |	0 |	0	| 0

Each row provides annotation of a *"bottom-up"* candidate mediating gene, the corresponding sentinel rsID is indicated in **column 1**; **columns 2-4** provide basic annotation of candidate mediators; **columns 5-10** are boolean variables, with 0=FALSE and 1=TRUE; **columns 11-12** are nominal variables referring to VEP IMPACT annotation with 0=MODIFIER or LOW IMPACT, 1=MODERATE IMPACT, 2=HIGH IMPACT.

**5.2. "OUTPUT_top_down_scores.txt"**

LEAD_rsID |	ensembl_id |	hgnc_symbol |	gene_biotype |	GO |	KEGG |	MGI |	orphanet |	reactome |	score
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
rs1801133 |	ENSG00000177000 |	MTHFR |	protein_coding |	1 |	1 |	0 |	1 |	1 |	4
rs1801133 |	ENSG00000175206 |	NPPA |	protein_coding |	1 |	0 |	1 |	0 |	0 |	2
rs1801133 |	ENSG00000204624 |	PTCHD2 |	protein_coding |	1 |	0 |	0 |	0 |	0 |	1

Each row provides annotation of a *"top-down"* candidate mediating gene based solely on top-down evidence, the corresponding sentinel rsID is indicated in **column 1**; **columns 2-4** provide basic annotation of candidate mediators; **columns 5-9** are boolean variables indicating whether the candidate mediator is metabolic-related according to the database in question, with 0=FALSE and 1=TRUE; **column 10** represents a simple score corresponding to the sum of columns 5-9 for each row; i.e., how many databases [out of 5] highlighted a given *top-down* candidate mediator as being involved in metabolism.

**5.3. "OUTPUT_cooccurring_candidates.txt"**

LEAD_rsID |	ensembl_id |	hgnc_symbol |	gene_biotype
--- | --- | --- | ---
rs1801133 |	ENSG00000177000 |	MTHFR |	protein_coding
rs532545 |	ENSG00000158825 |	CDA |	protein_coding
rs532545 |	ENSG00000158828 |	PINK1 |	protein_coding

Each row corresponds to a candidate mediating gene with both *bottom-up* and *top-down* evidence, the corresponding sentinel rsID is indicated in **column 1**; **columns 2-4** provide basic annotation of so-called *"cooccurring"* candidates.

**5.4. Supplementary output files**

In addition to the above three summary output files, ProGeM also writes additional data files to enable the user to probe for more in-depth information should they wish:

**- "OUTPUT_lead_cis_eQTLs.txt":** contains all cis-eQTL targets for all sentinel variants across each of the user-specified tissues from GTEx.

**- "OUTPUT_proxy_cis_eQTLs.txt":** contains all cis-eQTL targets for all proxy variants across each of the user-specified tissues from GTEx.

**- "OUTPUT_VEP_IMPACT_annotations.txt":** contains records from the user-provided VEP output predicted to have either a MODERATE or HIGH IMPACT on gene functioning.

**- "OUTPUT_nearest_genes.txt":** contains basic annotation of genes nearest to each sentinel variant, filtered on the basis of user-specified biotype(s) of interest.  This file contains the distance (bp) between each TSS and the corresponding sentinel variants.  Further, this file will be particularly useful if the user specified to pull out more than the three nearest genes, as the bottom-up summary will only annotate the first three.
