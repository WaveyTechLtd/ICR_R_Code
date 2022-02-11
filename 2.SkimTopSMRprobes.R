SkimSMRprobes_GliomaGWAS_eQTL <- function(
  # n.cases.outcome, # number of cases in the outcome GWAS
  # n.controls.outcome, # number of controls in the outcome GWAS
  dir.input, # directory containing all GWAS data, in meta-output format, should contain all files to consider and no others
  dir.output,
  nfiles # number of chromosome files to consider
) {
  # required libraries
  library(tidyverse)
  
  for (subtype in subtypes) {
  subtype = "nonGBM"  
    
    # Create an empty tibbles for the final combined data
    SMR.data_full_combined <- tibble()
    
    # First loop through all 22 chromosomes to produce a combined results file
    for (chr in nfiles) {
      
      # Read in correct chromosome file
      message(paste0("Reading in smr results for chr - ", chr))
      filename.input <- file.path(dir.input, paste0("GliomaGWAS_CAGE_SMR_results/",subtype,"/Glioma_",subtype,"_GWAS_CHR_", chr,"_CAGEeQTL_eQTL_SMR.msmr"))
      SMR.data <- read.delim(filename.input, sep="\t", header=T)                        
      
      # rbind SMR.data with combined table
      SMR.data_full_combined <- rbind(SMR.data_full_combined, SMR.data)
    }
    
    write.table(SMR.data_full_combined, 
                file=paste0(dir.output, paste0("GliomaGWAS_CAGE_SMR_results/",subtype,"/Glioma_",subtype,"_GWAS_CAGEeQTL_eQTL_SMR_results_allchr_combined.tsv")),
                quote=F,
                row.names=F,
                sep="\t")
  }
}

SkimSMRprobes_TelomereGWAS_eQTL_20_loci_only <- function(
  # n.cases.outcome, # number of cases in the outcome GWAS
  # n.controls.outcome, # number of controls in the outcome GWAS
  dir.input, # directory containing all GWAS data, in meta-output format, should contain all files to consider and no others
  dir.output,
  nfiles # number of chromosome files to consider
) {
  # required libraries
  library(tidyverse)
  
  # Pull in table with 20 Codd2020 loci
  # This table was manually curated, to include one probe location +/- 1MB of probe per row
  # Previously had list of all 20 probes used, but this lead to duplicates in results table when iterating through rows
  Codd_loci <- read.delim(file = paste0(dir.input, "data/Telomere_GWAS_Codd2020/Codd_loci_for_filtering_eQTL_data.tsv"), sep="\t", header=T)
  
  # list of unique chromosomes for simple filtering down to chromosomes of Codd Loci 
  chr_list <- c(1,3,4,5,6,7,10,11,14,16,19,20)
  
  # Create an empty tibbles for the final combined data
  SMR.data_full_combined <- tibble()
  SMR.data_trimmed_to_probe_locations_only <- tibble()
  
  # First loop through all 22 chromosomes to produce a combined results file
  for (chr in nfiles) {
    
    # Read in correct chromosome file
    message(paste0("Reading in smr results for chr - ", chr))
    filename.input <- file.path(dir.input, paste0("TelomereGWAS_CAGE_eQTL_SMR_results/Telomere_GWAS_CHR_", chr,"_CAGEeQTL_eQTL_SMR.msmr"))
    SMR.data <- read.delim(filename.input, sep="\t", header=T)                        
    
    # rbind SMR.data with combined table
    SMR.data_full_combined <- rbind(SMR.data_full_combined, SMR.data)
  }
  
  write.table(SMR.data_full_combined, 
              file=paste0(dir.output, "TelomereGWAS_CAGE_eQTL_SMR_results/Telomere_GWAS_CAGEL_eQTL_SMR_results_allchr_combined.tsv"),
              quote=F,
              row.names=F,
              sep="\t")
  
  # Loop through 20 Codd loci. filtering down to +/- 1MB of probe locus for each of the relevant chromosomes
  for (row in 1:nrow(Codd_loci)) {
    
    locus <-  Codd_loci[row, "locus"]
    chr <-  Codd_loci[row, "chr"]
    lower_pos <-  Codd_loci[row, "lower_pos"]
    upper_pos <-  Codd_loci[row, "upper_pos"]
    
    # Read in correct chromosome file
    message(paste0("Reading in smr for probe - ", locus))
    filename.input <- file.path(dir.input, paste0("TelomereGWAS_CAGE_eQTL_SMR_results/Telomere_GWAS_CHR_", chr,"_CAGEeQTL_eQTL_SMR.msmr"))
    SMR.data <- read.delim(filename.input, sep="\t", header=T)
    
    # Space for filtering steps ...
    SMR.data %>%
      filter(Probe_bp <= upper_pos) %>%
      filter(Probe_bp >= lower_pos) -> SMR.data.trim
    
    # rbind SMR.data with combined table
    SMR.data_trimmed_to_probe_locations_only <- rbind(SMR.data_trimmed_to_probe_locations_only, SMR.data.trim)
  }
  
  write.table(SMR.data_trimmed_to_probe_locations_only, 
              file=paste0(dir.output, "TelomereGWAS_CAGE_eQTL_SMR_results/Telomere_GWAS_CAGE_eQTL_SMR_20Coddloci_combined.tsv"),
              quote=F,
              row.names=F,
              sep="\t")
}

##################################################
# main
##################################################
if (sys.nframe() == 0L) {
  ### start section to edit
  nfiles = 1:22
  subtypes = c("ALL", "GBM", "nonGBM")
  
  dir.input <- "//rds.icr.ac.uk/DATA/DGE/DUDGE/MOPOPGEN/csaunders/SMR_Telomere_Length_Project/2020_reanalysis/" 
  dir.output <- "//rds.icr.ac.uk/DATA/DGE/DUDGE/MOPOPGEN/csaunders/SMR_Telomere_Length_Project/2020_reanalysis/" 
  ### end section to edit
  
  # reformat the heading of GWAS data
  SkimSMRprobes_TelomereGWAS_eQTL(dir.input, dir.output, nfiles)
}
