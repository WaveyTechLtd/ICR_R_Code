plot_regression_curves <- function(LD_clumping_filepath,
                                   LD_clumping_filename) {
  
  # First pull in the results data with all the different gradients of regression lines
  read_delim(paste0(LD_clumping_filepath,"/",LD_clumping_results), delim="\t") %>%
    select(id.exposure, method, b) -> LD_results
  
  # Pull in the LD_assoc stats with the individual beta's for exposures and GWASs
  read_delim(paste0(LD_clumping_filepath,"/",LD_clumping_assocstats), delim="\t") %>%
    select(Exposure, SNP, `Beta exposure`,`SE exposure`, `Beta GWAS`, `SE GWAS`) %>%
    rename(beta_exp=`Beta exposure`, se_exp=`SE exposure`, beta_gwas=`Beta GWAS`, se_gwas=`SE GWAS`) -> LD_assoc_stats
  
  # Filter down to exposure of interest only
  LD_results %>% filter(id.exposure == exposure) -> LD_results
  LD_assoc_stats %>% filter(Exposure == exposure) -> LD_assoc_stats
  
  # Fake data for geom_plot
  tibble("x" = seq(-1,1,0.01), "y" = seq(-1,1,0.01)) %>%
    rowwise() %>%
    mutate("IVW_RE" = unlist(x*(LD_results[1,"b"])),
           "IVW_FE" = unlist(x*(LD_results[2,"b"])),
           "Maximum_likelihood" = unlist(x*(LD_results[3,"b"])),
           "Weighted_Median"= unlist(x*(LD_results[5,"b"])),
           "Weighted_Mode" = unlist(x*(LD_results[7,"b"])),
           "MR_Egger" = unlist(x*(LD_results[8,"b"]))) %>%
    # Needs a filtering step so that 
    filter(x < 0.61, x > -0.42) -> regression_data

LD_assoc_stats %>%
  ggplot() + 
  geom_point(aes(x=beta_exp, y=beta_gwas))+
  geom_errorbar(aes(x = beta_exp, ymin = beta_gwas-se_gwas, ymax = beta_gwas+se_gwas), colour="gray10")+
  geom_errorbarh(aes(y = beta_gwas, xmin = beta_exp-se_exp, xmax = beta_exp+se_exp), colour="gray10")+
  
  # Plot the Regression Lines
  geom_line(aes(y=regression_data$IVW_RE, x=regression_data$x), colour="black", size=1.2)+
  geom_line(aes(y=regression_data$IVW_FE, x=regression_data$x), colour="gray", size=1.2)+
  geom_line(aes(y=regression_data$Maximum_likelihood, x=regression_data$x), colour="red", size=1.2)+
  geom_line(aes(y=regression_data$Weighted_Median, x=regression_data$x), colour="blue", size=1.2)+
  geom_line(aes(y=regression_data$Weighted_Mode, x=regression_data$x), colour="purple", size=1.2)+
  geom_line(aes(y=regression_data$MR_Egger, x=regression_data$x), colour="green", size=1.2)+

  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.2, 0.2)) +
  labs(x="Effect Size Exposure", y="Effect Size GWAS")+
  ggtitle("Scatter plot showing GLGC LDL effect \n sizes in both exposure and GWAS")

}

##################################################
# main
##################################################
if (sys.nframe() == 0L) {
  
  #Import Libraries
  library(tidyverse)
  
  ### start section to edit
  exposure="GLGC_LDL" # "GLGC_HDL" "GLGC_HDL" or "GLGC_TC"
  LD_clumping_filepath = "//rds.icr.ac.uk/DATA/DGE/DUDGE/MOPOPGEN/csaunders/CRC_statins_head_and_neck_SNPS/Final_LD_Clumping_results"
  LD_clumping_results = "mr_analyses.outcome_gwas.tsv"
  LD_clumping_assocstats = "association_statistics_formatted.outcome_gwas.tsv"
  
  ### end section to edit   
  
  plot_power_curves(exposure,
                    LD_clumping_filepath,
                    correl_correct_filepath,
                    correl_correct_filepath,
                    correl_correct_filepath,
                    n.cases.disease, 
                    n.controls.disease)
}
