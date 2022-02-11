# ICR_R_Code

Two simple code examples using R, written as part of my work at the Institute of Cancer Research, London.

(1) 2.SkimTopSMRprobes.R
Publication - Saunders et al 2022 NeuroOncology, DOI = https://doi.org/10.1093/neuonc/noab208

Part of a project investigating telomere length and glioma risk. Takes Summary Mendelian Randomisation (SMR) results data, loops through a list of significant hits, filtering  all the surrounding eQTL probes (both significant and non-singificant) in the initial dataset, ready for later analysis.

(2) plotting_regression_plots.R
Publication - Saunders et al 2020 British Journal of Cancer, DOI = https://doi.org/10.1038/s41416-020-01083-1

Simple tidyverse/ggplot2 script to plot linear regression lines for Mendelian Randomisation results, again looking at the relationship between telomere length and glioma risk, using a variety of different models. 
