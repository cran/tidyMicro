## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("CharlieCarpenter/tidyMicro")
#  library(tidyMicro); library(magrittr)

## ---- include = FALSE---------------------------------------------------------
library(tidyMicro); library(magrittr)

## ----data---------------------------------------------------------------------
## Loading OTU tables
data(phy, package = "tidyMicro")
data(cla, package = "tidyMicro")
data(ord, package = "tidyMicro")
data(fam, package = "tidyMicro")

## Loading meta data to merge with OTU table
data(clin, package = "tidyMicro")

## ---- echo = F----------------------------------------------------------------
cla[1:6, 1:6] %>% knitr::kable()

## ---- eval=FALSE--------------------------------------------------------------
#  
#  ## 1. Single OTU table
#  micro.set <- tidy_micro(otu_tabs = cla,        ## OTU Table
#                          tab_names = "Class",   ## OTU Names (Ranks)
#                          clinical = clin)       ## Clinical Data
#  
#  ## 2. Unnamed List
#  otu_tabs <- list(phy, cla, ord, fam)
#  tab_names <- c("Phylum", "Class", "Order", "Family")
#  
#  micro.set <- tidy_micro(otu_tabs = otu_tabs,   ## OTU Table
#                          tab_names = tab_names, ## OTU Names (Ranks)
#                          clinical = clin)       ## Clinical Data
#  
#  ## 3. Named List
#  otu_tabs <- list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#  
#  micro.set <- tidy_micro(otu_tabs = otu_tabs,  ## OTU Table
#                          clinical = clin)      ## Clinical Data

## ---- include = F-------------------------------------------------------------
otu_tabs <- list(Phylum = phy, Class = cla, Order = ord, Family = fam)

micro.set <- tidy_micro(otu_tabs = otu_tabs,  ## OTU Table
                        clinical = clin)      ## Clinical Data

## ---- echo = F----------------------------------------------------------------
micro.set[1:6, 1:12] %>% knitr::kable()

## -----------------------------------------------------------------------------
micro.set %<>% filter(day == 7)

## ----taxaSummary--------------------------------------------------------------
taxa_summary(micro.set, table = "Phylum") %>% 
  knitr::kable()

## ----PCA, fig.cap = "PCA Plot. Principle components of CLR transformed Family level OTU table.", message=F, warning=F, fig.width=6, fig.height=4, fig.width=6, fig.height=4----
micro.set %>% micro_pca(table = "Family",   ## Taxonomic table of interest
                       grp_var = bpd1, ## A factor variable for colors
                       legend_title = "BPD Severity")      

## ----PCoA, fig.cap = "PCoA Plot. Principle coordinates of Family level OTU table Bray-Curtis dissimilarity.", fig.width=6, fig.height=4----
bray_beta <- micro.set %>%  ## Calculating dissimilarity
  beta_div(table = "Family")
  
micro.set %>% 
  micro_pca(dist = bray_beta, ## Beta diversity
            grp_var = bpd1, ## A factor variable for colors
            legend_title = "BPD Severity")      

## -----------------------------------------------------------------------------
long_micro <- tidy_micro(otu_tabs = otu_tabs,  ## OTU tables
                         clinical = clin)      ## clinical Data

## ----ThreeMode, fig.cap="Three Mode PCA Plot. Principle components controlled for repeaated measures.", fig.width=6, fig.height=4----
long_micro %>% 
  three_mode(table = "Family", group = bpd1, subject = study_id, 
             time_var = day, main = "ThreeMode PCA", 
             subtitle = "3 Time Points", legend_title = "BPD Severity")

## ----pca3dTime, fig.cap="3D Time PCoA Plot. Plotting principle coordinants collapsing over time axis.", fig.width=6, fig.height=4----
long_micro %>%
  pca_3d(table = "Family", time_var = day, subject = study_id, 
         modes = "AC", type = "PCoA")

## ----pca3DSubj, fig.cap="3D Subject PCoA Plot. Plotting principle coordinants collapsing over subject axis.", fig.width=6, fig.height=4----
long_micro %>%
  pca_3d(table = "Family", time_var = day, subject = study_id, 
         modes = "CB", type = "PCoA")

## ----raBars, fig.cap = "Stack Bar Chart. Stacked bar charts of taxa relative abundance by BPD severity.", fig.width=6, fig.height=4----
ra_bars(micro.set,         ## Dataset
        table = "Phylum",  ## Table we want
        bpd1,              ## Variable of interest
        ylab = "% RA", 
        xlab = "BPD", 
        main = "Stacked Bar Charts")

## ----raBarsFilter, fig.cap = "'Taxa-Filtered' Stack Bar Chart. Stacked bar charts of taxa relative abundance using aggregated counts.", fig.width=6, fig.height=4----
ra_bars(micro.set,          ## Dataset
         table = "Phylum",  ## Table we want
         bpd1,              ## Variable of interest
         top_taxa = 3,     
         RA = 0,
         specific_taxa = c("Actinobacteria", "Bacteroidetes"),
         ylab = "% RA", xlab = "BPD", main = "Stacked Bar Charts")

## ----raBarsFactors, fig.cap = "Multiple Factor Stack Bar Chart. Stacked bar charts of taxa relative abundance by multiple factors.", fig.width=6, fig.height=4----

ra_bars(micro.set,           ## Dataset
        table = "Phylum",    ## Table we want
        bpd1, gender,        ## Variables of interest
        top_taxa = 6,
        ylab = "% RA", xlab = "BPD by Sex", 
        main = "Stacked Bar Charts") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----raBarsSubject, fig.cap = "Subject Stack Bar Chart. Subject level stacked bar charts of taxa relative abundance.", fig.width=6, fig.height=4----
ra_bars(micro.set,          ## Dataset
        table = "Phylum",   ## Table we want
        Lib,                ## Variable of interest
        top_taxa = 6,
        ylab = "% RA", xlab = "Library", main = "Stacked Bar Charts") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----raBarsCohort, fig.cap = "Cohort Stack Bar Chart. Overall microbiome community taxa relative abundance.", fig.width=6, fig.height=4----
ra_bars(micro.set,            ## Dataset
        table = "Phylum",     ## Table we want
        top_taxa = 6, 
        ylab = "% RA", main = "Stacked Bar Charts")

## ----raBarsManipulated, fig.cap = "Manipulated Stack Bar Chart. Stacked bar charts of taxa relative abundance by multiple factors.", fig.width=6, fig.height=4----
ra_bars(micro.set,            ## Dataset
        table = "Phylum",     ## Table we want
        top_taxa = 6,
        main = "Manipulated Stacked Bar Charts") +
  
  ## Additional geoms
  theme_dark() + 
  coord_flip() + 
  theme(legend.title = element_text(color = "blue", size = 20),
        legend.text = element_text(color = "red"))

## ----taxaBoxplot, fig.cap="Single Variable Box Plots. Box plots of Staphylococcaceae relative abundance by BPD severity.", fig.width=6, fig.height=4----
staph <- "Firmicutes/Bacilli/Bacillales/Staphylococcaceae"

taxa_boxplot(micro.set,      ## Our dataset
             taxa = staph,   ## Taxa we are interested in
             bpd1,           ## Variable of interest
             xlab = "BPD", 
             ylab = "Relative Abundance", 
             main = "Box Plot") 

## ----taxaBoxplotMultipleFactor, fig.cap="Multi-Variable Box Plots. Box plots of Staphylococcaceae relative abundance by BPD severity and sex.", fig.width=6, fig.height=4----
taxa_boxplot(micro.set,        ## Our dataset
             taxa = staph,     ## Taxa we are interested in
             y = clr,          ## Making Boxplot of CLR
             bpd1, gender,     ## Variables of interest
             ylab = "Staphylococcaceae CLR",
             main = "Box plot", subtitle = "Subtitle") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----corHeatmap, message=F, warning=F, fig.cap= "Lasagna Plot. Spearman orrelations between CLR taxa counts infant weight, gestational age (weeks).", fig.width=6, fig.height=4----
micro.set %>% cor_heatmap(table = "Class", weight, gestational_age)

## ----corRockyMtn, message=F, warning=F, fig.cap="Correlation Rocky Mountain Plot. Rocky mountain plot of Spearman's correlation betwen CLR taxa counts and class level taxa.", fig.width=8, fig.height=6----
micro.set %>% cor_rocky_mtn(table = "Class",
                            weight,
                            cor_label = 0.3)

## ----microAlpha---------------------------------------------------------------
micro_alpha <- alpha_div(micro.set, 
                         table = "Family",  ## Table of interest
                         min_depth = 5000,  ## Requires a Seq Depth of 5000 to be included
                         min_goods = 80)    ## Requires a Good's coverage of %80 

## ----alphaReg, eval = F-------------------------------------------------------
#  micro_alpha %>%
#    micro_alpha_reg(table = "Family", bpd1, gender) %>%
#    knitr::kable()

## ---- echo = F----------------------------------------------------------------
micro_alpha %>% micro_alpha_reg(table = "Family", bpd1, gender) %>% 
  mutate(Beta = round(Beta, 2), std.error = round(std.error, 2), t.stat = round(t.stat,4),
         p.value = round(p.value,4))

## ----betaDiv------------------------------------------------------------------
## Beta Diversity
bray <- beta_div(micro.set, table = "Class", method = "bray")

## ----betaHeatmap, fig.cap="Beta Diversity Heat Map. Heat map of Bray-Curtis beta diversity grouped by BPD severity.", fig.width=6, fig.height=4----
bray %>% 
  beta_heatmap(micro.set, bpd1)

## ----PERMANOVA----------------------------------------------------------------
micro_PERMANOVA(micro.set, ## micro_set to pull covariates from
                bray,      ## Beta diversity matrix (or any distance matrix)
                method = "bray",   ## method used to calculate the beta diversity
                bpd1, mom_ethncty_2)  ## Covariates

## ----taxaFilter---------------------------------------------------------------
## Taxa names "Bacteria" are essentially unclassified
exclude_taxa <- c("Unclassified", "Bacteria")

micro.filt <- micro.set %>% 
  otu_filter(prev_cutoff = 1,               ## Prevalence cutoff
             ra_cutoff = 0.1,               ## Relative abundance cutoff
             exclude_taxa = exclude_taxa)   ## Uninteresting taxa

## ----readInFilter, eval = F---------------------------------------------------
#  ## Named List
#  otu_tabs <- list(Phylum = phy, Class = cla, Ord = ord, Family = fam)
#  
#  tidy.filt <- tidy_micro(otu_tabs = otu_tabs,  ## OTU Table
#                          clinical = clin,      ## Clinical Data
#                          prev_cutoff = 1,      ## Prevalence cutoff
#                          ra_cutoff = 1,        ## Relative abundance cutoff
#                          exclude_taxa = exclude_taxa)   ## Uninteresting taxa

## ----nbMods, message=F, warning=F---------------------------------------------
nb_fam <- micro.filt %>%       ## micro_set
  nb_mods(table = "Family",    ## Rank of taxa we want to model
          bpd1)                ## The covariate in our model

## ----nbInt, eval = F----------------------------------------------------------
#  ## If we wanted the covariates to be bpd1+gender+bpd1*gender we just need input bpd1*gender.
#  
#  nb_int <- micro.filt %>%
#    nb_mods(table = "Class", bpd1*gender)

## ----nbConvergentSummary------------------------------------------------------
nb_fam$Convergent_Summary %>% .[1:6,] %>% knitr::kable()

## ----nbEstimateSummary--------------------------------------------------------
nb_fam$Estimate_Summary %>% .[1:4,] %>% knitr::kable()

## ----bbMods, message=F, warning=F, eval = F-----------------------------------
#  bb_fam <- micro.filt %>%       ## micro_set
#    bb_mods(table = "Phylum",    ## Table we want to model
#            bpd1)                ## The covariate in our model

## ----bbInt, eval = F----------------------------------------------------------
#  ## If we wanted the covariates to be bpd1+gender+bpd1*gender we just need input bpd1*gender.
#  
#  bb_int <- micro.filt %>%
#    bb_mods(table = "Class", bpd1*gender)

## ----nbFamBarCharts, warning=FALSE, fig.cap="Negative Binomial Bar Charts. Stacked bar charts of negative binomial estimated taxa RA by BPD severity.", fig.width=6, fig.height=4----
nb_fam %>% nb_bars(bpd1,                 ## Covariate of interest
                   top_taxa = 5,         ## How many named taxa we want
                   xlab = "Group", 
                   xaxis = c("1","2","3")) ## Labels

## ----manipulatedSBC, warning=FALSE, fig.width=6, fig.height=4-----------------
nb_fam$Model_Coef$Taxa %<>% 
  stringr::str_split("/") %>%  ## Splitting by "/" in taxa names
  lapply(function(x) x[length(x)]) %>% ## selecting piece after the last "/"
  unlist ## Unlisting to put back into data frame

## Reordering to put "Other" at the bottom of the legend
## our functions usually do this automatically, but we'll need to do 
## this externally since we are messing with the output
non.other <- nb_fam$Model_Coef %>% 
  filter(Taxa != "Other") %>% 
  arrange(Taxa) %>% 
  distinct(Taxa) %>%
  pull(Taxa)

nb_fam$Model_Coef$Taxa <- factor(nb_fam$Model_Coef$Taxa,
                                 levels = c(non.other, "Other"))

nb_fam %>% nb_bars(bpd1,                 ## Covariate of interest
                   top_taxa = 5,         ## How many named taxa we want
                   xlab = "BPD Severity", main = "Cleaner Legend") ## Labels

## ----nbRockyMtn, fig.cap="Negative Bionomial Rocky Mountian Plot. Rocky mountain plot created from log FDR adjusted p-values of beta coefficients. Log p-values from beta estimate > 0 are multiplied by -1.", fig.width=10, fig.height=6, warning=F, message=F----
## Order level models
nb_ord <- micro.filt %>%       ## micro_set
  nb_mods(table = "Order",     ## Rank of taxa we want to model
          bpd1)                ## The covariate in our model

nb_ord %>%
  micro_rocky_mtn(bpd1, ## Covariate of interest
                  xlab = "Taxa", main = "Rocky Mountain Plot",
                  subtitle = "Direction of bar indicates direction of relationship", facet_labels = c("Moderate", "Severe"))
  

## ----nbForest, warning = F, message = F, fig.wide = T, fig.cap="Forest Plot. A forest plot showing the beta estimates and 95% confidence intervals of selected covariate from each taxa model.", fig.width=6, fig.height=4----
## Class level models
nb_cla <- micro.filt %>%       ## micro_set
  nb_mods(table = "Class",     ## Rank of taxa we want to model
          bpd1)                ## The covariate in our model

nb_cla %>%
  micro_forest(bpd1, ## Covariate of interest
               main = "Forest Plot for BPD Severity")

## ----microHeatmap, warning=F, message=F, fig.cap="Heat Map of Parameter Estimates. Heat map (or lasagna plot) created from beta estimates of taxa models. '*' indicates FDR adjusted p-values < 0.05.", fig.width=6, fig.height=4----
nb_cla %>% micro_heatmap(top_taxa = 7)

## ----rankSum, message=F, warning=F--------------------------------------------
## Kurskal-Wallis on every taxa
micro.filt %>% 
  micro_rank_sum(table = "Family", grp_var = bpd1) %>% 
  knitr::kable()

## ----rankSumNB, message=F, warning=F------------------------------------------
## Kruskal-Wallis on non-convergent taxa
micro.filt %>% 
  micro_rank_sum(table = "Family", grp_var = bpd1, mod = nb_fam) %>% 
  knitr::kable()

## ----chisq, message=F, warning=F----------------------------------------------
## Chisq on every taxa
micro.filt %>% 
  micro_chisq(table = "Family", grp_var = bpd1, simulate.p.value = T) %>% 
  knitr::kable()

## ----chisqNB, message=F, warning=F--------------------------------------------
## Chisq on non-convergent taxa
micro.filt %>% 
  micro_chisq(table = "Family", grp_var = bpd1, mod = nb_fam) %>% 
  knitr::kable()

## -----------------------------------------------------------------------------
sessionInfo()

