#' @title Create stacked bar charts based on negative binomial model estimates
#' @name nb_bars
#' @description nb_bars takes the output from nb_mods and creates stacked bar charts of the estimated relative abundance for each taxa. The benefit of modeling each taxa before created stacked bar charts is the ability to control for potential confounders. The function will facet wrap interaction terms. Currently, only quant_style "discrete" can be used for an interaction between two quantitative variables
#' @param modsum The output from nb_mods
#' @param ... The covariate you'd like to plot. Can be an interaction term or main effect, but must be in the models created by nb_mods
#' @param range The range you'd like to plot over for a quantitative variable. Will default to the IQR
#' @param quant_style "continuous" will plot over the entire range specified; "discrete" will plot only the endpoints of the range specified. "continuous" by default. This option is ignored without a quantitative variable
#' @param top_taxa Only plot X taxa with the highest relative abundance. The rest will be aggregated into an "Other" category
#' @param RA Only plot taxa with a relative abundance higher than X. The rest will be aggregated into an "Other" category
#' @param specific_taxa Plot this specific taxa even if it doesn't meet the top_taxa or RA requirements
#' @param lines Logical; Add outlines around the different taxa colors in the stacked bar charts
#' @param xaxis Labels for the x-axis ticks. Most useful for categorical variables and defaults to the levels
#' @param main Plot title
#' @param subtitle Subtitle for the plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param facet_labels Labels for the facets created for interaction terms
#' @param facet_layout Rearrange the facets created for interaction terms
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Creating negative binomial models on filtered tidy_micro set
#' nb_fam <- set %>%
#' otu_filter(ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' nb_mods(table = "Family", bpd1)
#'
#' nb_fam %>%
#' nb_bars(bpd1, top_taxa = 9, xlab = "BPD Severity")
#' @export
nb_bars <- function(modsum, ..., range, quant_style = c("continuous", "discrete"),
                    top_taxa = 0, RA = 0, specific_taxa = NULL, lines = TRUE,
                    xaxis, main, subtitle, xlab, ylab, facet_labels = NULL, facet_layout = 1){

  if(!is.double(top_taxa)) stop("top_taxa must be an integer.")
  if(RA < 0 | RA > 100) stop("RA must be between 0 and 100.")
  if(top_taxa > 0 & RA > 0) stop("Can not plot when both top_taxa and RA are greater than 0.")
  if(top_taxa > length(unique(modsum$Convergent_Summary$Taxa))){
    stop("top_taxa must be equal to or less than the number of convergent taxa.")
  }
  if(modsum$Model_Type != "nb_mod") stop("'modsum' should be the output from 'nb_mods'")

  cc <- nb_type(modsum, ...)

  if(length(cc) == 0) stop("Variable specified for bar charts is not in original model.")

  ## Making leaving labels blank optional
  if(missing(main)) main <- NULL ; if(missing(xlab)) xlab <- NULL
  if(missing(subtitle)) subtitle <- NULL
  if(missing(ylab)) ylab <- "Relative Abundance (%)" ; if(missing(facet_labels)) facet_labels <- NULL
  if(missing(quant_style)) quant_style <- "continuous"
  if(facet_layout %nin% c(1,2)) stop("facet_layout must be either 1 or 2")

  if(cc == "categ") {
    if(missing(xaxis)) xaxis <- modsum$Model_Covs[,cov_str(...)] %>% levels

    nb_bars_categ(modsum, ..., cc = cc, top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle)
  } else if(cc == "quant"){

    ## Setting range to 1st and 3rd quartile if not specified
    if(missing(range)) {
      message("Range not specified and will be set to the 1st and 3rd quartile.\n")
      range <- modsum$Model_Covs[,cov_str(...)] %>%
        stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)
    }
    if(!is.numeric(range) | length(range) != 2) stop("range must be a numeric vector of length 2")
    if(missing(xaxis)) xaxis <- as.character(range)

    nb_bars_quant(modsum, ..., cc = cc, quant_style = quant_style, range = range, top_taxa = top_taxa,
                  lines = lines, RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle)
  } else if(cc == "c*c.int"){
    if(missing(xaxis)) xaxis <- modsum$Model_Covs[,cov_str(...)[1]] %>% levels

    nb_bars_ccint(modsum, ..., cc = cc, top_taxa = top_taxa,
                   RA = RA, specific_taxa = specific_taxa, main = main, lines = lines,
                   xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout)

  } else if(cc == "q*c.int"){

    ## Setting range to 1st and 3rd quartile if not specified
    if(missing(range)) {
      message("Range not specified and will be set to the 1st and 3rd quartile.\n")
      range <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select_if(is.numeric) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
        stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)
    }
    if(!is.numeric(range) | length(range) != 2) stop("range must be a numeric vector of length 2")
    if(missing(xaxis)) xaxis <- as.character(range)

    nb_bars_qcint(modsum, ..., cc = cc, range = range, quant_style = quant_style, top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout)
  } else if(cc == "q*q.int"){

    if(missing(range)){
      message("Ranges not specified and will be set to the 1st and 3rd quartile.\n")

      r1 <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select(1) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
        stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)

      r2 <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select(2) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
        stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)

      range <- rbind(r1,r2)
    }

    if(all(dim(range) != 2) | !is.matrix(range)) stop("range must be a 2x2 matrix")

    if(missing(xaxis)){ ## xaxis labels when missing
      if(facet_layout == 1) xaxis <- as.character(round(range[1,], 2))
      if(facet_layout == 2) xaxis <- as.character(round(range[2,], 2))
    }

    if(is.null(facet_labels)){ ## facet labels when missing
      if(facet_layout == 1) facet_labels <- as.character(round(range[2,], 2))
      if(facet_layout == 2) facet_labels <- as.character(round(range[1,], 2))
    }
    names(facet_labels) <- c("Low","High")

    nb_bars_qqint(modsum, ..., cc = cc, range = range, top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout)
  }
}
