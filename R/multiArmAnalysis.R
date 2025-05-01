#' Bayesian or Frequentist Analysis with Forest Plot Comparison for Multi-Arm Trial Designs
#'
#' This function fits Bayesian or frequentist and producing a forest plot across multiple
#' intervention groups for cluster randomized trials (CRT), multisite trials (MST) or simple randomized trials (SRT).
#'
#' @param method The model fitting method. Should be specified as a character string. Choices are:
#'   - `"crtBayes"`: Bayesian analysis of cluster randomised trials using vague priors.
#'   - `"crtFREQ"`: Analysis of cluster randomised trials using a multilevel model under a frequentist setting.
#'   - `"mstBayes"`: Bayesian analysis of multisite randomised trials using vague priors.
#'   - `"mstFREQ"`: Analysis of multisite randomised trials using a multilevel model under a frequentist setting.
#'   - `"srtBayes"`: Bayesian analysis of simple randomised trials using vague priors.
#'   - `"srtFREQ"`: Analysis of simple randomised trials under a frequentist setting.
#' @param data A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention.
#' @param outcome The name of the outcome (post-test) variable. 
#' @param interventions A string specifying the intervention variable.
#' @param Random The name of the clustering variable (e.g., schools or sites) for CRT and MST designs.
#' @param Nsim Number of MCMC iterations to be performed for Bayesian analysis. A minimum of 10,000 is recommended to ensure convergence.
#' @param Threshold The effect size threshold for posterior computation for Bayesian analysis (default = 0.05).
#' @param FREQoption The option for frequentist methods. Choices are "Default", "Permutation", or "Bootstrap". 
#' @param nPerm The number of permutations required to generate a permutated p-value.
#' @param nBoot The number of bootstraps required to generate bootstrap confidence intervals.
#' @param bootType method of bootstrapping including case re-sampling at student level "case(1)",case re-sampling at school level "case(2)", case re-sampling at both levels "case(1,2)" and residual bootstrapping using "residual". If not provided, default will be case re-sampling at student level.
#' @param covariates Additional covariates include in the model. It should be specified as a character vector.
#' @param intlabels Optional custom intervention labels for the plot.
#' @param intcolors Optional intervention colors for the plot. 
#' @param maintitle main title for the plot. 
#' @param xlabel Label for the x-axis. 
#' @param ylabel Label for the y-axis.
#' @param vlinecolor Color of the vertical reference line (default = "black"). 
#'
#' @return A `ggplot` object showing intervention effect sizes and their confidence intervals.
#'
#' @details
#' This function loops through each intervention, fits the requested statistical model, 
#' stores the results, and forest plot visualization for easy comparison. It allows flexible customization for plotting aesthetics.
#'
#' @seealso Functions from the \pkg{eefAnalytics} package:
#'   \code{\link[eefAnalytics]{crtBayes}}, \code{\link[eefAnalytics]{crtFREQ}}, \code{\link[eefAnalytics]{mstBayes}}, \code{\link[eefAnalytics]{mstFREQ}}, \code{\link[eefAnalytics]{srtBayes}}, \code{\link[eefAnalytics]{srtFREQ}}
#'
#' @examples
#' \donttest{
#' ### Bayesian analysis of cluster randomised trials ###
#' data(crt4armSimData)
#' multiArmAnalysis(method = "crtBayes", data = crt4armSimData, outcome = "posttest", 
#' interventions = "interventions", Random = "schools", Nsim = 10000, Threshold = 0.05, 
#' covariates = c("pretest"), intlabels = c("Intervention A", "Intervention B", "Intervention C"),
#  intcolors = c("Intervention A" = "blue", "Intervention B" = "green", "Intervention C" = "red"),
#' maintitle = "Forest plot of comparison of effect sizes", xlabel = "Hedges'g", 
#' ylabel = "Interventions", vlinecolor = "black")
#'
#' ###MLM analysis of multisite trials with residual bootstrap confidence intervals ###
#' data(mst4armSimData)
#' multiArmAnalysis(method = "mstFREQ", data = mst4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", nBoot = 1000, bootType="residual",
#' covariates = c("pretest"), intlabels = c("Intervention A", "Intervention B", "Intervention C"),
#' intcolors = c("Intervention A" = "blue", "Intervention B" = "green", "Intervention C" = "red"),
#' maintitle = "Forest plot of comparison of effect sizes ", xlabel = "Hedges'g",
#' ylabel = "Interventions", vlinecolor = "black")
#'
#' ###MLM analysis of multisite trials with permutation p-value###
#' data(mst4armSimData)
#' multiArmAnalysis(method = "mstFREQ", data = mst4armSimData, outcome = "posttest", 
#' interventions = "interventions", Random = "schools", nPerm = 1000, covariates = c("pretest"), 
#' intlabels = c("Intervention A", "Intervention B", "Intervention C"),
#' intcolors = c("Intervention A" = "blue", "Intervention B" = "green", "Intervention C" = "red"), 
#' maintitle = "Forest plot of comparison of effect sizes ",
#' xlabel = "Hedges'g", ylabel = "Interventions", vlinecolor = "black")
#'
#'
#' ###Bayesian analysis of simple randomised trials###
#' data(srt4armSimData)
#' multiArmAnalysis(method = "srtBayes", data = srt4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", Nsim = 10000, Threshold = 0.05, 
#' covariates = c("pretest"), intlabels = c("Int A", "Int B", "Int C"),
#' intcolors = c("Int A" = "#1F77B4", "Int B" = "#2CA02C", "Int C" = "#D62728"),
#' maintitle = "Forest plot of comparison of effect sizes ", xlabel = "Hedges'g", 
#' ylabel = "Interventions", vlinecolor = "black")
#'}
#'
#' @import ggplot2
#' @importFrom eefAnalytics crtBayes crtFREQ mstBayes mstFREQ srtBayes srtFREQ
#'
#' @export
multiArmAnalysis <- function(
    method = "crtBayes", 
    data, 
    outcome = "posttest", 
    interventions = "interventions", 
    Random = "schools", 
    Nsim = 10000, 
    Threshold = 0.05, 
    FREQoption = "Default", 
    nPerm = NULL, 
    nBoot = NULL, 
    bootType = NULL, 
    covariates = NULL, 
    maintitle = NULL,
    xlabel = NULL,
    ylabel = NULL,
    vlinecolor = "black",
    intlabels = NULL,
    intcolors = NULL
) {
  
  # Construct formula
  if (is.null(covariates) || length(covariates) == 0) {
    formula_str <- paste(outcome, "~", interventions)
  } else {
    formula_str <- paste(outcome, "~", interventions, "+", paste(covariates, collapse = " + "))
  }
  
  # Get unique interventions excluding control (0)
  intervention_col <- sort(unique(data[[interventions]][data[[interventions]] != 0]))
  output <- list()
  
  for (i in seq_along(intervention_col)) {
    intervention <- intervention_col[i]
    
    # Filter data for control vs current intervention
    intervention_data <- data[data[[interventions]] %in% c(intervention, 0), ]
    
    # Binary indicator: 1 for intervention, 0 for control
    intervention_data[[interventions]] <- ifelse(intervention_data[[interventions]] == intervention, 1, 0)
    
    # Update formula in case of different covariates
    formula_str <- paste(outcome, "~", interventions)
    if (!is.null(covariates) && length(covariates) > 0) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    
    # Fit the model
    if (method == "crtBayes") {
      output[[i]] <- crtBayes(
        as.formula(formula_str),
        random = Random,
        intervention = interventions,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "crtFREQ") {
      if (FREQoption == "Default") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data)
      } else if (FREQoption == "Permutation") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data, nPerm = nPerm)
      } else if (FREQoption == "Bootstrap") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data, nBoot = nBoot, type = bootType)
      }
      
    } else if (method == "mstBayes") {
      output[[i]] <- mstBayes(
        as.formula(formula_str),
        random = Random,
        intervention = interventions,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "mstFREQ") {
      if (FREQoption == "Default") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data)
      } else if (FREQoption == "Permutation") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data, nPerm = nPerm)
      } else if (FREQoption == "Bootstrap") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = interventions, data = intervention_data, nBoot = nBoot, type = bootType)
      }
      
    } else if (method == "srtBayes") {
      output[[i]] <- srtBayes(
        as.formula(formula_str),
        intervention = interventions,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "srtFREQ") {
      if (FREQoption == "Default") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = interventions, data = intervention_data)
      } else if (FREQoption == "Permutation") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = interventions, data = intervention_data, nPerm = nPerm)
      } else if (FREQoption == "Bootstrap") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = interventions, data = intervention_data, nBoot = nBoot)
      }
    }
  }
  
  # Example: intervention_col = c(1, 2, 3)
  
  # Handle defaults for intlabels and intcolors
  if (is.null(intlabels) || length(intlabels) != length(intervention_col)) {
    intlabels <- paste("Intervention", intervention_col)
  }
  
  if (is.null(intcolors) || length(intcolors) != length(intervention_col)) {
    default_colors <- c("blue", "green", "red", "purple", "orange", "brown", "pink")
    intcolors <- rep(default_colors, length.out = length(intervention_col))
  }
  

  # Generate final comparison plot
  result <- forestPlotMultiArms(
    output,
    modelNames = paste("Intervention", 1:length(intervention_col)),
    group = 1,
    intlabels = intlabels,
    intcolors = intcolors,
    maintitle = maintitle,
    xlabel = xlabel,
    ylabel = ylabel,
    vlinecolor = vlinecolor
  )
  
  return(result)
}
