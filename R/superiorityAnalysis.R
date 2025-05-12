#' Superiority Analysis Across Interventions for CRT, MST, or SRT Designs
#'
#' This function performs a Bayesian superiority analysis, comparing each intervention against a reference
#' intervention, across cluster randomized trials (CRT), multisite trials (MST) or simple randomized trials (SRT).
#'
#' @param method The trial design type: "crt", "mst", or "srt".
#' @param data A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention.
#' @param outcome The name of the outcome (post-test) variable. 
#' @param interventions A string specifying the intervention variable.
#' @param Random The name of the clustering variable (e.g., schools or sites) for CRT and MST designs.
#' @param Nsim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param refintervention The value of the intervention used as the reference group (default = 1).
#' @param Threshold The effect size threshold for posterior computation (default = 0.05).
#' @param SupThreshold The minimum posterior probability threshold to declare superiority (default = 0.8).
#' @param continuous_covariates A character vector specifying the names of continuous covariates.
#' @param categorical_covariates A character vector specifying the names of categorical covariates (converted to factors).
#'
#' @return A `data.frame` with columns:
#' \itemize{
#'   \item `Intervention`: Intervention group identifier.
#'   \item `ProbES`: Posterior probability of superiority over the reference intervention.
#'   \item `Superiority`: Label indicating `"Superior"`, `"Not Superior"`, or `"Reference"`.
#' }
#'
#' @details
#' The effect size is estimated against a reference intervention, which by default is intervention 1 but can be reassigned to any other intervention, including the control (refintervention = 0).
#' @seealso  \code{\link[eefAnalytics]{crtBayes}}, \code{\link[eefAnalytics]{mstBayes}}, \code{\link[eefAnalytics]{srtBayes}} functions from the \pkg{eefAnalytics} package
#'  
#' @examples
#' \donttest{
#' ###Superiority analysis of cluster randomized trial###
#' data(crt4armSimData)
#' superiorityAnalysis(method = "crt", data = crt4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", Nsim = 10000, refintervention = 2,
#' Threshold = 0.05, SupThreshold = 0.8,continuous_covariates = c("pretest"),
#' categorical_covariates = c("gender", "ethnicity"))
#' 
#' ###Superiority analysis of multisite trial###
#' data(mst4armSimData)
#' superiorityAnalysis(method = "mst", data = mst4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", Nsim = 10000, refintervention = 2,
#' Threshold = 0.05, SupThreshold = 0.8,continuous_covariates = c("pretest"),
#' categorical_covariates = c("gender", "ethnicity"))
#' 
#' ###Superiority analysis of simple randomized trial###
#' data(srt4armSimData)
#' superiorityAnalysis(method = "srt", data = srt4armSimData, outcome = "posttest",
#' interventions = "interventions", Nsim = 10000, refintervention = 2,
#' Threshold = 0.05, SupThreshold = 0.8,continuous_covariates = c("pretest"),
#' categorical_covariates = c("gender", "ethnicity"))
#' 
#'}
#'
#' @importFrom eefAnalytics crtBayes mstBayes srtBayes
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
superiorityAnalysis <- function(method = c("crt", "mst", "srt"),
                                data,
                                outcome = "posttest",
                                interventions = "interventions",
                                Random = "schools",
                                Nsim = 10000,
                                Threshold = 0.05,
                                refintervention = 1,
                                SupThreshold = 0.8,
                                continuous_covariates = NULL, categorical_covariates = NULL
) {
  method <- match.arg(method)
  
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  covariates <- c(continuous_covariates, categorical_covariates)
  
  
  # Get unique intervention_col excluding the control group (0)
  intervention_col <- sort(unique(data[[interventions]]))
  intervention_col <- intervention_col[intervention_col != 0 & intervention_col != refintervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(intervention_col))
  sup_decisions <- numeric(length(intervention_col))
  
  for (i in seq_along(intervention_col)) {
    intervention <- intervention_col[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- data %>%
      filter(.data[[interventions]] == intervention | .data[[interventions]] == refintervention)
    intervention_data[[interventions]] <- ifelse(intervention_data[[interventions]] == intervention, 1, 0)
    
    # Ensure factors in subset too
    if (!is.null(categorical_covariates)) {
      for (cat in categorical_covariates) {
        intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
      }
    }
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(outcome, "~", interventions)
    } else {
      paste(outcome, "~", interventions, "+", paste(covariates, collapse = " + "))
    }
    
    formula <- as.formula(formula_str)
    
    # Run appropriate Bayesian model
    result <- switch(method,
                     crt = crtBayes(
                       formula, random = Random, intervention = interventions,
                       nsim = Nsim, data = intervention_data, threshold = Threshold
                     ),
                     mst = mstBayes(
                       formula, random = Random, intervention = interventions,
                       nsim = Nsim, data = intervention_data, threshold = Threshold
                     ),
                     srt = srtBayes(
                       formula, intervention = interventions,
                       nsim = Nsim, data = intervention_data, threshold = Threshold
                     )
    )
    
    output[[as.character(intervention)]] <- result
    
    # Which posterior probability to extract
    prob_key <- if (method == "srt") "Cond" else "Total1"
    prob_es_values[i] <- as.numeric(result$ProbES[[1]][[prob_key]])
    sup_decisions[i] <- if (prob_es_values[i] > SupThreshold) 1 else 0
  }
  
  # Add the reference intervention
  reference_row <- data.frame(
    Intervention  = refintervention,
    ProbES = NA,
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Intervention  = intervention_col,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine and sort
  result <- rbind(reference_row, result)
  result <- result[order(result$Intervention ), ]
  
  return(result)
}
