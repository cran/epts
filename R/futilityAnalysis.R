#' Futility Analysis Across Interventions for CRT, MST, or SRT Designs
#'
#' This function performs a Bayesian futility analysis for each intervention group compared to control, 
#' across cluster randomized trials (CRT), multisite trials (MST) or simple randomized trials (SRT).
#'
#' @param method The trial design type: "crt", "mst", or "srt".
#' @param data A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention.
#' @param outcome The name of the outcome (post-test) variable. 
#' @param interventions A string specifying the intervention variable.
#' @param Random The name of the clustering variable (e.g., schools or sites) for CRT and MST designs.
#' @param Nsim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param Threshold The effect size threshold for posterior computation (default = 0.05).
#' @param FutThreshold The minimum posterior probability threshold for non-futility (default = 0.8).
#' @param covariates Additional covariates to include in the model. It should be specified as a character vector.
#'
#' @return A `data.frame` with columns:
#' \itemize{
#'   \item `Intervention`: Intervention group identifier.
#'   \item `Futility`: 1 if considered futile (posterior probability < FutThreshold), 0 otherwise.
#'   \item `ProbES`: Bayesian posterior probabilities that the observed effect size is greater than or equal to a pre-specified threshold
#' }
#'
#' @details
#' The function loops over each intervention, fits the appropriate Bayesian model (`crtBayes`, `mstBayes`, `srtBayes`),
#' extracts the posterior probability, and determines futility based on the specified probability threshold.
#'
#' @seealso  \code{\link[eefAnalytics]{crtBayes}}, \code{\link[eefAnalytics]{mstBayes}}, \code{\link[eefAnalytics]{srtBayes}} functions from the \pkg{eefAnalytics} package
#'
#' @examples
#' \donttest{
#' ###Futility analysis of cluster randomized trial###
#' data(crt4armSimData)
#' futilityAnalysis(method = "crt", data = crt4armSimData, outcome = "posttest", 
#' interventions = "interventions", Random = "schools", Nsim = 10000, 
#' Threshold = 0.05, FutThreshold = 0.8)
#'
#' ###Futility analysis of multisite trial###
#' data(mst4armSimData)
#' futilityAnalysis(method = "mst", data = mst4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", Nsim = 10000, 
#' Threshold = 0.05, FutThreshold = 0.8)
#'
#' ###Futility analysis of simple randomized trial###
#' data(srt4armSimData)
#' futilityAnalysis(method = "srt", data = srt4armSimData, outcome = "posttest",
#' interventions = "interventions", Nsim = 10000, Threshold = 0.05, FutThreshold = 0.8)
#'}
#'
#' @importFrom eefAnalytics crtBayes mstBayes srtBayes
#' @export
futilityAnalysis <- function(method = c("crt", "mst", "srt"),
    data,
    outcome = "posttest",
    interventions = "interventions",
    Random = "schools",
    Nsim = 10000,
    Threshold = 0.05,
    FutThreshold = 0.8,
    covariates = NULL) {
  method <- match.arg(method)
  
  intervention_col <- sort(unique(data[[interventions]][data[[interventions]] != 0]))
  prob_es_values <- numeric(length(intervention_col))
  futility_decisions <- numeric(length(intervention_col))
  output <- list()
  
  for (i in seq_along(intervention_col)) {
    intervention <- intervention_col[i]
    intervention_data <- data[data[[interventions]] %in% c(intervention, 0), ]
    intervention_data[[interventions]] <- ifelse(intervention_data[[interventions]] == intervention, 1, 0)
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(outcome, "~", interventions)
    } else {
      paste(outcome, "~", interventions, "+", paste(covariates, collapse = " + "))
    }
    
    formula <- as.formula(formula_str)
    
    # Call the appropriate method
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
    
    output[[i]] <- result
    
    prob_key <- if (method == "srt") "Cond" else "Total1"
    prob_es_values[i] <- as.numeric(result$ProbES[[1]][[prob_key]])
    futility_decisions[i] <- if (prob_es_values[i] < FutThreshold) 1 else 0
  }
  
  return(data.frame(
    Intervention = intervention_col,
    Futility = futility_decisions,
    ProbES = prob_es_values
  ))
}
