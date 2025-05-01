#' Plot Posterior Probabilities Across Thresholds for CRT, MST, or SRT Designs
#'
#' This function generates a Bayesian posterior probability plot across multiple thresholds 
#' for each intervention group in a clustered randomized trial (CRT), multisite trial (MST),  
#' or simple randomized trial (SRT).
#'
#' @param method The trial design type: "crt", "mst", or "srt".
#' @param data A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention.
#' @param outcome The name of the outcome (post-test) variable. 
#' @param interventions A string specifying the intervention variable.
#' @param Random The name of the clustering variable (e.g., schools or sites) for CRT and MST designs.
#' @param Nsim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param covariates Additional covariates names to include in the model. It should be specified as a character vector.
#' @param threshold_range The range of thresholds to evaluate. It should be specified as a numeric vector of length 2 (default = c(0, 1.0)).
#' @param VerticalLine Optional vertical reference line added at a threshold value. It should be specified as a numeric value.
#' @param VerticalLineColor The color of the vertical reference line. It should be specified as a character string (default = "#0000FF").
#' @param HorizontalLine Optional posterior probability cutoff for adding a horizontal reference line. It should be specified as a numeric value.
#' @param HorizontalLineColor The color of the horizontal reference line. It should be specified as a character string (default = "#FF0000").
#' @param maintitle The main title of the plot. 
#' @param xlabel The label for the x-axis. 
#' @param ylabel The label for the y-axis.
#' @param intcolors Optional intervention colors specified manually. It should be provided as a named character vector.
#' @param intlabels Optional intervention labels to use instead of default names. It should be specified as a character vector.
#' @param xbreaks Tick marks for the x-axis. Must be a numeric vector with values within the specified threshold_range (default = 0.1). 
#' @param ybreaks Tick marks for the y-axis. It should be specified as a numeric vector (default = seq(0, 1, by = 0.1)).
#'
#' @return A `ggplot` object that displays posterior probabilities across thresholds for each intervention.
#'
#' @details
#' The function uses `crtBayes()`, `mstBayes()`, or `srtBayes()` from eefAnalytics package depending on the `method`.
#'
#' @seealso  \code{\link[eefAnalytics]{crtBayes}}, \code{\link[eefAnalytics]{mstBayes}}, \code{\link[eefAnalytics]{srtBayes}} functions from the \pkg{eefAnalytics} package
#'
#' @examples
#' \donttest{
#'
#' ###Plot Posterior Probabilities of cluster randomized trial###
#' data(crt4armSimData)
#' plotPosteriorProbs(method = "crt",data = crt4armSimData, outcome = "posttest",
#' interventions = "interventions", Random = "schools", Nsim = 10000,
#' threshold_range = c(0, 0.1), VerticalLine = 0.05, HorizontalLine = 0.8,
#' VerticalLineColor= "purple", HorizontalLineColor= "black", 
#' intlabels = c("Intervention A", "Intervention B", "Intervention C"), 
#' intcolors = c("Intervention A" = "blue", "Intervention B" = "red", 
#' "Intervention C" = "green"), maintitle= "Posterior probability plot",
#' xlabel= "Threshold", ylabel= "Posterior probability",
#' xbreaks= 0.1, ybreaks= seq(0, 1, by = 0.1))
#' 
#' ###Plot Posterior Probabilities of multisite trial###
#' data(mst4armSimData)
#' plotPosteriorProbs(method = "ms",data = mst4armSimData, outcome = "posttest", 
#' interventions = "interventions", Random = "schools", Nsim = 10000,
#' threshold_range = c(0, 0.1), VerticalLine = 0.05, HorizontalLine = 0.8, 
#' VerticalLineColor= "purple", HorizontalLineColor= "black",
#' intlabels = c("Intervention A", "Intervention B", "Intervention C"), 
#' intcolors = c("Intervention A" = "blue", "Intervention B" = "red",
#' "Intervention C" = "green"), maintitle= "Posterior probability plot",
#' xlabel= "Threshold", ylabel= "Posterior probability",
#' xbreaks= 0.1, ybreaks= seq(0, 1, by = 0.1))
#'
#' ###Futility analysis of simple randomized trial###
#' data(srt4armSimData)
#' plotPosteriorProbs(method = "srt",data = srt4armSimData, outcome = "posttest",
#' interventions = "interventions", Nsim = 10000, threshold_range = c(0, 0.2),
#' VerticalLine = 0.05, HorizontalLine = 0.8, VerticalLineColor= "purple",
#' HorizontalLineColor= "black", intlabels = c("Intervention A", "Intervention B",
#' "Intervention C"), intcolors = c("Intervention A" = "#1F77B4", 
#' "Intervention B" = "#D62728", "Intervention C" = "#2CA02C"),
#' maintitle= "Posterior probability plot", xlabel= "Threshold", 
#' ylabel= "Posterior probability", xbreaks= 0.1, ybreaks= seq(0, 1, by = 0.1))
#'}
#'
#' @import ggplot2 
#' @import MCMCvis coda  mvtnorm 
#' @importFrom ggpubr theme_pubclean
#' @importFrom graphics abline barplot hist legend lines mtext par title
#' @importFrom methods is
#' @importFrom stats aggregate as.formula model.matrix na.omit rbinom rnorm sd setNames
#' @import ggplot2 
#' @importFrom eefAnalytics crtBayes mstBayes srtBayes
#' @export
plotPosteriorProbs <- function(method = c("crt", "mst", "srt"),
    data, 
    outcome = "posttest", 
    interventions = "interventions", 
    Random = "schools", 
    Nsim = 10000, 
    covariates = NULL, 
    VerticalLine = NULL, 
    VerticalLineColor = "#0000FF", 
    HorizontalLine = NULL, 
    HorizontalLineColor = "#FF0000", 
    threshold_range = c(0, 1.0),
    maintitle = "Posterior Probabilities Across Thresholds", 
    xlabel = "Threshold", 
    ylabel = "Posterior Probability", 
    intcolors = NULL, 
    intlabels = NULL,
    xbreaks = NULL, 
    ybreaks = seq(0, 1, by = 0.1)
) {
  method <- match.arg(method)
  
  if (method == "srt" && (is.null(covariates) || length(covariates) == 0)) {
    data$zero_covariate <- 0
    covariates <- "zero_covariate"
  }
  
  intervention_col <- sort(unique(data[[interventions]][data[[interventions]] != 0]))
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  output <- list()
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(intervention_col))
  
  for (i in seq_along(intervention_col)) {
    for (j in seq_along(thresholds)) {
      intervention <- intervention_col[i]
      intervention_data <- data[data[[interventions]] %in% c(intervention, 0), ]
      intervention_data[[interventions]] <- ifelse(intervention_data[[interventions]] == intervention, 1, 0)
      
      formula_str <- paste(outcome, "~", interventions)
      if (!is.null(covariates)) {
        formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      }
      
      if (method == "crt") {
        output[[j]] <- crtBayes(
          as.formula(formula_str),
          random = Random,
          intervention = interventions,
          nsim = Nsim,
          data = intervention_data,
          threshold = thresholds[j]
        )
        probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
      } else if (method == "mst") {
        output[[j]] <- mstBayes(
          as.formula(formula_str),
          random = Random,
          intervention = interventions,
          nsim = Nsim,
          data = intervention_data,
          threshold = thresholds[j]
        )
        probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
      } else if (method == "srt") {
        interventioncall <- paste0(interventions, 1)
        output[[j]] <- srtBayes(
          as.formula(formula_str),
          intervention = interventions,
          nsim = Nsim,
          data = intervention_data,
          threshold = thresholds[j]
        )
        probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[interventioncall]]["Cond"])
      }
    }
  }
  
  df <- data.frame(
    Threshold = rep(thresholds, length(intervention_col)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(intervention_col)), each = length(thresholds)))
  )
  
  if (!is.null(intlabels)) {
    levels(df$Intervention) <- intlabels
  }
  
  # Generate xbreaks dynamically if not provided
  # Auto-generate xbreaks if it's a single number (step) or NULL
  if (is.null(xbreaks)) {
    step <- 0.1
    xbreaks <- seq(threshold_range[1], threshold_range[2], by = step)
  } else if (length(xbreaks) == 1 && is.numeric(xbreaks)) {
    # If user passed a single number (e.g., xbreaks = 0.1), treat it as step size
    step <- xbreaks
    xbreaks <- seq(threshold_range[1], threshold_range[2], by = step)
  } else {
    # Otherwise, filter xbreaks to within threshold_range
    xbreaks <- xbreaks[xbreaks >= threshold_range[1] & xbreaks <= threshold_range[2]]
  }
  
  
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = maintitle, x = xlabel, y = ylabel) +
    scale_x_continuous(limits = range(xbreaks), breaks = xbreaks)  +
    scale_y_continuous(limits = c(0, 1), breaks = ybreaks) +
    theme_pubclean() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    coord_cartesian(clip = "off")
  
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = VerticalLineColor, size = 1)
  }
  if (!is.null(HorizontalLine)) {
    p <- p + geom_hline(yintercept = HorizontalLine, linetype = "dashed", color = HorizontalLineColor, size = 1)
  }
  if (!is.null(intcolors)) {
    p <- p + scale_color_manual(values = intcolors)
  }
  
  print(p)
}
