#' Add a New Intervention Group to Multisite Trial (MST) 
#'
#' This function adds a new intervention group to an existing Multisite Trial (MST) dataset. It fits a 
#' linear mixed-effects model to the original data, then uses its estimates to generate post-test
#' outcomes for the new group, incorporating random intercepts, slopes, and user-defined effect size.
#'
#' @param originalData A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention for MST design.
#' @param ns The number of schools to assign to the new intervention group.
#' @param np The number of pupils per new school.
#' @param es The standardized effect size for the new intervention group.
#' @param attritionrate The proportion of pupils in the new group to drop due to attrition.
#' @param intper Proportion of pupils per new school assigned to the Intervention group.
#' @param outcome A string specifying the name of the column containing outcome variable (e.g., post-test scores). 
#' @param interventions A string specifying the name of the intervention assignment column. 
#' @param schoolsID A string specifying the name of the school ID column. 
#' @param pupilsID A string specifying the name of the pupil ID column. 
#' @param continuous_covariates A character vector specifying the names of continuous covariates.
#' @param categorical_covariates A character vector specifying the names of categorical covariates (converted to factors). 
#'
#' @return A `data.frame` containing the combined dataset with the newly added intervention group and simulated outcomes.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Fits a linear mixed-effects model (`lmer`) with random slopes and intercepts using existing MST data.
#'   \item Simulates new schools and pupils, assigning intervention randomly by specified percentage. 
#'   \item Simulates attrition by removing post-test scores at random.
#' }
#' 
#'
#' @examples
#' data(mst4armSimData)
#' new_mst5armData <- mstAddIntervention(originalData = mst4armSimData, ns = 2, np = 100, es = 0.3,
#' intper = 0.5, attritionrate = 0.1, outcome = "posttest", interventions = "interventions",
#' schoolsID = "schools", pupilsID = "pupils", 
#' continuous_covariates = c("pretest"), categorical_covariates = c("gender", "ethnicity"))
#' head(new_mst5armData)
#'
#' @importFrom lme4 lmer fixef VarCorr
#' @seealso \code{\link[lme4]{lmer}} from the \pkg{lme4} package
#' @export
mstAddIntervention<- function(originalData, ns, np, 
                                 es, attritionrate, intper,
                                 outcome, interventions, 
                                 schoolsID, pupilsID, continuous_covariates, categorical_covariates) {
  

    all_covariates <- c(continuous_covariates, categorical_covariates)
    formula <- as.formula(paste(outcome, "~", paste(all_covariates, collapse = " + "), "+ (1 +", interventions, "|", schoolsID, ")"))
    model <- lmer(formula, data = originalData, REML = TRUE)
    
    sigma <- sigma(model)
    rand_sd <- attr(VarCorr(model)[[schoolsID]], "stddev")
    sigmab0 <- rand_sd[1]
    sigmab1 <- ifelse(length(rand_sd) > 1, rand_sd[2], 0)
    fixed <- fixef(model)
    
    new_treatment_num <- max(originalData[[interventions]], na.rm = TRUE) + 1
    new_school_ids <- seq(max(originalData[[schoolsID]], na.rm = TRUE) + 1, length.out = ns)
    
    new_data <- expand.grid(schls = new_school_ids, ppls = 1:np)
    names(new_data)[1] <- schoolsID
    new_data[[pupilsID]] <- seq(max(originalData[[pupilsID]], na.rm = TRUE) + 1, length.out = nrow(new_data))
    
    # Assign intervention
    new_data[[interventions]] <- unlist(lapply(seq_along(new_school_ids), function(i) {
      treated <- round(np * intper)
      c(rep(new_treatment_num, treated), rep(0, np - treated))
    }))
    
    # Generate covariates
    for (cov in categorical_covariates) {
      freq <- table(originalData[[cov]])
      probs <- freq / sum(freq)
      new_data[[cov]] <- factor(sample(names(probs), nrow(new_data), replace = TRUE, prob = probs), levels = names(probs))
    }
    for (cov in continuous_covariates) {
      mu <- mean(originalData[[cov]], na.rm = TRUE)
      sd_val <- sd(originalData[[cov]], na.rm = TRUE)
      new_data[[cov]] <- rnorm(nrow(new_data), mu, sd_val)
    }
    
    # Random effects
    school_intercepts <- rnorm(ns, 0, sigmab0)
    slopes <- rnorm(ns, 0, sigmab1)
    new_data$school_effects <- rep(school_intercepts, each = np)
    new_data$treatment_slopes <- rep(slopes, each = np)
    new_data$individual_residuals <- rnorm(nrow(new_data), 0, sigma)
    
    new_data[[outcome]] <- fixed[1] + new_data$school_effects + 
      new_data$treatment_slopes * new_data[[interventions]] + 
      new_data$individual_residuals + 
      (new_data[[interventions]] == new_treatment_num) * es * sqrt(sigma^2 + sigmab0^2 + sigmab1^2)
    
    # Apply attrition
    idx <- sample(seq_len(nrow(new_data)), round(nrow(new_data) * attritionrate))
    new_data[[outcome]][idx] <- NA
    
    # Final structure
    # Combine datasets
    # Ensure new_data has all columns in originalData
    missing_in_new <- setdiff(names(originalData), names(new_data))
    for (col in missing_in_new) {
      new_data[[col]] <- NA
    }
    
    # Ensure originalData has all columns in new_data
    missing_in_existing <- setdiff(names(new_data), names(originalData))
    for (col in missing_in_existing) {
      originalData[[col]] <- NA
    }
    
    # Reorder columns to match
    new_data <- new_data[names(originalData)]
    
    # Combine safely
    combined <- rbind(originalData, new_data)
    
    # Remove temp/internal columns
    combined <- combined[, !names(combined) %in% c("school_effects", "treatment_slopes", "individual_residuals", "schls", "ppls")]
    
    # Reorder columns
    key_cols <- c(pupilsID, schoolsID, interventions, outcome)
    covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
    final_cols <- intersect(c(key_cols, covariate_cols), names(combined))
    
    combined <- combined[, final_cols, drop = FALSE]
    
    return(combined)
    
  }
  