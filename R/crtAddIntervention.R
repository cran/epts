#' Add a New Intervention Group to Clustered Randomized Trial (CRT) 
#'
#' This function adds a new intervention group to an existing CRT dataset. It models 
#' post-test outcomes using fixed and random effects estimated from the original data and incorporates 
#' user-specified effect size and attrition for the new intervention.
#'
#' @param originalData A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention for CRT design.
#' @param ns The number of schools to assign to the new intervention group.
#' @param np The number of pupils per new school.
#' @param es The standardized effect size for the new intervention group.
#' @param attritionrate The proportion of pupils in the new group to drop due to attrition.
#' @param outcome A string specifying the name of the column containing outcome variable (e.g., post-test scores). 
#' @param interventions A string specifying the name of the intervention assignment column. 
#' @param schoolsID A string specifying the name of the school ID column. 
#' @param pupilsID A string specifying the name of the pupil ID column. 
#' @param continuous_covariates A character vector specifying the names of continuous covariates in the model.
#' @param categorical_covariates A character vector specifying the names of categorical covariates in the model (converted to factors).
#' @return A `data.frame` combining the original and new intervention group, including post-test outcomes 
#' simulated for the new intervention based on the estimated mixed model.
#'
#' @details
#' The function performs the following:
#' \itemize{
#'   \item Fits a linear mixed-effects model (`lmer`) to the original dataset using provided covariates.
#'   \item Applies the specified effect size (`es`) and generates new post-test scores.
#'   \item Simulates attrition by removing post-test scores at random.
#' }
#'
#' @examples
#' data(crt4armSimData)
#' new_crt5armData <- crtAddIntervention(originalData = crt4armSimData, ns = 2,
#' np = 100, es = 0.3, attritionrate = 0.1, outcome = "posttest", interventions = "interventions", 
#' schoolsID = "schools", pupilsID = "pupils", 
#' continuous_covariates = c("pretest"), categorical_covariates = c("gender", "ethnicity"))
#' head(new_crt5armData)
#' 
#' @importFrom lme4 lmer fixef VarCorr
#' @seealso \code{\link[lme4]{lmer}} from the \pkg{lme4} package
#' @export
# Function to add a new treatment to CRT data
crtAddIntervention <- function(originalData, ns, np, 
                               es, attritionrate, 
                               outcome, interventions, 
                               schoolsID, pupilsID, continuous_covariates, categorical_covariates) {
  
 
    # Check for CRT structure
    school_intervention_check <- aggregate(originalData[[interventions]], 
                                           by = list(originalData[[schoolsID]]), 
                                           FUN = function(x) length(unique(x)))
    if (any(school_intervention_check$x > 1)) {
      stop("Dataset is not a CRT. Each school must have only one unique intervention.")
    }
    
    # Fit LMM
    all_covariates <- c(continuous_covariates, categorical_covariates)
    formula <- as.formula(paste(outcome, "~", paste(all_covariates, collapse = " + "), "+ (1 |", schoolsID, ")"))
    lmer_model <- lmer(formula, data = originalData, REML = TRUE)
    
    fixed_effects <- fixef(lmer_model)
    sigma <- sigma(lmer_model)
    school_sd <- as.numeric(VarCorr(lmer_model)[[schoolsID]])^0.5
    
    new_treatment_num <- max(originalData[[interventions]], na.rm = TRUE) + 1
    new_school_ids <- seq(from = max(originalData[[schoolsID]], na.rm = TRUE) + 1, length.out = ns)
    
    new_data <- expand.grid(
      schls = new_school_ids,
      ppls = 1:np
    )
    
    max_existing_pupil_id <- max(originalData[[pupilsID]], na.rm = TRUE)
    new_data[[pupilsID]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
    new_data[[schoolsID]] <- rep(new_school_ids, each = np)
    new_data[[interventions]] <- new_treatment_num
    
    # Generate categorical covariates
    for (covariate in categorical_covariates) {
      freq_table <- table(originalData[[covariate]])
      probs <- freq_table / sum(freq_table)
      new_data[[covariate]] <- factor(sample(names(probs), size = nrow(new_data), replace = TRUE, prob = probs), levels = names(probs))
    }
    
    # Generate continuous covariates
    for (covariate in continuous_covariates) {
      mu <- mean(originalData[[covariate]], na.rm = TRUE)
      sd_val <- sd(originalData[[covariate]], na.rm = TRUE)
      new_data[[covariate]] <- rnorm(nrow(new_data), mean = mu, sd = sd_val)
    }
    
    # School-level and individual residuals
    new_data$school_effects <- rep(rnorm(ns, mean = 0, sd = school_sd), each = np)
    new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
    
    # Baseline post-test score
    new_data[[outcome]] <- fixed_effects[1]
    
    # Add categorical effects (model.matrix)
    if (length(categorical_covariates) > 0) {
      mat <- model.matrix(~ . - 1, data = new_data[, categorical_covariates, drop = FALSE])
      for (col in colnames(mat)) {
        if (col %in% names(fixed_effects)) {
          new_data[[outcome]] <- new_data[[outcome]] + fixed_effects[col] * mat[, col]
        }
      }
    }
    
    # Add continuous effects
    for (cov in continuous_covariates) {
      if (cov %in% names(fixed_effects)) {
        new_data[[outcome]] <- new_data[[outcome]] + fixed_effects[cov] * new_data[[cov]]
      }
    }
    
    # Add treatment effect
    new_data[[outcome]] <- new_data[[outcome]] +
      (new_data[[interventions]] == new_treatment_num) * es * sqrt(sigma^2 + school_sd^2) +
      new_data$school_effects + new_data$individual_residuals
    
    # Remove temp columns
    new_data <- new_data[, !names(new_data) %in% c("school_effects", "individual_residuals")]
    
    # Apply attrition
    idx <- sample(seq_len(nrow(new_data)), round(nrow(new_data) * attritionrate))
    new_data[[outcome]][idx] <- NA
    
    # Align with existing data
    missing_cols <- setdiff(names(originalData), names(new_data))
    for (col in missing_cols) new_data[[col]] <- NA
    missing_cols2 <- setdiff(names(new_data), names(originalData))
    for (col in missing_cols2) originalData[[col]] <- NA
    new_data <- new_data[names(originalData)]
    
    # Combine
    combined <- rbind(originalData, new_data)
    
    # Drop any auxiliary columns like schls and ppls if they exist
    combined <- combined[, !names(combined) %in% c("schls", "ppls")]
    
    # Reorder columns (optional)
    key_cols <- c(pupilsID, schoolsID, interventions, outcome)
    covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
    final_cols <- c(key_cols, covariate_cols)
    
    # Only include columns that actually exist
    final_cols <- intersect(final_cols, names(combined))
    combined <- combined[, final_cols, drop = FALSE]
    
    return(combined)
    
  }