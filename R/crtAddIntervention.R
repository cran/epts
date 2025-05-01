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
#' @param covariates A string vector specifying the names of additional covariates (both categorical and continuous) used in the model. 
#'
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
#' new_crt5armData <- crtAddIntervention(originalData = crt4armSimData, ns = 2, np = 100, es = 0.3, 
#' attritionrate = 0.1, outcome = "posttest", interventions = "interventions", schoolsID = "schools", 
#' pupilsID = "pupils", covariates = c("pretest"))
#' head(new_crt5armData)
#' 
#' @importFrom lme4 lmer fixef VarCorr
#' @seealso \code{\link[lme4]{lmer}} from the \pkg{lme4} package
#' @export
# Function to add a new treatment to CRT data
crtAddIntervention <- function(originalData, ns, np, 
                               es, attritionrate, 
                               outcome, interventions, 
                               schoolsID, pupilsID, covariates) {
  
  # --- CRT Design Check ---
  # Each school should have only one unique intervention
  school_intervention_check <- aggregate(originalData[[interventions]], 
                                         by = list(originalData[[schoolsID]]), 
                                         FUN = function(x) length(unique(x)))
  if (any(school_intervention_check$x > 1)) {
    stop("Design check failed: The dataset is not a CRT. Each school must have exactly one unique intervention.")
  }
  
  
  # Step 1: Fit a Linear Mixed-Effects Model using lmer
  formula <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + "), "+ (1 |", schoolsID, ")"))
  lmer_model <- lmer(formula, data = originalData, REML = TRUE)
  
  # Extract fixed effect estimates
  fixed_effects <- fixef(lmer_model)
  
  # Extract variance components
  sigma <- sigma(lmer_model)  # Residual standard deviation
  school_sd <- as.numeric(VarCorr(lmer_model)[[schoolsID]])^0.5  # School-level variance
  
  # Assign the new treatment number as the next available number
  new_treatment_num <- max(originalData[[interventions]], na.rm = TRUE) + 1
  
  # Ensure school IDs continue sequentially
  max_existing_school_id <- max(originalData[[schoolsID]], na.rm = TRUE)
  new_school_ids <- seq(from = max_existing_school_id + 1, length.out = ns)
  
  # Generate the new dataset structure
  new_data <- expand.grid(
    schls = new_school_ids,
    ppls = 1:np
  )
  
  # Ensure pupil IDs start sequentially after the highest existing ID
  max_existing_pupil_id <- max(originalData[[pupilsID]], na.rm = TRUE)
  new_data[[pupilsID]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
  
  # Assign schools column properly
  new_data[[schoolsID]] <- rep(new_school_ids, each = np)
  
  # Assign the new intervention to all new pupils
  new_data[[interventions]] <- new_treatment_num
  
  # Convert categorical covariates to factors in both datasets
  for (covariate in covariates) {
    if (is.character(originalData[[covariate]]) || is.factor(originalData[[covariate]])) {
      originalData[[covariate]] <- as.factor(originalData[[covariate]])
      new_data[[covariate]] <- as.factor(sample(levels(originalData[[covariate]]), nrow(new_data), replace = TRUE))
    }
  }
  
  # Generate numerical covariates (including pretest)
  for (covariate in covariates) {
    if (is.numeric(originalData[[covariate]])) {
      covariate_mean <- mean(originalData[[covariate]], na.rm = TRUE)
      covariate_sd <- sd(originalData[[covariate]], na.rm = TRUE)
      new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
    }
  }
  
  # Generate school-level random effects
  new_school_effects <- rnorm(ns, mean = 0, sd = school_sd)
  
  # Assign school effects to students
  new_data$school_effects <- rep(new_school_effects, each = np)
  
  # Generate individual-level residuals
  new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Compute baseline post-test score (excluding new treatment effect)
  new_data[[outcome]] <- fixed_effects[1]  # Intercept
  
  # Convert categorical variables to dummy variables using model.matrix
  cat_vars <- covariates[sapply(originalData[covariates], is.factor) | sapply(originalData[covariates], is.character)]
  if (length(cat_vars) > 0) {
    cat_matrix <- model.matrix(~ . - 1, data = new_data[, cat_vars, drop = FALSE])
    col_names <- colnames(cat_matrix)
    matched_cols <- intersect(col_names, names(fixed_effects))
    
    for (col in matched_cols) {
      new_data[[outcome]] <- new_data[[outcome]] + fixed_effects[col] * cat_matrix[, col]
    }
  }
  
  # Add continuous covariate effects dynamically
  num_vars <- setdiff(covariates, cat_vars)
  for (covariate in num_vars) {
    if (covariate %in% names(fixed_effects)) {
      new_data[[outcome]] <- new_data[[outcome]] + fixed_effects[covariate] * new_data[[covariate]]
    }
  }
  
  # Apply the effect size **only to the new treatment group**
  effect_size_adjustment <- (new_data[[interventions]] == new_treatment_num) * es * sqrt(sigma^2 + school_sd^2)
  
  # Final post-test scores
  new_data[[outcome]] <- new_data[[outcome]] + effect_size_adjustment + new_data$school_effects + new_data$individual_residuals
  
  # Remove auxiliary columns
  new_data <- new_data[, !names(new_data) %in% c("school_effects", "individual_residuals")]
  
  # Apply attrition randomly by setting some post-test scores to NA
  attrition_size <- round(nrow(new_data) * attritionrate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[outcome]][attrition_idx] <- NA
  
  # Ensure all necessary columns exist in both datasets before combining
  missing_in_new <- setdiff(names(originalData), names(new_data))
  missing_in_existing <- setdiff(names(new_data), names(originalData))
  
  # Add missing columns as NA to new data
  for (col in missing_in_new) {
    new_data[[col]] <- NA
  }
  
  # Add missing columns as NA to existing data
  for (col in missing_in_existing) {
    originalData[[col]] <- NA
  }
  
  # Reorder new_data columns to match originalData
  new_data <- new_data[names(originalData)]
  
  # Combine the new treatment group with the existing data
  final_combined_data <- rbind(originalData, new_data)
  
  # Keep only relevant columns
  final_combined_data <- final_combined_data[, c(pupilsID, schoolsID, outcome, interventions, covariates)]
  
  return(final_combined_data)
}
