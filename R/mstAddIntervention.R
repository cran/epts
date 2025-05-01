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
#' @param covariates The names of additional covariates (both categorical and continuous) used in the model. 
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
#' schoolsID = "schools", pupilsID = "pupils", covariates = c("pretest"))
#' head(new_mst5armData)
#'
#' @importFrom lme4 lmer fixef VarCorr
#' @seealso \code{\link[lme4]{lmer}} from the \pkg{lme4} package
#' @export
mstAddIntervention<- function(originalData, ns, np, 
                                 es, attritionrate, intper,
                                 outcome, interventions, 
                                 schoolsID, pupilsID, covariates) {
  
  
  # --- MST Design Check ---
  # At least one school should have more than one unique intervention
  school_intervention_check <- aggregate(originalData[[interventions]], 
                                         by = list(originalData[[schoolsID]]), 
                                         FUN = function(x) length(unique(x)))
  if (all(school_intervention_check$x <= 1)) {
    stop("Design check failed: The dataset is not an MST. At least one school must contain more than one intervention group.")
  }
  
  # Fit hierarchical model with random slopes and intercepts
  formula <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + "), "+ (1 +", interventions, "|", schoolsID, ")"))
  lmer_model <- lmer(formula, data = originalData, REML = TRUE)
  
  # Extract fixed effects and variance components
  fixed_effects <- fixef(lmer_model)
  random_effects <- VarCorr(lmer_model)
  
  sigma <- sigma(lmer_model)  # Residual standard deviation
  sigmab0 <- attr(random_effects[[schoolsID]], "stddev")[1]  # Random intercept SD
  
  # Extract random slope SD (handle missing cases)
  sigmab1 <- ifelse(length(attr(random_effects[[schoolsID]], "stddev")) > 1, 
                    attr(random_effects[[schoolsID]], "stddev")[2], 
                    0)  # If missing, set to 0
  
  # Assign new treatment number
  new_treatment_num <- ifelse(is.na(max(originalData[[interventions]], na.rm = TRUE)), 
                              1, 
                              max(originalData[[interventions]], na.rm = TRUE) + 1)
  
  # Assign new school IDs
  last_school_id <- max(originalData[[schoolsID]], na.rm = TRUE)
  last_school_id <- ifelse(is.na(last_school_id), 0, last_school_id)
  new_school_ids <- seq(from = last_school_id + 1, length.out = ns)
  
  # Create new data structure
  new_data <- expand.grid(
    schls = new_school_ids,
    ppls = 1:np
  )
  
  # Rename the school column
  names(new_data)[names(new_data) == "schls"] <- schoolsID
  
  # Assign unique pupil IDs
  max_existing_pupil_id <- max(originalData[[pupilsID]], na.rm = TRUE)
  max_existing_pupil_id <- ifelse(is.na(max_existing_pupil_id), 0, max_existing_pupil_id)
  new_data[[pupilsID]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
  
  # Assign treatment within schools
  new_data[[interventions]] <- unlist(lapply(new_school_ids, function(x) {
    treatment_group_size <- round(np * intper)
    group_assignment <- c(rep(new_treatment_num, treatment_group_size), 
                          rep(0, np - treatment_group_size))
    sample(group_assignment)
  }))
  
  # Generate school-level random effects
  new_school_intercepts <- rnorm(ns, mean = 0, sd = sigmab0)
  new_treatment_slopes <- rnorm(ns, mean = 0, sd = sigmab1)
  
  # Handle covariates properly
  for (covariate in covariates) {
    if (covariate %in% names(originalData)) {
      unique_values <- unique(na.omit(originalData[[covariate]]))
      
      if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
        prob_1 <- mean(originalData[[covariate]] == 1, na.rm = TRUE)
        new_data[[covariate]] <- rbinom(nrow(new_data), size = 1, prob = prob_1)
        
      } else if (is.numeric(originalData[[covariate]]) && length(unique_values) > 10) {
        covariate_mean <- mean(originalData[[covariate]], na.rm = TRUE)
        covariate_sd <- sd(originalData[[covariate]], na.rm = TRUE)
        new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
        
      } else {
        value_counts <- table(originalData[[covariate]])
        category_probs <- value_counts / sum(value_counts)
        new_data[[covariate]] <- sample(names(category_probs), size = nrow(new_data), replace = TRUE, prob = category_probs)
      }
    } else {
      new_data[[covariate]] <- NA
    }
  }
  
  # Assign school effects to students
  new_data$school_effects <- rep(new_school_intercepts, each = np)
  new_data$treatment_slopes <- rep(new_treatment_slopes, each = np)
  new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Compute post-test scores
  new_data[[outcome]] <- fixed_effects[1] +  # Intercept
    new_data$school_effects + 
    new_data$treatment_slopes * new_data[[interventions]] + 
    new_data$individual_residuals
  
  # Apply effect size **only to the new treatment group**
  effect_size_adjustment <- (new_data[[interventions]] == new_treatment_num) * 
    (es * sqrt(sigma^2 + sigmab0^2 + sigmab1^2))
  
  new_data[[outcome]] <- new_data[[outcome]] + effect_size_adjustment
  
  # Apply attrition
  attrition_size <- round(nrow(new_data) * attritionrate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[outcome]][attrition_idx] <- NA
  
  # Merge data
  new_data <- new_data[, names(originalData)]
  final_combined_data <- rbind(originalData, new_data)
  
  return(final_combined_data)
}