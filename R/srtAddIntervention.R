#' Add a New Intervention Group to Simple Randomized Trial (SRT) Data
#'
#' This function adds a new intervention group to an existing SRT dataset by generating 
#' new participant-level data. 
#'
#' @param originalData A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention for CRT design.
#' @param np The number of new participants to generate for the new intervention group.
#' @param es The standardized effect size for the new intervention group.
#' @param attritionrate The proportion of pupils in the new group to drop due to attrition.
#' @param outcome A string specifying the name of the column containing outcome variable (post-test scores). 
#' @param interventions A string specifying the name of the intervention assignment column. 
#' @param id A string specifying the name of the participant ID column.
#' @param covariates A string vector specifying the names of additional covariates (both categorical and continuous) used in the model.
#'
#' @return A `data.frame` combining the original dataset with the newly simulated intervention group.
#'
#'
#' @examples
#' data(srt4armSimData)
#' new_srt5armData <- srtAddIntervention(originalData = srt4armSimData, np = 100, 
#' es = 0.3, attritionrate = 0.1, outcome = "posttest", interventions = "interventions",
#' id = "ID", covariates = c("pretest"))
#' head(new_srt5armData)
#'
#' @export
srtAddIntervention <- function(originalData, np, es, attritionrate, 
                                 outcome, interventions, id, covariates) {
  
  # Extract standard deviation from existing data
  sigma <- sd(originalData[[outcome]], na.rm = TRUE)
  
  # Determine new treatment number
  new_treatment_num <- ifelse(is.na(max(originalData[[interventions]], na.rm = TRUE)), 
                              1, 
                              max(originalData[[interventions]], na.rm = TRUE) + 1)
  
  # Select existing ID column
  last_pupil_id <- max(originalData[[id]], na.rm = TRUE)
  last_pupil_id <- ifelse(is.na(last_pupil_id), 0, last_pupil_id)
  new_pupil_ids <- seq(from = last_pupil_id + 1, length.out = np)
  
  # Create new dataset, keeping only selected variables
  new_data <- originalData[1:np, ]  # Duplicate the structure of existing data
  new_data[[id]] <- new_pupil_ids  # Assign new IDs from existing column
  new_data[[interventions]] <- new_treatment_num  # Assign new intervention
  
  # Handle covariates (keep same distribution as existing data)
  for (covariate in covariates) {
    if (covariate %in% names(originalData)) {
      unique_values <- unique(na.omit(originalData[[covariate]]))
      
      if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
        # Binary Variable: Sample using the same proportion from existing data
        prob_1 <- mean(originalData[[covariate]] == 1, na.rm = TRUE)
        new_data[[covariate]] <- rbinom(nrow(new_data), size = 1, prob = prob_1)
        
      } else if (is.numeric(originalData[[covariate]]) && length(unique_values) > 10) {
        # Continuous Variable: Generate using normal distribution
        covariate_mean <- mean(originalData[[covariate]], na.rm = TRUE)
        covariate_sd <- sd(originalData[[covariate]], na.rm = TRUE)
        new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
        
      } else {
        # Categorical Variable: Sample based on observed proportions
        value_counts <- table(originalData[[covariate]])
        category_probs <- value_counts / sum(value_counts)
        new_data[[covariate]] <- sample(names(category_probs), size = nrow(new_data), replace = TRUE, prob = category_probs)
      }
    } else {
      new_data[[covariate]] <- NA  # If covariate not in existing data, set NA
    }
  }
  
  # Generate post-test scores based on effect size and existing distribution
  effect_size_adjustment <- es * sigma
  new_data[[outcome]] <- mean(originalData[[outcome]], na.rm = TRUE) + effect_size_adjustment + 
    rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Apply attrition randomly by setting some post-test scores to NA
  attrition_size <- round(nrow(new_data) * attritionrate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[outcome]][attrition_idx] <- NA
  
  # Ensure consistency with existing data structure (remove unwanted columns)
  new_data <- new_data[names(originalData)]
  
  # Combine with existing dataset
  final_combined_data <- rbind(originalData, new_data)
  
  return(final_combined_data)
}
