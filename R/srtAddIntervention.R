#' Add a New Intervention Group to Simple Randomized Trial (SRT) Data
#'
#' This function adds a new intervention group to an existing SRT dataset by generating 
#' new participant-level data. 
#'
#' @param existing_data A data frame containing the variables including outcome, predictors, the clustering variable, and the intervention for CRT design.
#' @param np The number of new participants to generate for the new intervention group.
#' @param es The standardized effect size for the new intervention group.
#' @param attritionrate The proportion of pupils in the new group to drop due to attrition.
#' @param outcome A string specifying the name of the column containing outcome variable (post-test scores). 
#' @param interventions A string specifying the name of the intervention assignment column. 
#' @param id A string specifying the name of the participant ID column.
#' @param continuous_covariates A character vector specifying the names of continuous covariates.
#' @param categorical_covariates A character vector specifying the names of categorical covariates (converted to factors).
#' @return A `data.frame` combining the original dataset with the newly simulated intervention group.
#'
#'
#' @examples
#' data(srt4armSimData)
#' new_srt5armData <- srtAddIntervention(existing_data = srt4armSimData, np = 100, 
#' es = 0.3, attritionrate = 0.1, outcome = "posttest", interventions = "interventions",
#' id = "ID", continuous_covariates = c("pretest"), categorical_covariates = c("gender", "ethnicity"))
#' head(new_srt5armData)
#'
#' @export
srtAddIntervention <- function(existing_data, np, es, attritionrate, 
                                   outcome, interventions, id, 
                                   continuous_covariates, categorical_covariates) {
    
    sigma <- sd(existing_data[[outcome]], na.rm = TRUE)
    new_treatment_num <- max(existing_data[[interventions]], na.rm = TRUE) + 1
    new_ids <- seq(max(existing_data[[id]], na.rm = TRUE) + 1, length.out = np)
    
    template <- existing_data[1:np, , drop = FALSE]
    template[[id]] <- new_ids
    template[[interventions]] <- new_treatment_num
    
    for (cov in categorical_covariates) {
      freq <- table(existing_data[[cov]])
      probs <- freq / sum(freq)
      template[[cov]] <- factor(sample(names(probs), size = np, replace = TRUE, prob = probs), levels = names(probs))
    }
    
    for (cov in continuous_covariates) {
      mu <- mean(existing_data[[cov]], na.rm = TRUE)
      sd_val <- sd(existing_data[[cov]], na.rm = TRUE)
      template[[cov]] <- rnorm(np, mu, sd_val)
    }
    
    template[[outcome]] <- mean(existing_data[[outcome]], na.rm = TRUE) + 
      es * sigma + 
      rnorm(np, 0, sigma)
    
    idx <- sample(seq_len(nrow(template)), round(nrow(template) * attritionrate))
    template[[outcome]][idx] <- NA
    
    # Combine datasets
    combined <- rbind(existing_data, template)
    
    # Remove internal columns (if any exist)
    combined <- combined[, !names(combined) %in% c("school_effects", "treatment_slopes", "individual_residuals", "schls", "ppls")]
    
    # Reorder columns
    key_cols <- c(id, interventions, outcome)  # no schools_col in SRT
    covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
    final_cols <- intersect(c(key_cols, covariate_cols), names(combined))
    
    combined <- combined[, final_cols, drop = FALSE]
    
    return(combined)
    
  }