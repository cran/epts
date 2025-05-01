#' Simulate Multisite Trial (MST) Data
#'
#' This function simulates a multiple intervention arms Multisite Trial (MST) data. The model 
#' includes intervention and pre-test scores as covariates.
#'
#' @param ni The number of intervention groups excluding the control group.
#' @param tpi The proportions (in percent) of total participants assigned to each group, with the first value for the control group. It should be specified as a numeric vector of length ni + 1.
#' @param np The number of pupils per school.
#' @param ns The number of schools.
#' @param sigma The standard deviation of the individual-level error.
#' @param sigmab0 The standard deviation of random intercepts at the school level.
#' @param sigmab1 The standard deviation of random slopes for the intervention effect.
#' @param sigmaPret The standard deviation of the pretest scores.
#' @param B0 The intercept of the model.
#' @param B1 The coefficient for the pretest score in the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The attrition proportions for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
#'
#' @return A `data.frame` containing:
#' \describe{
#'   \item{pupils}{Pupil ID}
#'   \item{schools}{School ID}
#'   \item{interventions}{Intervention group assignment (0 = control, 1 to `ni` = intervention groups)}
#'   \item{pretest}{Pretest score}
#'   \item{posttest}{Posttest score (NA if attrited)}
#' }
#'
#' @examples
#' mstdata <- mstDataSimulation(ni = 3, ns = 10, np = 100, tpi = c(30, 30, 20, 20),
#' sigma = 1, sigmab0 = 0.5, sigmab1 = 0.5, sigmaPret = 1, B0 = 0, B1 = 0.5,
#' es = c(0.2, 0.3, 0.1), seed = 1234, attritionrates = c(0.1, 0.1, 0.1, 0)) 
#' head(mstdata)
#'
#' @export
mstDataSimulation <- function(ni, ns, np, tpi, sigma, sigmab0, sigmab1, sigmaPret, B0, B1, es, seed, attritionrates) {
  # Error checking: ensure tpi has length ni + 1 (control + treatment groups)
  if (length(tpi) != ni + 1) {
    stop("Error: 'tpi' must have length ni + 1 (first value for control group, remaining for treatment groups).")
  }
  
  # Error checking: ensure the sum of tpi is 100
  if (sum(tpi) != 100) {
    stop("Error: The sum of 'tpi' must be 100%.")
  }
  
  # Error checking: ensure es has length equal to ni
  if (length(es) != ni) {
    stop("Error: 'es' must have length ni (one for each treatment group).")
  }
  
  # Error checking: ensure attritionrates has length equal to ni + 1
  if (length(attritionrates) != ni + 1) {
    stop("Error: 'attritionrates' must have length ni + 1 (one for control group and one for each treatment group).")
  }
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Total number of participants
  total_participants <- nrow(data)
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- floor(total_participants * (tpi[1] / 100))
  control_pupils <- sample(1:total_participants, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups to the remaining pupils
  remaining_pupils <- setdiff(1:total_participants, control_pupils)
  
  # Calculate the number of participants for each treatment group except the last one
  treatment_sizes <- floor(total_participants * (tpi[2:(ni+1)] / 100))
  
  # Adjust the last treatment group to take the remaining participants
  treatment_sizes[ni] <- length(remaining_pupils) - sum(treatment_sizes[1:(ni-1)])
  
  for (i in 1:ni) {
    treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
    data[[interventions]][treated_pupils] <- i  # Assign treatment group i to selected pupils
    remaining_pupils <- setdiff(remaining_pupils, treated_pupils)  # Remove assigned pupils
  }
  
  # Step 4: Generate pre-test scores for everyone
  prets <- "pret"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)
  
  # Initialize post-test scores
  posts <- "post"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 5: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:ni) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attritionrates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 6: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  
  # Cluster-level random effects
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab0)  # Random intercepts for schools
  b1_i_full <- rnorm(ns, mean = 0, sd = sigmab1)  # Random slopes for treatment effect
  
  # Get the random effects for the relevant schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]
  b1_i <- b1_i_full[data$schools[non_attrition_idx]]
  
  # Step 7: Use model.matrix to create treatment effect dummy variables matrix
  unique_interventions <- unique(data[non_attrition_idx, interventions])
  
  # Create the model matrix without intercept
  treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
  
  # Remove the first column
  treatment_matrix <- treatment_matrix[, -1]
  
  treatment_effects <-  sweep(treatment_matrix, 2, es* sqrt(sigmab0^2 + sigmab1^2 + sigma^2), `*`)
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects) 
  
  random_slope <-  treatment_matrix*b1_i
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_random_slope_effect <- rowSums(random_slope) 
  
  # Step 8: Compute post-test scores for non-attrited participants
  data[non_attrition_idx, posts] <- B0 + 
    B1 * data[non_attrition_idx, prets] + 
    total_treatment_effect + 
    b_i + total_random_slope_effect + 
    e_ij
  
  # Only keep relevant columns (pupils, schools, interventions, pret, post)
  data <- data[, c("pupils", "schools", interventions, prets, posts)]
  
  return(data)
}

