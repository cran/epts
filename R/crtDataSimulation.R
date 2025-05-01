#' Simulate Clustered Randomized Trial (CRT) Data
#'
#' This function simulates a multiple intervention arms CRT data. The model
#' includes intervention and pre-test scores as covariates.
#' 
#'
#' @param ni The number of intervention groups excluding the control group.
#' @param nstreated The number of schools in each group, including the control group. It should be specified as an integer vector of length ni + 1.
#' @param np The number of pupils per school.
#' @param ns The total number of schools.
#' @param sigma The standard deviation of the individual-level error.
#' @param ICC The intra-class correlation coefficient.
#' @param sigmaPret The standard deviation of the pretest scores.
#' @param B0 The intercept of the model.
#' @param B1 The coefficient for the pretest covariate in the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The proportion of attrition for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
#'
#' @return A `data.frame` containing:
#' \describe{
#'   \item{pupils}{Unique pupil ID}
#'   \item{schools}{School ID}
#'   \item{interventions}{Intervention group (0 = control, 1 to `ni` for interventions)}
#'   \item{pretest}{Simulated pretest scores}
#'   \item{posttest}{Simulated posttest scores (NA if attrited)}
#' }
#'
#' @examples
#' crtdata <- crtDataSimulation(ni = 3, ns = 10, np = 100, nstreated = c(2, 3, 2, 3), 
#' sigma = 1, ICC = 0.3, sigmaPret = 1, B0 = 1.45, B1 = 1.7, es = c(0.1, 0.2, 0.5), 
#' seed = 1234, attritionrates = c(0, 0.1, 0.2, 0.1))
#' head(crtdata)
#'
#' @export
# Define the CRT data simulation function
crtDataSimulation <- function(ni, ns, np, nstreated, sigma, ICC, sigmaPret, B0, B1, es, seed, attritionrates) {
  
  # Error checking
  if (length(nstreated) != ni + 1) {
    stop("Error: The length of nstreated must be number of interventions + 1 (including control group).")
  }
  if (length(es) != ni) {
    stop("Error: The length of es must be equal to number of interventions.")
  }
  if (length(attritionrates) != ni + 1) {
    stop("Error: The length of attritionrates must be number of interventions + 1 (including control group).")
  }
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Step 2: Initialize treatment assignment and assign treatment groups to schools
  interventions <- "interventions"
  unique_schools <- unique(data$schools)
  
  # Initialize treatment column (0 = control)
  data[[interventions]] <- 0
  
  # Check if the total number of schools assigned matches the number of schools (ns)
  if (sum(nstreated) != ns) {
    stop("Error: The total number of schools in nstreated must be equal to the total number of schools (ns).")
  }
  
  # Assign schools to control group and treatment groups
  remaining_schools <- unique_schools
  
  # Step 2a: Assign schools to the control group
  control_schools <- sample(remaining_schools, nstreated[1])
  remaining_schools <- setdiff(remaining_schools, control_schools)  # Remove control schools from remaining
  
  # Step 2b: Assign schools to each treatment group
  for (i in 1:ni) {
    treated_schools <- sample(remaining_schools, nstreated[i + 1])
    data[[interventions]][data$schools %in% treated_schools] <- i  # Assign treatment group i
    remaining_schools <- setdiff(remaining_schools, treated_schools)  # Remove assigned schools
  }
  
  # Step 3: Generate pre-test scores for everyone
  prets <- "pretest"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)  # Pre-test scores for everyone
  
  # Initialize post-test scores
  posts <- "posttest"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 4: Handle attrition for control and treatment groups
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
  
  # Step 5: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  sigmab <- sqrt(ICC / (1 - ICC) * sigma^2)  # School-level variance for all schools
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab)  # Cluster-level random effects for all schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]  # Get the random effects for relevant schools
  
  # Step 6: Create treatment effects for each group
  treatment_effects <- sapply(1:ni, function(i) as.numeric(data[non_attrition_idx, interventions] == i))
  
  treatment_effects_sizes <- sweep(treatment_effects, 2, es* sqrt(sigmab^2 + sigma^2), `*`)
  
  # Calculate the total treatment effect by combining individual treatment effects with corresponding effect sizes
  total_treatment_effect <- rowSums(treatment_effects_sizes)
  
  # Step 7: Compute post-test scores using the combined formula
  data[non_attrition_idx, posts] <- B0 + B1 * data[non_attrition_idx, prets] + total_treatment_effect + b_i + e_ij
  
  # Only keep relevant columns (pupils, schools, interventions, pretest, posttest)
  data <- data[, c("pupils", "schools", interventions, prets, posts)]
  
  return(data)
}