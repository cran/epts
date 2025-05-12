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
#' @param B0 The intercept of the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The attrition proportions for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
#' @param covariates List of covariate specifications. Each element should be a list with the following fields:
#'   \describe{
#'     \item{name}{Character. Name of the covariate.}
#'     \item{type}{Character. Either \code{"continuous"} or \code{"categorical"}.}
#'     \item{sd}{Numeric. Standard deviation (only for continuous covariates).}
#'     \item{coefficient}{Numeric. Coefficient (only for continuous covariates).}
#'     \item{levels}{Character vector. Category levels (only for categorical covariates).}
#'     \item{probs}{Numeric vector. Sampling probabilities (must sum to 1) (categorical only).}
#'     \item{reference}{Character. Reference category (categorical only).}
#'     \item{coefficients}{Named list of numeric values. Coefficients for each non-reference level.}
#'   }
#'
#' @return A `data.frame` containing:
#' \describe{
#'   \item{pupils}{Pupil ID}
#'   \item{schools}{School ID}
#'   \item{interventions}{Intervention group assignment (0 = control, 1 to `ni` = intervention groups)}
#'   \item{covariates}{Simulated covariates}
#'   \item{posttest}{Posttest score (NA if attrited)}
#' }
#'
#' @examples
#'covariates <- list(
#'  list(name = "pretest", type = "continuous", sd = 1, coefficient = 1.7),
#'  list(name = "gender", type = "categorical", levels = c("Male", "Female"), 
#'  probs = c(0.3, 0.7), reference = "Male", coefficients = list(B = -0.5)),
#'  list(name = "ethnicity", type = "categorical", levels = c("White", "Black", "Asian"),
#'  probs = c(0.3, 0.3, 0.4), reference = "White", coefficients = list(B = 1.02, C = 1.3))
#')
#'
#' mstdata <- mstDataSimulation(ni = 3, ns = 10, np = 100, tpi = c(30, 30, 20, 20),
#' sigma = 1, sigmab0 = 0.5, sigmab1 = 0.5, B0 = 1.45, 
#' es = c(0.2, 0.3, 0.1), seed = 1234, attritionrates = c(0.1, 0.1, 0.1, 0), covariates = covariates) 
#' head(mstdata)
#'
#' @export
mstDataSimulation <- function(ni, tpi, np, ns, sigma, sigmab0, sigmab1, B0, es, seed, attritionrates, covariates) {
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
  remaining_n <- length(remaining_pupils)
  
  if (ni == 1) {
    # Only one treatment group: assign all remaining pupils to it
    data[[interventions]][remaining_pupils] <- 1
  } else {
    # More than one treatment group
    treatment_props <- tpi[2:(ni+1)] / sum(tpi[2:(ni+1)])  # Normalize proportions
    treatment_sizes <- round(remaining_n * treatment_props)
    
    # Adjust last group to absorb rounding error
    treatment_sizes[ni] <- remaining_n - sum(treatment_sizes[1:(ni-1)])
    
    for (i in 1:ni) {
      treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
      data[[interventions]][treated_pupils] <- i
      remaining_pupils <- setdiff(remaining_pupils, treated_pupils)
    }
  }
  
  
  for (cov in covariates) {
    if (cov$type == "continuous") {
      data[[cov$name]] <- rnorm(nrow(data), mean = 0, sd = cov$sd)
    } else {
      data[[cov$name]] <- sample(cov$levels, nrow(data), replace = TRUE, prob = cov$probs)
    }
  }
  
  # Convert categorical covariates to numeric (0 = reference)
  for (cov in covariates) {
    if (cov$type == "categorical") {
      levels_ordered <- c(cov$reference, setdiff(cov$levels, cov$reference))
      data[[cov$name]] <- match(data[[cov$name]], levels_ordered) - 1
    }
  }
  
  
  # Initialize post-test scores
  posts <- "posttest"
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
  
  
  if (length(unique_interventions) > 1) {
    treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
    treatment_matrix <- treatment_matrix[, -1, drop = FALSE]  # Remove control column, keep matrix structure
    
    if (is.null(dim(treatment_matrix))) {
      treatment_matrix <- matrix(treatment_matrix, ncol = 1)
    }
    
    treatment_effects <- sweep(treatment_matrix, 2, es * sqrt(sigmab0^2 + sigmab1^2 + sigma^2), `*`)
    total_treatment_effect <- rowSums(treatment_effects)
    random_slope <- treatment_matrix * b1_i
    total_random_slope_effect <- rowSums(random_slope)
  } else {
    total_treatment_effect <- rep(0, length(non_attrition_idx))
    total_random_slope_effect <- rep(0, length(non_attrition_idx))
  }
  
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects) 
  
  random_slope <-  treatment_matrix*b1_i
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_random_slope_effect <- rowSums(random_slope)
  
  cov_effects <- rep(0, nrow(data))
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (cov$type == "continuous") {
        cov_effects <- cov_effects + cov$coefficient * data[[cov$name]]
      } else {
        for (lvl in cov$levels) {
          if (lvl != cov$reference) {
            cov_effects <- cov_effects + ifelse(data[[cov$name]] == lvl, cov$coefficients[[lvl]], 0)
          }
        }
      }
    }
  }
  
  # Step 8: Compute post-test scores for non-attrited participants
  data[non_attrition_idx, posts] <- B0 + 
    cov_effects[non_attrition_idx] + 
    total_treatment_effect + 
    b_i + total_random_slope_effect + 
    e_ij
  return(data)
}