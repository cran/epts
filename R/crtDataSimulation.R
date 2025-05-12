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
#' @param B0 The intercept of the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The proportion of attrition for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
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
#'   \item{pupils}{Unique pupil ID}
#'   \item{schools}{School ID}
#'   \item{interventions}{Intervention group (0 = control, 1 to `ni` for interventions)}
#'   \item{covariates}{Simulated covariates}
#'   \item{posttest}{Simulated posttest scores (NA if attrited)}
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
#' crtdata <- crtDataSimulation(ni = 3, ns = 10, np = 100, nstreated = c(2, 3, 2, 3), 
#' sigma = 1, ICC = 0.1, B0 = 1.45, es = c(0.1, 0.2, 0.5), 
#' seed = 1234, attritionrates = c(0, 0.1, 0.2, 0.1), covariates = covariates)
#' head(crtdata)
#'
#' @export
# Define the CRT data simulation function
crtDataSimulation <- function(ni, nstreated, np, ns, sigma, ICC, B0, es, seed, attritionrates, covariates) {
  set.seed(seed)
  
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
  
  # Validate covariates structure
  if (!is.list(covariates)) stop("covariates must be a list.")
  for (cov in covariates) {
    if (!is.list(cov)) stop("Each covariate must be a list.")
    if (is.null(cov$name) || is.null(cov$type)) stop("Each covariate must have 'name' and 'type'.")
    if (cov$type == "continuous") {
      if (is.null(cov$sd) || is.null(cov$coefficient)) stop("Continuous covariates must have 'sd' and 'effect'.")
    } else if (cov$type == "categorical") {
      if (is.null(cov$levels) || is.null(cov$probs) || is.null(cov$reference) || is.null(cov$coefficients)) {
        stop("Categorical covariates must have 'levels', 'probs', 'reference', and 'coefficients'.")
      }
    } else {
      stop("Covariate type must be either 'continuous' or 'categorical'.")
    }
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
  
  # Generate covariates (all individuals, no NA)
  for (cov in covariates) {
    if (cov$type == "continuous") {
      data[[cov$name]] <- rnorm(nrow(data), mean = 0, sd = cov$sd)
    } else if (cov$type == "categorical") {
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
  
  # Step 4: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:ni) {
    group_idx <- which(data$interventions == i)
    drop <- round(length(group_idx) * attritionrates[i + 1])
    if (drop > 0) {
      drop_ids <- sample(group_idx, drop)
      data$posttest[drop_ids] <- NA
      non_attrition_idx <- setdiff(non_attrition_idx, drop_ids)
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
  
  
  cov_effects <- rep(0, nrow(data))
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (cov$type == "continuous") {
        cov_effects <- cov_effects + cov$coefficient * data[[cov$name]]
      } else if (cov$type == "categorical") {
        for (lvl in cov$levels) {
          if (lvl != cov$reference) {
            cov_effects <- cov_effects + ifelse(data[[cov$name]] == lvl, cov$coefficients[[lvl]], 0)
          }
        }
      }
    }
  }
  
  # Apply only to non-attrited rows
  # Step 7: Compute post-test scores using the combined formula
  data[non_attrition_idx, posts] <- B0 + cov_effects[non_attrition_idx] + total_treatment_effect + b_i + e_ij
  
  return(data)
}
