#' Simulate Simple Randomized Trial (SRT) Data
#'
#' This function simulates a Simple Randomized Trial (SRT), with multiple intervention arms, 
#' pre-test and post-test scores, and individual-level attrition. No clustering or hierarchical structure is assumed.
#'
#' @param ni The number of intervention groups excluding the control group.
#' @param tpi The proportions (in percent) assigned to each group, with the first value for the control group followed by the intervention groups. Must sum to 100. It should be specified as a numeric vector of length ni + 1.
#' @param np The total number of participants.
#' @param sigma The standard deviation of individual-level error for the post-test score.
#' @param B0 The intercept term in the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The attrition rates for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
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
#'   \item{ID}{Participant ID}
#'   \item{interventions}{Intervention assignment (0 = control, 1 to `ni` = intervention groups)}
#'   \item{covariates}{Simulated covariates}
#'   \item{posttest}{Posttest score (NA if participant attrited)}
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
#' srtdata <- srtDataSimulation(ni = 3, np = 1000, tpi = c(30, 30, 20, 20),
#' sigma = 1, B0 = 1.45,  es = c(0.2, 0.3, 0.1), seed = 1234,
#' attritionrates = c(0.1, 0.1, 0.1, 0), covariates = covariates) 
#' head(srtdata)
#'
#' @export
srtDataSimulation <- function(ni, tpi, np, sigma, B0, es, seed, attritionrates, covariates) {
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
  
  # Step 1: Create the base data structure with participants
  data <- data.frame(ID = 1:np)  # Assign unique IDs
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- round(np * (tpi[1] / 100))
  control_pupils <- sample(1:np, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups
  remaining_pupils <- setdiff(1:np, control_pupils)
  remaining_n <- length(remaining_pupils)
  
  if (ni == 1) {
    # Only one treatment group: assign everyone left
    data[[interventions]][remaining_pupils] <- 1
  } else {
    # Multiple treatment groups
    treatment_props <- tpi[2:(ni+1)] / sum(tpi[2:(ni+1)])  # Normalize proportions
    treatment_sizes <- round(remaining_n * treatment_props)
    
    # Adjust last group to take any rounding error
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
  
  # Step 6: Generate post-test scores for non-attrited participants
  
  # Check how many unique intervention groups remain after attrition
  unique_interventions <- unique(data[non_attrition_idx, interventions])
  
  if (length(unique_interventions) > 1) {
    # Create treatment matrix with dummy variables (no intercept)
    treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
    
    # Remove the first column (control group is baseline)
    if (ncol(treatment_matrix) >= 1) {
      treatment_matrix <- treatment_matrix[, -1, drop = FALSE]
    }
    
    if (is.null(dim(treatment_matrix))) {
      treatment_matrix <- matrix(treatment_matrix, ncol = 1)
    }
    
    # Apply treatment effects using only individual-level variance (sigma)
    treatment_effects <- sweep(treatment_matrix, 2, es * sigma, `*`)
    
    # Calculate the total treatment effect
    total_treatment_effect <- rowSums(treatment_effects)
  } else {
    # No treatment groups left after attrition
    total_treatment_effect <- rep(0, length(non_attrition_idx))
  }
  
  # Generate individual-level errors for non-attrited individuals
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)
  
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
  
  # Compute post-test scores
  data[non_attrition_idx, posts] <- B0 + 
    cov_effects[non_attrition_idx] + 
    total_treatment_effect + 
    e_ij
  return(data)
}