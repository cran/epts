#' Simulate Simple Randomized Trial (SRT) Data
#'
#' This function simulates a Simple Randomized Trial (SRT), with multiple intervention arms, 
#' pre-test and post-test scores, and individual-level attrition. No clustering or hierarchical structure is assumed.
#'
#' @param ni The number of intervention groups excluding the control group.
#' @param tpi The proportions (in percent) assigned to each group, with the first value for the control group followed by the intervention groups. Must sum to 100. It should be specified as a numeric vector of length ni + 1.
#' @param np The total number of participants.
#' @param sigma The standard deviation of individual-level error for the post-test score.
#' @param sigmaPret The standard deviation of the pretest scores.
#' @param B0 The intercept term in the model.
#' @param B1 The coefficient for the pretest in the model.
#' @param es The standardized effect sizes for each intervention group. It should be specified as a numeric vector.
#' @param seed The random seed for reproducibility.
#' @param attritionrates The attrition rates for each group, including the control group. It should be specified as a numeric vector of length ni + 1.
#'
#' @return A `data.frame` containing:
#' \describe{
#'   \item{ID}{Participant ID}
#'   \item{interventions}{Intervention assignment (0 = control, 1 to `ni` = intervention groups)}
#'   \item{pretest}{Pretest score}
#'   \item{posttest}{Posttest score (NA if participant attrited)}
#' }
#'
#' @examples
#' srtdata <- srtDataSimulation(ni = 2, np = 300, tpi = c(40, 30, 30),
#' sigma = 1, sigmaPret = 1, B0 = 0, B1 = 0.6, es = c(0.2, 0.3), 
#' seed = 101, attritionrates = c(0.1, 0.05, 0.05))
#' head(srtdata)
#'
#' @export
srtDataSimulation <- function(ni, np, tpi, sigma, sigmaPret, B0, B1, es, seed, attritionrates) {
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
  control_size <- floor(np * (tpi[1] / 100))
  control_pupils <- sample(1:np, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups to the remaining participants
  remaining_pupils <- setdiff(1:np, control_pupils)
  
  # Calculate the number of participants for each treatment group except the last one
  treatment_sizes <- floor(np * (tpi[2:(ni+1)] / 100))
  
  # Adjust the last treatment group to take the remaining participants
  treatment_sizes[ni] <- length(remaining_pupils) - sum(treatment_sizes[1:(ni-1)])
  
  for (i in 1:ni) {
    treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
    data[[interventions]][treated_pupils] <- i  # Assign treatment group i to selected participants
    remaining_pupils <- setdiff(remaining_pupils, treated_pupils)  # Remove assigned participants
  }
  
  # Step 4: Generate pre-test scores for everyone
  prets <- "pretest"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)
  
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
  # Create treatment matrix
  treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
  
  # Remove the first column (control group is baseline)
  if (ncol(treatment_matrix) > 1) {
    treatment_matrix <- treatment_matrix[, -1, drop = FALSE]
  }
  
  # Apply treatment effects using only individual-level variance (sigma)
  treatment_effects <- sweep(treatment_matrix, 2, es * sigma, `*`)
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects)
  
  # Generate individual-level errors for non-attrited individuals
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)
  
  # Compute post-test scores
  data[non_attrition_idx, posts] <- B0 + 
    B1 * data[non_attrition_idx, prets] + 
    total_treatment_effect + 
    e_ij
  
  # Only keep relevant columns (ID, interventions, pret, post)
  data <- data[, c("ID", interventions, prets, posts)]
  
  return(data)
}
