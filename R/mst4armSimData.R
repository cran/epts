#' Simulated 4-Arm Multisite Trial (MST) Data
#'
#' A simulated multisite trial dataset containing 10 schools and 1,000 pupils.
#' This is a 4-arm trial design with one control group and three intervention groups.
#'
#' @format A data frame with 1,000 rows and 7 variables:
#' \describe{
#'   \item{pupils}{Identifier for each pupil}
#'   \item{schools}{Identifier for each school}
#'   \item{interventions}{Treatment assignment coded as 0 for control and 1–3 for intervention groups}
#'   \item{pretest}{Pre-test scores}
#'   \item{gender}{Binary gender}
#'   \item{ethnicity}{Ethnicity (3-level categorical)}
#'   \item{posttest}{Post-test scores}
#' }
#' @source Simulated
"mst4armSimData"