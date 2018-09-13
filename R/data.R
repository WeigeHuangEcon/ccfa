#' Intergenerational Mobility data from the PSID
#'
#' A dataset with 500 observations of matched parent's family income and child's family income data that also contains information on the education level of the family head (this is primary earner in the family)
#'
#' @format A data frame with 500 rows and 3 columns
#' \describe{
#'   \item{lcfincome}{log of child's family income}
#'   \item{lfincome}{log of parent's family income}
#'   \item{HEDUC}{head of family's education; less than HS, HS, or COLLEGE}
#' }
#' @source subset of PSID data used in Callaway and Huang (2017)
"igm"

