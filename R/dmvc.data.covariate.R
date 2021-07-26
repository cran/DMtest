#' Example covariate data for DMVC
#'
#' Covariate data from 334 TCGA-COAD samples

#'
#' @docType data
#'
#' @usage data(covariate)
#'
#' @format An object of class \code{"matrix"} with with 334 rows and 3 vaiables.
#' \describe{
#'   \item{group}{Whether the sample is normal or tumor, normal:0, tumor:1}
#'   \item{gender}{Female or Male}
#'   \item{age}{age (31--90)}
#' }
#'
#' @keywords datasets
#'
#' @references Dai et al.
#'
#'
#' @examples
#' data(covariate)
"covariate"
