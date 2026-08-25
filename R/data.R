#' Daily Hot-Hour Counts at KCWC Weather Station
#' 
#' Daily aggregated weather measurements for KCWC station. The main response is the
#' number of hot hours in each day, defined as the number of hourly temperature measurements
#' exceeding 95 degrees Fahrenheit.
#'
#' @format A data frame with 2665 rows and 9 variables.
#' \describe{
#'   \item{date}{Calendar date of the daily observation.}
#'   \item{hot}{Number of hot hours observed during the day.}
#'   \item{humidity}{Daily average humidity.}
#'   \item{t}{Time index from 1 to 2665.}
#'   \item{year}{Year extracted from \code{date}.}
#'   \item{month}{Month extracted from \code{date}.}
#'   \item{dow}{Day of the week extracted from \code{date}.}
#'   \item{sin365}{Annual sine term, computed as \eqn{\sin(2\pi t/365.25)}.}
#'   \item{cos365}{Annual cosine term, computed as \eqn{\cos(2\pi t/365.25)}.}
#' }
#' @usage data("KCWC")
#' @docType data
#' @keywords datasets
#' @name KCWC
NULL

#' Daily Cold-Hour Counts at KCMR Weather Station
#' 
#' Daily aggregated weather measurements for KCMR station. The main response is the
#' number of cold hours in each day, defined as the number of hourly temperature measurements
#' below 28 degrees Fahrenheit.
#'
#' @format A data frame with 2677 rows and 9 variables.
#' \describe{
#'   \item{date}{Calendar date of the daily observation.}
#'   \item{cold}{Number of cold hours observed during the day.}
#'   \item{humidity}{Daily average humidity.}
#'   \item{t}{Time index from 1 to 2677.}
#'   \item{year}{Year extracted from \code{date}.}
#'   \item{month}{Month extracted from \code{date}.}
#'   \item{dow}{Day of the week extracted from \code{date}.}
#'   \item{sin365}{Annual sine term, computed as \eqn{\sin(2\pi t/365.25)}.}
#'   \item{cos365}{Annual cosine term, computed as \eqn{\cos(2\pi t/365.25)}.}
#' }
#' @usage data("KCMR")
#' @docType data
#' @keywords datasets
#' @name KCMR
NULL

#' Weekly Campylobacter Case Counts Across Germany from 2001 to 2024
#' 
#' Weekly Campylobacter case counts reported across German districts from 2001 to 2024. Each year contains 52 weeks.
#' Each row corresponds to one reporting week, and the district-level columns contain
#' weekly case counts for individual districts. The data use 52 reporting weeks
#' per year.
#'
#' @format A data frame with 1248 rows and 411 variables.
#' @usage data("campyl")
#' @docType data
#' @keywords datasets
#' @name campyl
NULL

#' Weekly Rotavirus Case Counts Across Germany from 2001 to 2025
#' 
#' Weekly Rotavirus case counts reported across German districts from 2001 to 2025. Each year contains 52 weeks.
#' Each row corresponds to one reporting week, and the district-level columns contain
#' weekly case counts for individual districts. The data use 52 reporting weeks
#' per year.
#'
#' @format A data frame with 1300 rows and 411 variables.
#' @usage data("rota")
#' @docType data
#' @keywords datasets
#' @name rota
NULL