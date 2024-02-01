#' Dataset on Annual Homicide Rates
#'
#' A dataset of homicide rates across 9 states from 1994 to 2008. Missouri repealed its permit-to-purchase (PTP) law in 2007.
#'
#' @format A dataframe of 7 variables with 135 observations (9 states, 15 years).
#' @name ptpdata
#' @docType data
#' @details
#'     details for dataset2
#' \describe{
#'   \item{state}{the name of the state; states present include Arkansas, Illinois, Iowa, Kansas, Kentucky, Missouri, Nebraska, Oklahoma, and Tennessee.}
#'   \item{year}{the year of each observation. Years range from 1994 to 2008.}
#'   \item{deaths}{the number of homicide deaths in the corresponding state in the corresponding year.}
#'   \item{crude_rate}{the homicide rate in the corresponding state in the corresponding year.}
#'   \item{age_adj_rate}{the age-adjusted homicide rate in the corresponding state in the corresponding year.}
#'   \item{group}{whether each observation belongs to the "to-be-treated" group; 1 for Missouri and 0 for all other states.}
#'   \item{treat}{whether each observation is treated; 1 for Missouri in 2008 and 0 otherwise.}
#' }
#'
"ptpdata"