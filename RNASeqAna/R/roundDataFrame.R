#' A roundDataFrame Function
#'
#' This function is to round numeric value in a dataframe
#' @param x a dataframe
#' @param digits the digits will be round
#' @keywords math 
#' @export
#' @examples roundDataFrame()
#' roundDataFrame()


roundDataFrame <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  return(x)
}
