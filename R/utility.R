#' The inverse %in% operator
#'
#' Matching operator (x not in y) oposed to the %in%-operator (x in y)
#'
#' @param x	 vector
#' @param y	 vector of same type as x
#'
#' @return A logical vector
#'
#' @export
#' @concept utility
#'
#' @examples
#' \dontrun{
#' Notx=x%!in%y
#' }
#'

'%!in%' <- function(x,y)!('%in%'(x,y))
