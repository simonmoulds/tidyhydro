
#' Compute the baseflow index
#'
#' @param x A tibble.
#' @param method A string specifying the method to use.
#' @param min_segment_length An integer specifying the
#'   minimum number of consecutive time points for which to
#'   perform baseflow separation.
#' @param ... Additional arguments passed to `baseflow_sep`
#'
#' @return A tibble
#' @export
#'
#' @examples
#' \dontrun{
#' print("Hello, world")
#' }
baseflow_index <- function(x, method = "LH", min_segment_length = 30, ...) {
  x <- x |> baseflow_sep(method, min_segment_length, ...)
  bfi <- sum(x$Qb, na.rm = TRUE) / sum(x$Q, na.rm = TRUE)
  return(bfi)
}
