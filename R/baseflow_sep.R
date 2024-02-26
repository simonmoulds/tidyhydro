

.baseflow_sep_lh <- function(strflow, ...) {
  # Initialization
  f1 <- 0.95
  f2 <- (1 + f1) / 2
  qs <- strflow
  qb <- matrix(rep(strflow, 3), ncol = 3, byrow = TRUE)
  qs[1] <- strflow[1] * 0.5
  qb[1, 1] <- strflow[1] - qs[1]
  qb[1, 2] <- qb[1, 1]
  qb[1, 3] <- qb[1, 1]

  # First pass (forward)
  for (i in 2:length(strflow)) {
    qs[i] <- f1 * qs[i-1] + f2 * (strflow[i] - strflow[i - 1])
    qs[i] <- max(0, qs[i])
    qb[i, 1] <- strflow[i] - qs[i]
    qb[i, 1] <- max(0, qb[i, 1])
    qb[i, 1] <- min(strflow[i], qb[i, 1])
  }

  # Second pass (backward)
  qb[length(strflow)-1, 2] <- qb[length(strflow) - 1, 1]
  for (i in (length(strflow)-2):1) {
    qs[i] <- f1 * qs[i + 1] + f2 * (qb[i, 1] - qb[i + 1, 1])
    qs[i] <- max(0, qs[i])

    qb[i, 2] <- qb[i, 1] - qs[i]
    qb[i, 2] <- max(0, qb[i, 2])
    qb[i, 2] <- min(qb[i, 1], qb[i, 2])
  }

  # Third pass (forward)
  qb[length(strflow) - 1, 3] <- qb[length(strflow) - 1, 1]
  for (i in 2:length(strflow)) {
    qs[i] <- f1 * qs[i - 1] + f2 * (qb[i, 2] - qb[i - 1, 2])
    qs[i] <- max(0, qs[i])

    qb[i, 3] <- qb[i, 2] - qs[i]
    qb[i, 3] <- max(0, qb[i, 3])
    qb[i, 3] <- min(qb[i, 2], qb[i, 3])
  }

  ## Take the third pass
  qb <- as.numeric(qb[,3])
  return(qb)
}


get_baseflow_sep_fun <- function(method = "LH", ...) {
  valid_methods <- c("LH")
  if (!method %in% valid_methods) {
    stop(paste0("Invalid method `", method, "`"))
  }
  if (method == "LH") {
    return(.baseflow_sep_lh)
  }
}


#' Baseflow separation
#'
#' @param x A tibble.
#' @param method A string specifying the method to use.
#' @param min_segment_length An integer specifying the
#'   minimum number of consecutive time points for which to
#'   perform baseflow separation.
#' @param ... Additional arguments passed to `baseflow_sep`
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A tibble
#' @export
#'
#' @examples
#' \dontrun{
#' print("Hello, world")
#' }
baseflow_sep <- function(x, method = "LH", min_segment_length = 30, ...) {

  baseflow_sep_fun <- get_baseflow_sep_fun(method)

  ## Make a copy of the original dataframe
  x_orig <- x

  ## Remove NA values from start/end of timeseries
  na_index <- which(complete.cases(x))
  x <- x[min(na_index):max(na_index), ]
  na_index <- which(complete.cases(x))

  ## The baseflow separation filter will not work if there are NA values in the
  ## timeseries. To get around this we divide the timeseries into continuous
  ## segments and separate the baseflow separately for each.

  ## Divide the timeseries into segments
  missing <- as.numeric(complete.cases(x))
  run_length <- rle(missing)
  vals <- run_length$values
  vals[vals > 0] <- seq_len(sum(vals))
  run_length$values <- vals
  segments <- inverse.rle(run_length)
  x <- x |>
    mutate(segment = segments) |>
    filter(.data$segment > 0)
    group_by(.data$segment) |>
    filter(n() >= min_segment_length) |>
    mutate(Qb = baseflow_sep_fun(.data$Q, ...)) |>
    mutate(Qs = .data$Q - .data$Qb) |>
    ungroup() |>
    dplyr::select(-all_of(c("Q", "segment")))

  x <- x_orig |> left_join(x, by = "date")
  return(x)
}
