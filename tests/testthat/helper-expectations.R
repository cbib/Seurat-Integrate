expect_nan <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")

  # 2. Call expect()
  act$na <- is.na(act$val)
  expect(
    act$na,
    sprintf("%s is %s, not NA.", act$lab, act$val)
  )

  # 3. Invisibly return the value
  invisible(act$val)
}

expect_inf <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")

  # 2. Call expect()
  act$inf <- is.infinite(act$val)
  expect(
    act$inf,
    sprintf("%s is %s, not (-)Inf.", act$lab, act$val)
  )

  # 3. Invisibly return the value
  invisible(act$val)
}

expect_squared <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")

  # 2. Call expect()
  act$nrow <- nrow(act$val)
  act$ncol <- ncol(act$val)
  expect(
    act$nrow == act$ncol,
    sprintf("%s is not squared (%d columns and %d rows).", act$lab, act$nrow, act$ncol)
  )

  # 3. Invisibly return the value
  invisible(act$val)
}

