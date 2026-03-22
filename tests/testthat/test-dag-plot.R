testthat::test_that("leo_dag_write requires explicit output directories", {
  testthat::skip_if_not_installed("dagitty")
  testthat::skip_if_not_installed("ggdag")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("readr")

  dag_obj <- dagitty::dagitty("dag { X -> Y; Z -> X; Z -> Y }")
  x <- leo.ukb::leo_dag(dag_obj, exposure_name = "X", outcome_name = "Y", dag_name = "toy_dag")

  testthat::expect_error(leo.ukb::leo_dag_write(x), "`figure_dir` must be provided explicitly")
  testthat::expect_error(
    leo.ukb::leo_dag_write(x, figure_dir = file.path(tempdir(), "figure")),
    "`output_dir` must be provided explicitly"
  )
})

testthat::test_that("leo_dag_write writes files to user-provided directories", {
  testthat::skip_if_not_installed("dagitty")
  testthat::skip_if_not_installed("ggdag")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("readr")

  dag_obj <- dagitty::dagitty("dag { X -> Y; Z -> X; Z -> Y }")
  x <- leo.ukb::leo_dag(dag_obj, exposure_name = "X", outcome_name = "Y", dag_name = "toy_dag")
  out_dir <- file.path(tempdir(), paste0("leo_dag_test_", as.integer(Sys.time())))
  figure_dir <- file.path(out_dir, "figure")
  output_dir <- file.path(out_dir, "output")

  res <- leo.ukb::leo_dag_write(x, figure_dir = figure_dir, output_dir = output_dir)

  testthat::expect_true(file.exists(res$pdf))
  testthat::expect_true(file.exists(res$adjustment_sets))
  testthat::expect_true(file.exists(res$roles))
  testthat::expect_true(startsWith(res$pdf, figure_dir))
  testthat::expect_true(startsWith(res$adjustment_sets, output_dir))
  testthat::expect_true(startsWith(res$roles, output_dir))
})
