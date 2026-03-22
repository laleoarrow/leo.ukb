make_lung_cox_df <- function() {
  stats::na.omit(
    dplyr::transmute(
      survival::lung,
      outcome = as.integer(status == 2),
      outcome_censor = time / 365.25,
      age = age,
      sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
      ecog_group = factor(ph.ecog, levels = 0:3, labels = c("ECOG0", "ECOG1", "ECOG2", "ECOG3"))
    )
  )
}

make_custom_event_df <- function() {
  set.seed(20260322)
  data.frame(
    outcome = c(2, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0),
    outcome_censor = c(5, 8, 12, 16, 7, 20, 15, 6, 18, 22, 9, 25, 4, 11, 19, 10, 14, 21, 13, 24, 17, 23, 26, 27),
    exposure = rnorm(24),
    sex = factor(rep(c("Male", "Female"), each = 12))
  )
}

make_mediation_df <- function() {
  set.seed(123)
  n <- 300
  data.frame(
    outcome = rbinom(n, 1, plogis(-0.4 + 0.7 * rbinom(n, 1, 0.5))),
    outcome_censor = rexp(n, rate = 0.08) + 0.1,
    exposure = factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("No", "Yes")),
    mediator = factor(rbinom(n, 1, 0.45), levels = c(0, 1), labels = c("Low", "High")),
    age = rnorm(n, 60, 8),
    bmi = rnorm(n, 26, 4),
    smoking = factor(sample(c("Never", "Ever"), n, replace = TRUE))
  )
}

test_that("leo_cox keeps named model columns and tidy models aligned", {
  lung_df <- make_lung_cox_df()
  model <- list("Crude" = NULL, "Model A" = c("sex"), "Model B" = c("sex", "ecog_group"))

  res <- leo_cox(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "age",
    x_cov = model,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox")
  expect_equal(unique(res$result_tidy$model), names(model))
  expect_true(all(c("Crude HR", "Model A HR", "Model B HR") %in% names(res$result)))
  expect_true("Person-years" %in% names(res$result))
  expect_false("Person-time" %in% names(res$result))
  expect_equal(res$data_info$n_after_followup, res$data_info$n_after_complete_case)
})

test_that("leo_cox_interaction returns one row per model with finite interaction p values", {
  lung_df <- make_lung_cox_df()
  model <- list("Crude" = NULL, "Model A" = c("ecog_group"))

  res <- leo_cox_interaction(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "age",
    x_inter = "sex",
    x_cov = model,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_interaction")
  expect_equal(res$result$Model, names(model))
  expect_equal(length(res$fit_main), length(model))
  expect_equal(length(res$fit_inter), length(model))
  expect_true("Person-years" %in% names(res$result))
  expect_true(all(is.finite(res$result_tidy$p_interaction)))
  expect_equal(length(unique(res$result_tidy$n)), 1L)
})

test_that("leo_cox_subgroup propagates custom event_value into nested subgroup fits", {
  df <- make_custom_event_df()

  expect_no_error(
    res <- leo_cox_subgroup(
      df = df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_subgroup = "sex",
      x_cov = list("Crude" = NULL),
      event_value = 2,
      verbose = FALSE
    )
  )

  expect_s3_class(res, "leo_cox_subgroup")
  expect_true(all(res$result_tidy$case_n > 0))
  expect_true("Person-years" %in% names(res$result))
  expect_true(all(c("Crude P for interaction", "Crude P for heterogeneity") %in% names(res$result)))
})

test_that("leo_cox_subgroup supports multiple subgroup variables and retains interaction output", {
  lung_df <- make_lung_cox_df()
  model <- list("Crude" = NULL, "Model A" = c("ecog_group"))

  res <- leo_cox_subgroup(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "age",
    x_subgroup = c("sex", "ecog_group"),
    x_cov = model,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_subgroup")
  expect_true(all(c("sex", "ecog_group") %in% unique(res$result_tidy$subgroup)))
  expect_true(nrow(res$interaction) >= 2)
  expect_true(any(!is.na(res$result_tidy$p_interaction)))
})

test_that("leo_cox_mediation fits binary mediator models across multiple covariate models", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  model <- list("Model A" = c("age"), "Model B" = c("age", "bmi"))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = model,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_equal(sort(unique(res$result$Model)), sort(names(model)))
  expect_equal(nrow(res$result), 14L)
  expect_true("Person-years" %in% names(res$result))
  expect_true(all(c("Effect code", "Effect", "Scale", "Exposure contrast", "Mediator reference") %in% names(res$result)))
  expect_true(all(res$result$`Mediator model` == "logistic"))
  expect_true(all(c("cde", "pnde", "tnie", "tnde", "pnie", "te", "pm") %in% res$result$`Effect code`))
  expect_true(all(c(
    "Controlled direct effect",
    "Pure natural direct effect",
    "Total natural indirect effect",
    "Total natural direct effect",
    "Pure natural indirect effect",
    "Total effect",
    "Proportion mediated"
  ) %in% res$result$Effect))
})

test_that("leo_cox_mediation requests explicit c_cond for non-numeric covariates", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "smoking",
      verbose = FALSE
    ),
    "Please supply c_cond explicitly"
  )
})
