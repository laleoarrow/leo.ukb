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

make_add_interaction_df <- function() {
  set.seed(20260323)
  n <- 400
  exposure <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("No", "Yes"))
  sex <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Female", "Male"))
  exposure_num <- as.integer(exposure == "Yes")
  sex_num <- as.integer(sex == "Male")
  event_prob <- stats::plogis(-1.4 + 0.5 * exposure_num + 0.4 * sex_num + 0.7 * exposure_num * sex_num)
  data.frame(
    outcome = rbinom(n, 1, event_prob),
    outcome_censor = stats::rexp(n, rate = 0.08) + 0.1,
    exposure = exposure,
    sex = sex,
    age = stats::rnorm(n, 60, 8)
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
  expect_equal(names(res$formula_main), names(model))
  expect_equal(names(res$formula_inter), names(model))
  expect_true("Person-years" %in% names(res$result))
  expect_false("result_tidy" %in% names(res))
  expect_true(all(nzchar(res$result$`P for interaction`)))
})

test_that("leo_cox_interaction auto-detects continuous interaction variables", {
  lung_df <- make_lung_cox_df()

  res <- leo_cox_interaction(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "sex",
    x_inter = "age",
    x_exp_type = "categorical",
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_interaction")
  expect_equal(res$result$`Exposure class`, "Binary")
  expect_equal(res$result$`Interaction class`, "Continuous")
})

test_that("leo_cox_add_interaction returns additive interaction summaries for binary exposures", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_df()
  model <- list("Crude" = NULL, "Model A" = c("age"))

  res <- leo_cox_add_interaction(
    df = df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_inter = "sex",
    x_cov = model,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_add_interaction")
  expect_equal(res$result$Model, names(model))
  expect_equal(length(res$fit_main), length(model))
  expect_equal(length(res$fit_inter), length(model))
  expect_equal(names(res$formula_main), names(model))
  expect_equal(names(res$formula_inter), names(model))
  expect_type(res$backend, "list")
  expect_equal(names(res$backend), names(model))
  expect_true(all(c(
    "Reference group", "Recode applied",
    "Multiplicative interaction", "Multiplicative 95% CI", "P for interaction",
    "RERI", "RERI 95% CI", "AP", "AP 95% CI", "S", "S 95% CI"
  ) %in% names(res$result)))
})

test_that("leo_cox_add_interaction rejects non-binary interaction variables", {
  skip_if_not_installed("interactionR")
  lung_df <- make_lung_cox_df()

  expect_error(
    leo_cox_add_interaction(
      df = lung_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "sex",
      x_inter = "ecog_group",
      verbose = FALSE
    ),
    "must contain exactly 2 levels"
  )
})

test_that("leo_cox_add_interaction requires all four joint exposure groups", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_df()
  df <- df[!(df$exposure == "Yes" & df$sex == "Male"), , drop = FALSE]

  expect_error(
    leo_cox_add_interaction(
      df = df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_inter = "sex",
      verbose = FALSE
    ),
    "All four joint exposure groups"
  )
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
  expect_s3_class(res$interaction, "leo_cox_interaction")
  expect_true(all(res$result$`Case N` > 0))
  expect_true("Person-years" %in% names(res$result))
  expect_true("P for interaction" %in% names(res$result))
  expect_false("P for heterogeneity" %in% names(res$result))
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
  expect_true(all(c("sex", "ecog_group") %in% unique(res$result$Subgroup)))
  expect_type(res$interaction, "list")
  expect_true(all(c("sex", "ecog_group") %in% names(res$interaction)))
  expect_true(all(vapply(res$interaction, inherits, logical(1), what = "leo_cox_interaction")))
  expect_true(any(!is.na(res$result$`P for interaction`)))
})

test_that("leo_cox_subgroup can optionally append heterogeneity p values", {
  lung_df <- make_lung_cox_df()

  res <- leo_cox_subgroup(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "age",
    x_subgroup = "sex",
    x_cov = list("Crude" = NULL),
    add_heterogeneity = TRUE,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_subgroup")
  expect_true("P for heterogeneity" %in% names(res$result))
})

test_that("leo_cox_subgroup preserves original labels for continuous exposures", {
  lung_df <- make_lung_cox_df()
  model <- list("Crude" = NULL, "Model A" = c("ecog_group"))

  res <- leo_cox_subgroup(
    df = lung_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "age",
    x_subgroup = "sex",
    x_cov = model,
    verbose = FALSE
  )

  expect_true(all(res$result$Exposure == "age"))
  expect_true(all(res$result$Outcome == "outcome"))
  expect_true(all(res$result$`Exposure level` == "Per unit increase"))
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

test_that("leo_cox_mediation keeps readable labels and detailed outputs", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    verbose = FALSE
  )

  expect_true(all(res$result$`Exposure contrast` == "No -> Yes"))
  expect_true(all(res$result$`Mediator reference` == "High"))
  expect_true(all(c("result", "result_detail", "result_tidy", "evaluation", "fit") %in% names(res)))
  expect_identical(res$result_detail, res$result_tidy)
  expect_true(all(c("exposure_contrast", "mediator_reference", "c_cond") %in% names(res$evaluation)))
})

test_that("leo_cox_mediation accepts named c_cond for non-numeric covariates", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "smoking",
    c_cond = list(smoking = "Never"),
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_equal(res$evaluation$c_cond, "smoking=Never")
  expect_true(all(res$result$`Exposure contrast` == "No -> Yes"))
})

test_that("leo_cox_mediation rejects factor mediators when mediator_model is linear", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      mediator_model = "linear",
      verbose = FALSE
    ),
    "Mediator must be truly numeric"
  )
})
