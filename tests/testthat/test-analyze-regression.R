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
  age <- rnorm(n, 60, 8)
  bmi <- rnorm(n, 26, 4)
  smoking <- factor(sample(c("Never", "Ever"), n, replace = TRUE))
  exposure_num <- rbinom(n, 1, 0.5)
  mediator_num <- rbinom(n, 1, plogis(-0.4 + 0.9 * exposure_num + 0.015 * (age - 60) + 0.35 * (smoking == "Ever")))
  hazard_rate <- exp(-2.7 + 0.45 * exposure_num + 0.55 * mediator_num + 0.015 * (age - 60) + 0.02 * (bmi - 26) + 0.20 * (smoking == "Ever"))
  event_time <- rexp(n, rate = hazard_rate)
  censor_time <- rexp(n, rate = 0.08)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmax(pmin(event_time, censor_time), 0.1),
    exposure = factor(exposure_num, levels = c(0, 1), labels = c("No", "Yes")),
    mediator = factor(mediator_num, levels = c(0, 1), labels = c("Low", "High")),
    age = age,
    bmi = bmi,
    smoking = smoking
  )
}

make_mediation_rare_df <- function() {
  set.seed(20260327)
  n <- 400
  age <- rnorm(n, 60, 8)
  exposure_num <- rbinom(n, 1, 0.5)
  mediator_num <- rbinom(n, 1, plogis(-1.0 + 1.0 * exposure_num + 0.01 * (age - 60)))
  hazard_rate <- exp(-6.0 + 0.9 * exposure_num + 1.0 * mediator_num + 0.01 * (age - 60))
  event_time <- rexp(n, rate = hazard_rate)
  censor_time <- rexp(n, rate = 0.08)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmax(pmin(event_time, censor_time), 0.1),
    exposure = factor(exposure_num, levels = c(0, 1), labels = c("No", "Yes")),
    mediator = factor(mediator_num, levels = c(0, 1), labels = c("Low", "High")),
    age = age
  )
}

make_mediation_multimodel_df <- function() {
  set.seed(20260328)
  n <- 200
  age <- rnorm(n, 60, 8)
  bmi <- rnorm(n, 26, 4)
  smoking <- factor(sample(c("Never", "Ever"), n, replace = TRUE))
  exposure_num <- rbinom(n, 1, 0.5)
  mediator_num <- rbinom(n, 1, plogis(-0.5 + 0.9 * exposure_num + 0.01 * (age - 60)))
  hazard_rate <- exp(-6.0 + 0.7 * exposure_num + 0.8 * mediator_num + 0.01 * (age - 60))
  event_time <- rexp(n, rate = hazard_rate)
  censor_time <- rexp(n, rate = 0.08)
  outcome <- as.integer(event_time <= censor_time)
  non_event_idx <- which(outcome == 0L)
  bmi[sample(non_event_idx, size = 120)] <- NA_real_
  data.frame(
    outcome = outcome,
    outcome_censor = pmax(pmin(event_time, censor_time), 0.1),
    exposure = factor(exposure_num, levels = c(0, 1), labels = c("No", "Yes")),
    mediator = factor(mediator_num, levels = c(0, 1), labels = c("Low", "High")),
    age = age,
    bmi = bmi,
    smoking = smoking
  )
}

make_mediation_numeric_binary_cov_df <- function() {
  set.seed(20260329)
  med_df <- make_mediation_df()
  med_df$sex_num <- sample(c(0, 1), nrow(med_df), replace = TRUE)
  med_df
}

make_mediation_pm_unstable_df <- function() {
  set.seed(20260330)
  n <- 500
  age <- rnorm(n, 60, 8)
  exposure_num <- rbinom(n, 1, 0.5)
  mediator_num <- rbinom(n, 1, plogis(-0.4 + 0.9 * exposure_num + 0.01 * (age - 60)))
  hazard_rate <- exp(-2.1 + 0.015 * (age - 60))
  event_time <- rexp(n, rate = hazard_rate)
  censor_time <- rexp(n, rate = 0.08)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmax(pmin(event_time, censor_time), 0.1),
    exposure = factor(exposure_num, levels = c(0, 1), labels = c("No", "Yes")),
    mediator = factor(mediator_num, levels = c(0, 1), labels = c("Low", "High")),
    age = age
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

make_add_interaction_preventive_df <- function() {
  set.seed(20260324)
  n <- 600
  exposure <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("No", "Yes"))
  sex <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Female", "Male"))
  exposure_num <- as.integer(exposure == "Yes")
  sex_num <- as.integer(sex == "Male")
  hazard_rate <- exp(-2.4 - 0.9 * exposure_num + 0.25 * sex_num - 0.1 * exposure_num * sex_num)
  event_time <- stats::rexp(n, rate = hazard_rate)
  censor_time <- stats::rexp(n, rate = 0.08)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmin(event_time, censor_time),
    exposure = exposure,
    sex = sex,
    age = stats::rnorm(n, 60, 8)
  )
}

make_add_interaction_joint_protective_df <- function() {
  set.seed(20260325)
  n <- 3000
  exposure <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("No", "Yes"))
  sex <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Female", "Male"))
  exposure_num <- as.integer(exposure == "Yes")
  sex_num <- as.integer(sex == "Male")
  hazard_rate <- exp(-4 + 0.8 * exposure_num + 0.7 * sex_num - 2.0 * exposure_num * sex_num)
  event_time <- stats::rexp(n, rate = hazard_rate)
  censor_time <- stats::rexp(n, rate = 0.03)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmin(event_time, censor_time),
    exposure = exposure,
    sex = sex
  )
}

make_add_interaction_nested_name_df <- function() {
  set.seed(20260326)
  n <- 800
  exp <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("No", "Yes"))
  exp2 <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Low", "High"))
  exp_num <- as.integer(exp == "Yes")
  exp2_num <- as.integer(exp2 == "High")
  hazard_rate <- exp(-2.8 - 0.8 * exp_num + 0.35 * exp2_num + 0.15 * exp_num * exp2_num)
  event_time <- stats::rexp(n, rate = hazard_rate)
  censor_time <- stats::rexp(n, rate = 0.05)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmin(event_time, censor_time),
    exp = exp,
    exp2 = exp2
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

test_that("leo_cox_add_interaction keeps the double-negative group as reference when requested", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_preventive_df()
  model <- list("Crude" = NULL, "Model A" = c("age"))

  res <- suppressWarnings(leo_cox_add_interaction(
    df = df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_inter = "sex",
    x_cov = model,
    reference_group = "double_negative",
    verbose = FALSE
  ))

  expect_equal(res$result$Model, names(model))
  expect_true(all(res$result$`Reference group` == "exposure=No, sex=Female"))
  expect_true(all(!res$result$`Recode applied`))
})

test_that("leo_cox_add_interaction auto mode follows the package preventive recode convention", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_preventive_df()

  res <- suppressWarnings(leo_cox_add_interaction(
    df = df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_inter = "sex",
    reference_group = "auto",
    verbose = FALSE
  ))

  expect_identical(res$result$`Recode applied`, TRUE)
  expect_identical(res$result$`Reference group`, "exposure=Yes, sex=Female")
})

test_that("leo_cox_add_interaction auto mode does not recode when only the joint group is protective", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_joint_protective_df()

  res <- suppressWarnings(leo_cox_add_interaction(
    df = df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_inter = "sex",
    reference_group = "auto",
    verbose = FALSE
  ))

  expect_identical(res$result$`Recode applied`, FALSE)
  expect_identical(res$result$`Reference group`, "exposure=No, sex=Female")
})

test_that("leo_cox_add_interaction auto mode matches coefficient names exactly", {
  skip_if_not_installed("interactionR")
  df <- make_add_interaction_nested_name_df()

  expect_no_error(
    res <- leo_cox_add_interaction(
      df = df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exp",
      x_inter = "exp2",
      reference_group = "auto",
      verbose = FALSE
    )
  )
  expect_s3_class(res, "leo_cox_add_interaction")
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
  expect_true("Total follow-up" %in% names(res$result))
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

test_that("leo_cox_mediation auto-selects survCox for rare event outcomes", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_rare_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    verbose = FALSE
  )

  expect_true(all(res$result$`Outcome model` == "survCox"))
  expect_identical(unique(res$evaluation$yreg), "survCox")
})

test_that("leo_cox_mediation rare-event helper yields a positive mediated signal on the Cox scale", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_rare_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    verbose = FALSE
  )

  te_row <- dplyr::filter(res$result, `Effect code` == "te")
  tnie_row <- dplyr::filter(res$result, `Effect code` == "tnie")
  expect_gt(te_row$Estimate[1], 1)
  expect_gt(tnie_row$Estimate[1], 1)
})

test_that("leo_cox_mediation auto-selects survAFT_weibull for non-rare event outcomes", {
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

  expect_true(all(res$result$`Outcome model` == "survAFT_weibull"))
  expect_identical(unique(res$evaluation$yreg), "survAFT_weibull")
})

test_that("leo_cox_mediation labels AFT mediation effects as time ratios", {
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

  non_pm_rows <- res$result$`Effect code` != "pm"
  expect_true(all(res$result$Scale[non_pm_rows] == "Time ratio"))
  expect_true(all(res$result$Scale[!non_pm_rows] == "Proportion"))
})

test_that("leo_cox_mediation flags unstable proportion mediated estimates when total effect crosses the null", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_pm_unstable_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    verbose = FALSE
  )

  expect_true(all(c("pm_unstable", "pm_note") %in% names(res$evaluation)))
  expect_true(all(res$evaluation$pm_unstable))
  expect_true(all(grepl("Total effect 95% CI includes the null value", res$evaluation$pm_note, fixed = TRUE)))
})

test_that("leo_cox_mediation still warns when survCox is forced for non-rare outcomes", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  expect_warning(
    res <- leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      yreg = "survCox",
      verbose = FALSE
    ),
    "rare-event approximation"
  )
  expect_true(all(res$result$`Outcome model` == "survCox"))
  expect_identical(unique(res$evaluation$yreg), "survCox")
})

test_that("leo_cox_mediation requests explicit x_cov_cond for non-numeric covariates", {
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
    "Please supply x_cov_cond explicitly"
  )
})

test_that("leo_cox_mediation requests explicit x_cov_cond for binary numeric covariates", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_numeric_binary_cov_df()

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "sex_num",
      verbose = FALSE
    ),
    "binary numeric covariates"
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
  expect_true(all(c("exposure_contrast", "mediator_reference", "x_exp_a0", "x_exp_a1", "x_med_cde", "x_cov_cond") %in% names(res$evaluation)))
})

test_that("leo_cox_mediation accepts named x_cov_cond for non-numeric covariates", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "smoking",
    x_cov_cond = list(smoking = "Never"),
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_equal(res$evaluation$x_cov_cond, "smoking=Never")
  expect_true(all(res$result$`Exposure contrast` == "No -> Yes"))
})

test_that("leo_cox_mediation accepts renamed public evaluation parameters", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_cont <- stats::rnorm(nrow(med_df))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator_cont",
    x_cov = c("age", "smoking"),
    x_med_cde = stats::median(med_df$mediator_cont, na.rm = TRUE),
    x_cov_cond = list(age = 60, smoking = "Never"),
    mreg = "linear",
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
})

test_that("leo_cox_mediation accepts factor-level inputs for x_exp_a0, x_exp_a1, and x_med_cde", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_rare_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    x_exp_a0 = "No",
    x_exp_a1 = "Yes",
    x_med_cde = "High",
    verbose = FALSE
  )

  expect_true(all(res$result$`Exposure contrast` == "No -> Yes"))
  expect_true(all(res$result$`Mediator reference` == "High"))
  expect_true(all(res$result_detail$x_exp_a0 == 0))
  expect_true(all(res$result_detail$x_exp_a1 == 1))
  expect_true(all(res$result_detail$x_med_cde == 1))
})

test_that("leo_cox_mediation renames public evaluation fields in evaluation and result_detail", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_cont <- stats::rnorm(nrow(med_df))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator_cont",
    x_cov = c("age", "smoking"),
    x_med_cde = stats::median(med_df$mediator_cont, na.rm = TRUE),
    x_cov_cond = list(age = 60, smoking = "Never"),
    mreg = "linear",
    verbose = FALSE
  )

  expect_true(all(c("x_exp_a0", "x_exp_a1", "x_med_cde", "x_cov_cond") %in% names(res$evaluation)))
  expect_true(all(c("x_exp_a0", "x_exp_a1", "x_med_cde", "x_cov_cond") %in% names(res$result_detail)))
  expect_false(any(c("a0", "a1", "m_cde", "c_cond") %in% names(res$evaluation)))
  expect_false(any(c("a0", "a1", "m_cde") %in% names(res$result_detail)))
  expect_false(any(c("control_n", "person_year") %in% names(res$result_detail)))
  expect_true("exposure_contrast_value" %in% names(res$evaluation))
  expect_false("mediator_reference_value" %in% names(res$evaluation))
  expect_false(any(c("exposure_contrast_label", "mediator_reference_label") %in% names(res$result_detail)))
})

test_that("leo_cox_mediation no longer accepts old evaluation parameter names", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_cont <- stats::rnorm(nrow(med_df))

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator_cont",
      x_cov = "age",
      a0 = 0,
      a1 = 1,
      m_cde = 0,
      c_cond = 60,
      mreg = "linear",
      verbose = FALSE
    ),
    "unused argument"
  )
})

test_that("leo_cox_mediation accepts mreg as the public mediator model parameter", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_cont <- stats::rnorm(nrow(med_df))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator_cont",
    x_cov = "age",
    mreg = "linear",
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_true(all(res$result$`Mediator model` == "linear"))
  expect_identical(unique(res$evaluation$mreg), "linear")
})

test_that("leo_cox_mediation no longer accepts mediator_model as a public parameter", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_cont <- stats::rnorm(nrow(med_df))

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator_cont",
      x_cov = "age",
      mediator_model = "linear",
      verbose = FALSE
    ),
    "unused argument"
  )
})

test_that("leo_cox_mediation rejects factor mediators when mreg is linear", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      mreg = "linear",
      verbose = FALSE
    ),
    "Mediator must be truly numeric"
  )
})

test_that("leo_cox_mediation auto mode treats integer-coded score mediators as linear", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_tri <- sample(c(0, 1, 2), nrow(med_df), replace = TRUE)

  expect_message(
    res <- leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator_tri",
      x_cov = "age",
      verbose = TRUE
    ),
    "ignore this warning if expected"
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_true(all(res$result$`Mediator model` == "linear"))
})

test_that("leo_cox_mediation auto mode treats broader integer score mediators as linear", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_many <- sample(1:12, nrow(med_df), replace = TRUE)

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator_many",
    x_cov = "age",
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation")
  expect_true(all(res$result$`Mediator model` == "linear"))
})

test_that("leo_cox_mediation auto-selects yreg and complete cases separately for each model", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_multimodel_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = list("Model A" = "age", "Model B" = c("age", "bmi")),
    verbose = FALSE
  )

  eval_a <- dplyr::filter(res$evaluation, model == "Model A")
  eval_b <- dplyr::filter(res$evaluation, model == "Model B")
  expect_gt(eval_a$n, eval_b$n)
  expect_identical(eval_a$yreg, "survCox")
  expect_identical(eval_b$yreg, "survAFT_weibull")
})

test_that("pkgdown build script checks the current Step 3.0 and Step 3.1 article titles", {
  script_path <- testthat::test_path("..", "..", "scripts", "build_pkgdown_site.R")
  if (!file.exists(script_path)) testthat::skip("pkgdown source script is not installed in the built package.")
  script_lines <- readLines(script_path, warn = FALSE, encoding = "UTF-8")
  script_text <- paste(script_lines, collapse = "\n")

  expect_match(script_text, "Survival Analysis Step 3\\.0: Mediation Analysis")
  expect_match(script_text, "Survival Analysis Step 3\\.1: Mediation Case Study")
  expect_match(script_text, "生存分析 Step 3\\.0：中介分析")
  expect_match(script_text, "生存分析 Step 3\\.1：中介分析案例教程")
  expect_match(script_text, "survival-analysis-step3-1\\.html")
  expect_match(script_text, "survival-analysis-step3-1-zh\\.html")
  expect_false(grepl("Survival Analysis Step 3: Mediation Analysis", script_text, fixed = TRUE))
  expect_false(grepl("生存分析 Step 3：中介分析", script_text, fixed = TRUE))
})

test_that("leo_cox_mediation auto mode still rejects multi-level categorical mediators", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_group <- factor(sample(c("Low", "Mid", "High"), nrow(med_df), replace = TRUE))

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator_group",
      x_cov = "age",
      verbose = FALSE
    ),
    "multi-level categorical mediators"
  )
})

test_that("leo_cox_mediation_plot returns a recorded plot for a fitted mediation result", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = list("Model A" = "age", "Model B" = c("age", "bmi")),
    verbose = FALSE
  )

  p_en <- leo_cox_mediation_plot(
    x = res,
    model = "Model B",
    exposure_label = "Metabolic\nrisk",
    mediator_label = "Inflammation",
    outcome_label = "Incident\noutcome",
    language = "en",
    palette = "jco"
  )
  p_zh <- leo_cox_mediation_plot(
    x = res,
    model = "Model B",
    exposure_label = "代谢风险",
    mediator_label = "炎症",
    outcome_label = "结局事件",
    language = "zh",
    palette = "jco"
  )

  expect_s3_class(p_en, "recordedplot")
  expect_s3_class(p_zh, "recordedplot")
})
