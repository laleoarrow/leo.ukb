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

make_mediation_character_df <- function() {
  set.seed(20260331)
  n <- 320
  age <- rnorm(n, 60, 8)
  exposure_chr <- rep(c("Unexposed", "Exposed"), length.out = n)
  exposure_num <- as.integer(exposure_chr == "Exposed")
  mediator_chr <- ifelse(
    rbinom(n, 1, plogis(-0.6 + 1.0 * exposure_num + 0.02 * (age - 60))) == 1,
    "High",
    "Low"
  )
  hazard_rate <- exp(-5.7 + 0.8 * exposure_num + 0.7 * (mediator_chr == "High") + 0.015 * (age - 60))
  event_time <- rexp(n, rate = hazard_rate)
  censor_time <- rexp(n, rate = 0.08)
  data.frame(
    outcome = as.integer(event_time <= censor_time),
    outcome_censor = pmax(pmin(event_time, censor_time), 0.1),
    exposure = exposure_chr,
    mediator = mediator_chr,
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
    simplify = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox")
  expect_equal(unique(res$result_tidy$model), names(model))
  expect_true(all(c("Crude HR", "Model A HR", "Model B HR") %in% names(res$result)))
  expect_true("Person-years" %in% names(res$result))
  expect_false("Person-time" %in% names(res$result))
  expect_equal(res$data_info$n_after_followup, res$data_info$n_after_complete_case)
})

test_that("leo_cox verbose messages include analysis spec and numbered model progress", {
  lung_df <- make_lung_cox_df()
  model <- list("Crude" = NULL, "Model A" = c("sex"))

  messages <- testthat::capture_messages(
    leo_cox(
      df = lung_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "age",
      x_cov = model,
      simplify = FALSE,
      verbose = TRUE
    )
  )

  joined <- paste(messages, collapse = "\n")
  expect_match(joined, "Analysis spec: exposure = age, outcome = outcome, model\\(s\\) = Crude, Model A")
  expect_match(joined, "Model \\[1/2\\]: Crude")
  expect_match(joined, "Model \\[2/2\\]: Model A")
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

test_that(".prepare_cox_regression_data parses factor follow-up times by their displayed numeric values", {
  df <- data.frame(
    outcome = c(1, 0, 1, 0),
    outcome_censor = factor(c("10", "20", "30", "40")),
    exposure = c(0.3, -0.2, 0.8, -0.4)
  )

  prep <- leo.ukb:::.prepare_cox_regression_data(
    df = df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    verbose = FALSE
  )

  expect_equal(prep$data$time, c(10, 20, 30, 40))
})

test_that(".prepare_cox_regression_data maps logical events to TRUE/FALSE by default and rejects unmatched factor labels", {
  logical_df <- data.frame(
    outcome = c(TRUE, FALSE, TRUE, FALSE),
    outcome_censor = c(5, 7, 9, 11),
    exposure = c(0.2, -0.1, 0.4, -0.3)
  )

  logical_prep <- leo.ukb:::.prepare_cox_regression_data(
    df = logical_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    verbose = FALSE
  )

  expect_equal(logical_prep$data$event, c(1L, 0L, 1L, 0L))

  factor_df <- data.frame(
    outcome = factor(c("Case", "Control", "Case", "Control")),
    outcome_censor = c(5, 7, 9, 11),
    exposure = c(0.2, -0.1, 0.4, -0.3)
  )

  expect_error(
    leo.ukb:::.prepare_cox_regression_data(
      df = factor_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      verbose = FALSE
    ),
    "Please set event_value explicitly"
  )

  factor_prep <- leo.ukb:::.prepare_cox_regression_data(
    df = factor_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    event_value = "Case",
    verbose = FALSE
  )

  expect_equal(factor_prep$data$event, c(1L, 0L, 1L, 0L))
})

test_that(".prepare_cox_regression_data rejects unmatched numeric event_value instead of silently creating all-zero events", {
  numeric_df <- data.frame(
    outcome = c(1, 0, 1, 0),
    outcome_censor = c(5, 7, 9, 11),
    exposure = c(0.2, -0.1, 0.4, -0.3)
  )

  expect_error(
    leo.ukb:::.prepare_cox_regression_data(
      df = numeric_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      event_value = 2,
      verbose = FALSE
    ),
    "Please set event_value explicitly"
  )

  expect_error(
    leo.ukb:::.prepare_cox_regression_data(
      df = numeric_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      event_value = "Case",
      verbose = FALSE
    ),
    "Please set event_value explicitly"
  )
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

test_that("leo_cox_subgroup verbose messages include numbered subgroup and level progress", {
  lung_df <- make_lung_cox_df()

  messages <- testthat::capture_messages(
    leo_cox_subgroup(
      df = lung_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "age",
      x_subgroup = "sex",
      x_cov = list("Crude" = NULL),
      add_interaction = FALSE,
      verbose = TRUE
    )
  )

  joined <- paste(messages, collapse = "\n")
  expect_match(joined, "Analysis spec: exposure = age, outcome = outcome, subgroup\\(s\\) = sex, model\\(s\\) = Crude")
  expect_match(joined, "Subgroup \\[1/1\\]: sex")
  expect_match(joined, "Level \\[1/2\\] in sex")
  expect_match(joined, "Level \\[2/2\\] in sex")
})

test_that("leo_cox_mediation regmedint route returns the current structured outputs", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  model <- list("Model A" = c("age"), "Model B" = c("age", "bmi"))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = model,
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation_regmedint")
  expect_s3_class(res, "leo_cox_mediation")
  expect_identical(names(res), c("result", "result_detail", "result_tidy", "evaluation", "method"))
  expect_identical(res$method, "regmedint")
  expect_equal(sort(unique(res$result$Model)), sort(names(model)))
  expect_true("Total follow-up" %in% names(res$result))
  expect_equal(sort(unique(res$result$`Effect code`)), sort(c("cde", "pnde", "tnie", "tnde", "pnie", "te", "pm")))
  expect_true(all(res$result$`Mediator model` == "logistic"))
  expect_true(all(grepl("^age=", res$result$`Covariate profile`)))
  expect_false("fit" %in% names(res))
})

test_that("leo_cox_mediation regmedint simplify modes and evaluation labels match the current contract", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res_wide <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "regmedint",
    simplify = "wide",
    verbose = FALSE
  )
  res_tidy <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "regmedint",
    simplify = "tidy",
    verbose = FALSE
  )
  res_full <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res_wide, "data.frame")
  expect_s3_class(res_tidy, "data.frame")
  expect_s3_class(res_full, "leo_cox_mediation")
  expect_identical(res_full$result_detail, res_full$result_tidy)
  expect_true(all(res_full$result$`Exposure contrast` == "No -> Yes"))
  expect_identical(unique(stats::na.omit(res_full$result$`Mediator reference`)), "High")
  expect_true(all(c("exposure_contrast", "mediator_reference", "x_exp_a0", "x_exp_a1", "x_med_cde", "x_cov_cond") %in% names(res_full$evaluation)))
})

test_that("leo_cox_mediation regmedint uses model-specific complete cases for yreg selection", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_multimodel_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = list("Model A" = "age", "Model B" = c("age", "bmi")),
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  eval_a <- dplyr::filter(res$evaluation, model == "Model A")
  eval_b <- dplyr::filter(res$evaluation, model == "Model B")
  expect_gt(eval_a$n, eval_b$n)
  expect_identical(eval_a$yreg, "survCox")
  expect_identical(eval_b$yreg, "survAFT_weibull")
})

test_that("leo_cox_mediation regmedint requires explicit profiles for non-numeric covariates", {
  skip_if_not_installed("regmedint")

  expect_error(
    leo_cox_mediation(
      df = make_mediation_df(),
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "smoking",
      method = "regmedint",
      verbose = FALSE
    ),
    "Please supply x_cov_cond explicitly"
  )

  expect_error(
    leo_cox_mediation(
      df = make_mediation_numeric_binary_cov_df(),
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "sex_num",
      method = "regmedint",
      verbose = FALSE
    ),
    "binary numeric covariates"
  )

  res <- leo_cox_mediation(
    df = make_mediation_df(),
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "smoking",
    x_cov_cond = list(smoking = "Never"),
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_identical(res$evaluation$x_cov_cond[1], "smoking=Never")
  expect_true(all(res$result$`Covariate profile` == "smoking=Never"))
})

test_that("leo_cox_mediation regmedint accepts explicit evaluation values for factor and character inputs", {
  skip_if_not_installed("regmedint")

  res_factor <- leo_cox_mediation(
    df = make_mediation_rare_df(),
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    x_exp_a0 = "No",
    x_exp_a1 = "Yes",
    x_med_cde = "High",
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )
  res_character <- leo_cox_mediation(
    df = make_mediation_character_df(),
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    x_exp_a0 = "Unexposed",
    x_exp_a1 = "Exposed",
    x_med_cde = "High",
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_true(all(res_factor$result$`Exposure contrast` == "No -> Yes"))
  expect_identical(unique(stats::na.omit(res_factor$result$`Mediator reference`)), "High")
  expect_true(all(res_factor$result_detail$x_med_cde_internal == 1))
  expect_identical(res_character$evaluation$exposure_contrast[1], "Unexposed -> Exposed")
  expect_identical(res_character$evaluation$mediator_reference[1], "High")
  expect_true(all(res_character$result_detail$x_exp_a0 == 0))
  expect_true(all(res_character$result_detail$x_exp_a1 == 1))
})

test_that("leo_cox_mediation regmedint auto mode keeps score mediators linear and rejects categorical mediators", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  med_df$mediator_many <- sample(1:12, nrow(med_df), replace = TRUE)
  med_df$mediator_group <- factor(sample(c("Low", "Mid", "High"), nrow(med_df), replace = TRUE))

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator_many",
    x_cov = "age",
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_true(all(res$result$`Mediator model` == "linear"))
  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator_group",
      x_cov = "age",
      method = "regmedint",
      verbose = FALSE
    ),
    "multi-level categorical mediators"
  )
})

test_that("leo_cox_mediation mediation route returns average-effect outputs with no default named profile", {
  skip_if_not_installed("mediation")
  med_df <- make_mediation_df()

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "mediation",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "leo_cox_mediation_mediation")
  expect_s3_class(res, "leo_cox_mediation")
  expect_identical(names(res), c("result", "result_detail", "result_tidy", "evaluation", "method"))
  expect_identical(res$method, "mediation")
  expect_equal(sort(unique(res$result$`Effect code`)), sort(c("acme", "ade", "ate", "pm")))
  expect_identical(unique(res$evaluation$yreg), "survreg")
  expect_true(all(res$result$`Outcome model` == "survreg"))
  expect_true(all(is.na(res$evaluation$x_cov_cond)))
  expect_true(all(is.na(res$result$`Covariate profile`)))
  expect_false(isTRUE(res$evaluation$mediation_boot[1]))
  expect_identical(res$evaluation$mediation_sims[1], 1000L)
  expect_true(is.na(res$evaluation$mediation_boot_ci_type[1]))
  expect_false(isTRUE(res$evaluation$mediation_use_speed[1]))
  expect_false("fit" %in% names(res))
})

test_that("leo_cox_mediation mediation route accepts explicit profiles and inference controls", {
  skip_if_not_installed("mediation")
  med_df <- make_mediation_df()

  res_profile <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = c("age", "smoking"),
    x_cov_cond = list(age = 60, smoking = "Never"),
    method = "mediation",
    simplify = FALSE,
    verbose = FALSE
  )
  res_boot <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "mediation",
    mediation_boot = TRUE,
    mediation_sims = 40,
    mediation_boot_ci_type = "perc",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_identical(res_profile$evaluation$x_cov_cond[1], "age=60.000; smoking=Never")
  expect_true(all(res_profile$result$`Covariate profile` == "age=60.000; smoking=Never"))
  expect_true(isTRUE(res_boot$evaluation$mediation_boot[1]))
  expect_identical(res_boot$evaluation$mediation_sims[1], 40L)
  expect_identical(res_boot$evaluation$mediation_boot_ci_type[1], "perc")
  expect_false(isTRUE(res_boot$evaluation$mediation_use_speed[1]))
})

test_that("leo_cox_mediation keeps engine-specific interaction defaults and validation", {
  skip_if_not_installed("mediation")
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()

  res_default <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    simplify = FALSE,
    verbose = FALSE
  )
  res_reg <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "regmedint",
    simplify = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res_default, "leo_cox_mediation_mediation")
  expect_false(isTRUE(unique(res_default$evaluation$interaction)))
  expect_true(isTRUE(unique(res_reg$evaluation$interaction)))
  expect_warning(
    res_int <- leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      method = "mediation",
      interaction = TRUE,
      simplify = FALSE,
      verbose = FALSE
    ),
    "averaged ACME/ADE/ATE/PM only|treatment-specific"
  )
  expect_true(isTRUE(unique(res_int$evaluation$interaction)))
})

test_that("leo_cox_mediation rejects unsupported engine combinations under the current contract", {
  med_df <- make_mediation_df()

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      method = "foo",
      verbose = FALSE
    ),
    "method must be either 'regmedint' or 'mediation'"
  )

  skip_if_not_installed("mediation")
  skip_if_not_installed("regmedint")

  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      method = "mediation",
      yreg = "survCox",
      simplify = FALSE,
      verbose = FALSE
    ),
    "must be 'auto'|internally fits `survival::survreg"
  )
  expect_error(
    leo_cox_mediation(
      df = med_df,
      y_out = c("outcome", "outcome_censor"),
      x_exp = "exposure",
      x_med = "mediator",
      x_cov = "age",
      method = "regmedint",
      mediation_boot = TRUE,
      verbose = FALSE
    ),
    "only available when method = 'mediation'"
  )
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

test_that("Step 3.0 tutorials use event-time simulation and inline-derived interpretation", {
  vignette_en <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3.Rmd")
  vignette_zh <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3-zh.Rmd")
  if (!file.exists(vignette_en) || !file.exists(vignette_zh)) testthat::skip("Step 3.0 vignette sources are not installed in the built package.")
  text_en <- paste(readLines(vignette_en, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_zh <- paste(readLines(vignette_zh, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  expect_match(text_en, "hazard_rate <- exp\\(")
  expect_match(text_en, "event_time <- rexp\\(")
  expect_match(text_en, "outcome = as.integer\\(event_time <= censor_time\\)")
  expect_false(grepl("outcome = rbinom", text_en, fixed = TRUE))
  expect_false(grepl("te = 1.234", text_en, fixed = TRUE))
  expect_false(grepl("tnie = 1.019", text_en, fixed = TRUE))
  expect_false(grepl("9.7%", text_en, fixed = TRUE))

  expect_match(text_zh, "hazard_rate <- exp\\(")
  expect_match(text_zh, "event_time <- rexp\\(")
  expect_match(text_zh, "outcome = as.integer\\(event_time <= censor_time\\)")
  expect_false(grepl("outcome = rbinom", text_zh, fixed = TRUE))
  expect_false(grepl("te = 1.234", text_zh, fixed = TRUE))
  expect_false(grepl("tnie = 1.019", text_zh, fixed = TRUE))
  expect_false(grepl("9.7%", text_zh, fixed = TRUE))

  expect_true(grepl("Compare the built-in journal palettes", text_en, fixed = TRUE))
  expect_true(grepl("内置期刊风格配色示例", text_zh, fixed = TRUE))
  for (pal in c("jama", "jco", "lancet", "nejm")) {
    expect_true(grepl(paste0("\"", pal, "\""), text_en, fixed = TRUE))
    expect_true(grepl(paste0("\"", pal, "\""), text_zh, fixed = TRUE))
  }
})

test_that("Step tutorials track current output contracts and reference labeling", {
  step1_en <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step1.Rmd")
  step1_zh <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step1-zh.Rmd")
  step2_en <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step2.Rmd")
  step2_zh <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step2-zh.Rmd")
  step3_en <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3.Rmd")
  step3_zh <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3-zh.Rmd")
  step31_en <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3-1.Rmd")
  step31_zh <- testthat::test_path("..", "..", "vignettes", "survival-analysis-step3-1-zh.Rmd")
  paths <- c(step1_en, step1_zh, step2_en, step2_zh, step3_en, step3_zh, step31_en, step31_zh)
  if (!all(file.exists(paths))) testthat::skip("One or more tutorial sources are not installed in the built package.")

  text_step1_en <- paste(readLines(step1_en, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step1_zh <- paste(readLines(step1_zh, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step2_en <- paste(readLines(step2_en, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step2_zh <- paste(readLines(step2_zh, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step3_en <- paste(readLines(step3_en, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step3_zh <- paste(readLines(step3_zh, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step31_en <- paste(readLines(step31_en, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text_step31_zh <- paste(readLines(step31_zh, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  expect_true(grepl("ECOG0 \\(Ref\\)", text_step1_en))
  expect_true(grepl("ECOG0 \\(Ref\\)", text_step1_zh))
  expect_true(grepl("Hazard ratio (log scale)", text_step1_en, fixed = TRUE))
  expect_true(grepl("Hazard ratio（对数尺度）", text_step1_zh, fixed = TRUE))

  expect_true(grepl("add_heterogeneity = TRUE", text_step2_en, fixed = TRUE))
  expect_true(grepl("add_heterogeneity = TRUE", text_step2_zh, fixed = TRUE))
  expect_true(grepl("P for heterogeneity", text_step2_en, fixed = TRUE))
  expect_true(grepl("P for heterogeneity", text_step2_zh, fixed = TRUE))

  expect_true(grepl("x_cov_cond = NA", text_step3_en, fixed = TRUE))
  expect_true(grepl("x_cov_cond = NA", text_step3_zh, fixed = TRUE))
  expect_true(grepl("`Covariate profile`", text_step3_en, fixed = TRUE))
  expect_true(grepl("`Covariate profile`", text_step3_zh, fixed = TRUE))
  expect_true(grepl("average-effect", text_step3_en, fixed = TRUE))
  expect_true(grepl("average-effect", text_step3_zh, fixed = TRUE))
  expect_true(grepl("x_cov_cond = list(age = 60, smoking = \"Never\")", text_step3_en, fixed = TRUE))
  expect_true(grepl("x_cov_cond = list(age = 60, smoking = \"Never\")", text_step3_zh, fixed = TRUE))
  expect_false(grepl("wrapper can default to median values when it needs", text_step3_en, fixed = TRUE))
  expect_false(grepl("wrapper 可以用中位数自动构造 `x_cov_cond`", text_step3_zh, fixed = TRUE))

  expect_true(grepl("`Mediator reference`", text_step31_en, fixed = TRUE))
  expect_true(grepl("`Mediator reference`", text_step31_zh, fixed = TRUE))
  expect_true(grepl("`Outcome model`", text_step31_en, fixed = TRUE))
  expect_true(grepl("`Outcome model`", text_step31_zh, fixed = TRUE))
})

test_that("leo_cox_mediation_plot can add or suppress regmedint notes", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file, width = 8.8, height = 5.8)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = list("Model A" = "age", "Model B" = c("age", "bmi")),
    method = "regmedint",
    verbose = FALSE
  )
  res$evaluation$pm_unstable[res$evaluation$model == "Model B"] <- TRUE

  assign("captured_text", character(), envir = .GlobalEnv)
  suppressMessages(trace(
    what = graphics:::text.default,
    tracer = quote({
      assign(
        "captured_text",
        c(get("captured_text", envir = .GlobalEnv), paste(labels, collapse = "\n")),
        envir = .GlobalEnv
      )
    }),
    print = FALSE
  ))
  on.exit({
    suppressMessages(untrace(graphics:::text.default))
    rm(captured_text, envir = .GlobalEnv)
  }, add = TRUE)

  p_with_note <- leo_cox_mediation_plot(
    x = res,
    model = "Model B",
    exposure_label = "Metabolic\nrisk",
    mediator_label = "Inflammation",
    outcome_label = "Incident\noutcome",
    language = "en",
    palette = "jama"
  )
  captured_with_note <- get("captured_text", envir = .GlobalEnv)

  assign("captured_text", character(), envir = .GlobalEnv)
  p_without_note <- leo_cox_mediation_plot(
    x = res,
    model = "Model B",
    exposure_label = "Metabolic\nrisk",
    mediator_label = "Inflammation",
    outcome_label = "Incident\noutcome",
    language = "en",
    palette = "jama",
    add_note = FALSE
  )
  captured_without_note <- get("captured_text", envir = .GlobalEnv)

  expect_s3_class(p_with_note, "recordedplot")
  expect_s3_class(p_without_note, "recordedplot")
  expect_true(any(grepl("Exposure contrast: No -> Yes", captured_with_note, fixed = TRUE)))
  expect_true(any(grepl("CDE mediator reference: High", captured_with_note, fixed = TRUE)))
  expect_true(any(grepl("Covariate profile: age=", captured_with_note, fixed = TRUE)))
  expect_true(any(grepl("TE 95% CI includes the null; interpret cautiously", captured_with_note, fixed = TRUE)))
  expect_false(any(grepl("Exposure contrast: No -> Yes", captured_without_note, fixed = TRUE)))
  expect_false(any(grepl("CDE mediator reference: High", captured_without_note, fixed = TRUE)))
  expect_false(any(grepl("Covariate profile:", captured_without_note, fixed = TRUE)))
  expect_false(any(grepl("TE 95% CI includes the null; interpret cautiously", captured_without_note, fixed = TRUE)))
})

test_that("leo_cox_mediation_plot uses average-effect labels and only shows mediation profiles when explicit", {
  skip_if_not_installed("mediation")
  med_df <- make_mediation_df()
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file, width = 8.8, height = 5.8)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  res_default <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "mediation",
    simplify = FALSE,
    verbose = FALSE
  )
  res_profile <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = c("age", "smoking"),
    x_cov_cond = list(age = 60, smoking = "Never"),
    method = "mediation",
    simplify = FALSE,
    verbose = FALSE
  )

  assign("captured_text", character(), envir = .GlobalEnv)
  suppressMessages(trace(
    what = graphics:::text.default,
    tracer = quote({
      assign(
        "captured_text",
        c(get("captured_text", envir = .GlobalEnv), paste(labels, collapse = "\n")),
        envir = .GlobalEnv
      )
    }),
    print = FALSE
  ))
  on.exit({
    suppressMessages(untrace(graphics:::text.default))
    rm(captured_text, envir = .GlobalEnv)
  }, add = TRUE)

  p_default <- leo_cox_mediation_plot(
    x = res_default,
    exposure_label = "Metabolic\nrisk",
    mediator_label = "Inflammation",
    outcome_label = "Incident\noutcome",
    language = "en",
    palette = "jama"
  )
  captured_default <- get("captured_text", envir = .GlobalEnv)

  assign("captured_text", character(), envir = .GlobalEnv)
  p_profile <- leo_cox_mediation_plot(
    x = res_profile,
    exposure_label = "Metabolic\nrisk",
    mediator_label = "Inflammation",
    outcome_label = "Incident\noutcome",
    language = "en",
    palette = "jama"
  )
  captured_profile <- get("captured_text", envir = .GlobalEnv)

  expect_s3_class(p_default, "recordedplot")
  expect_s3_class(p_profile, "recordedplot")
  expect_true(any(grepl("Indirect effect \\(ACME", captured_default)))
  expect_true(any(grepl("Direct effect \\(ADE", captured_default)))
  expect_true(any(grepl("Total effect \\(ATE", captured_default)))
  expect_true(any(grepl("Exposure contrast: No -> Yes", captured_default, fixed = TRUE)))
  expect_false(any(grepl("Covariate profile:", captured_default, fixed = TRUE)))
  expect_false(any(grepl("CDE mediator reference", captured_default, fixed = TRUE)))
  expect_true(any(grepl("Covariate profile: age=60.000; smoking=Never", captured_profile, fixed = TRUE)))
})

test_that("leo_cox_mediation_plot validates add_note", {
  skip_if_not_installed("regmedint")
  med_df <- make_mediation_df()
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file, width = 8.8, height = 5.8)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  res <- leo_cox_mediation(
    df = med_df,
    y_out = c("outcome", "outcome_censor"),
    x_exp = "exposure",
    x_med = "mediator",
    x_cov = "age",
    method = "regmedint",
    verbose = FALSE
  )

  expect_error(
    leo_cox_mediation_plot(x = res, add_note = NA, language = "en", palette = "jama"),
    "add_note must be TRUE or FALSE"
  )
  expect_error(
    leo_cox_mediation_plot(x = res, add_note = c(TRUE, FALSE), language = "en", palette = "jama"),
    "add_note must be TRUE or FALSE"
  )
})
