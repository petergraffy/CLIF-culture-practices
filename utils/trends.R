# Descriptive temporal models, fit separately at each site.
# Negative-binomial GAM for detection counts; quasi-binomial GAM for tested fractions.
# Smooth time and annual harmonics avoid forcing a constant log-linear trend.
fit_temporal_model <- function(data, denominator_var, count_var = "n_detection_events", proportion = FALSE, season_adjusted = TRUE) {
  empty <- function(status) list(summary = tibble::tibble(model_status = status, model_family = if (proportion) "quasibinomial_gam" else "negative_binomial_gam", log_annual_effect = NA_real_, log_annual_se = NA_real_, theta = NA_real_, long_term_edf = NA_real_, smoothing_parameter = NA_real_, monthly_irr = NA_real_, annual_irr = NA_real_, annual_percent_change = NA_real_, annual_ci_low = NA_real_, annual_ci_high = NA_real_, p_value = NA_real_, temporal_p_value = NA_real_, residual_lag1 = NA_real_, residual_dependence_flag = NA, n_model_months = 0L, model_warning = NA_character_), predictions = tibble::tibble(calendar_month = as.POSIXct(character(), tz = "UTC"), fitted_value = numeric(), ci_low = numeric(), ci_high = numeric()))
  if ("n_observed_culture_events" %in% names(data)) data <- dplyr::filter(data, n_observed_culture_events > 0)
  dat <- data %>% dplyr::filter(is.finite(.data[[denominator_var]]), .data[[denominator_var]] > 0, is.finite(.data[[count_var]])) %>%
    dplyr::arrange(calendar_month) %>% dplyr::mutate(y = .data[[count_var]], denominator = .data[[denominator_var]], time_years = as.numeric(difftime(calendar_month, min(calendar_month), units = "days")) / 365.25,
      season_sin = sin(2 * pi * (lubridate::month(calendar_month) - 1) / 12), season_cos = cos(2 * pi * (lubridate::month(calendar_month) - 1) / 12), log_denominator = log(denominator))
  if (nrow(dat) < 24 || (!proportion && sum(dat$y) < 25) || (proportion && sum(dat$denominator) < 30)) return(empty("insufficient_data"))
  if (proportion && (sum(dat$y) < 10 || sum(dat$denominator - dat$y) < 10)) return(empty("insufficient_outcome_variation"))
  if (any(dat$y < 0) || (proportion && any(dat$y > dat$denominator))) stop("Invalid counts in temporal model.")
  k <- min(8L, max(4L, floor(nrow(dat) / 12)))
  s <- mgcv::s
  formula <- if (proportion) cbind(y, denominator - y) ~ s(time_years, k = k, bs = "cr") + season_sin + season_cos else y ~ s(time_years, k = k, bs = "cr") + season_sin + season_cos + offset(log_denominator)
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(mgcv::gam(formula, family = if (proportion) quasibinomial() else mgcv::nb(), method = "REML", data = dat), warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
  if (inherits(fit, "error")) return(empty(paste0("failed: ", conditionMessage(fit))))
  # mgcv can set the inner convergence flag even when outer optimization fails.
  outer_failed <- !is.null(fit$outer.info$conv) && !identical(fit$outer.info$conv, "full convergence")
  convergence_warning <- any(grepl("iteration limit|without full convergence|failed to converge|not converged|step failure", warnings, ignore.case = TRUE))
  if (!isTRUE(fit$converged) || outer_failed || convergence_warning) {
    rejected <- empty("not_converged")
    rejected$summary$model_warning <- if(length(warnings)) paste(unique(warnings),collapse="; ") else if(outer_failed) as.character(fit$outer.info$conv) else NA_character_
    return(rejected)
  }
  years <- max(dat$time_years) - min(dat$time_years)
  # Standardize both endpoints to the same seasonal setting. Offsets cancel.
  endpoints <- dat[c(1, nrow(dat)), ]; endpoints$season_sin <- 0; endpoints$season_cos <- 0; endpoints$log_denominator <- 0
  xp <- predict(fit, endpoints, type = "lpmatrix")
  contrast <- (xp[2, ] - xp[1, ]) / years
  beta <- sum(contrast * coef(fit)); se <- sqrt(as.numeric(contrast %*% vcov(fit, unconditional = TRUE) %*% contrast))
  effects <- exp(c(beta,beta-1.96*se,beta+1.96*se))
  if (!is.finite(beta) || !is.finite(se) || se<=0 || any(!is.finite(effects) | effects<=0)) {
    rejected <- empty("invalid_endpoint_uncertainty")
    rejected$summary$model_warning <- "Endpoint effect or uncertainty is non-finite or underflows; excluded from inference."
    return(rejected)
  }
  month_number <- lubridate::year(dat$calendar_month) * 12 + lubridate::month(dat$calendar_month)
  adjacent <- which(diff(month_number) == 1)
  residual <- as.numeric(residuals(fit, type = "pearson"))
  lag1 <- if (length(adjacent) >= 12 && sd(residual[adjacent]) > 0 && sd(residual[adjacent + 1]) > 0) cor(residual[adjacent], residual[adjacent + 1]) else NA_real_
  temporal_p <- unname(summary(fit)$s.table[1, ncol(summary(fit)$s.table)])
  pred <- dat; pred$log_denominator <- 0
  # Match endpoint contrasts: hold both annual harmonics at their cycle mean.
  # Keep seasonal terms in the fitted model and its residual diagnostics.
  if (season_adjusted) {pred$season_sin <- 0; pred$season_cos <- 0}
  pr <- predict(fit, pred, type = "link", se.fit = TRUE, unconditional = TRUE)
  pr$fit <- as.numeric(pr$fit); pr$se.fit <- as.numeric(pr$se.fit)
  transform <- if (proportion) plogis else function(x) exp(x) * 100
  predictions <- tibble::tibble(calendar_month = dat$calendar_month, fitted_value = transform(pr$fit), ci_low = transform(pr$fit - 1.96 * pr$se.fit), ci_high = transform(pr$fit + 1.96 * pr$se.fit))
  list(summary = tibble::tibble(model_status = "estimated", model_family = if (proportion) "quasibinomial_gam" else "negative_binomial_gam", log_annual_effect = beta, log_annual_se = se, theta = if (proportion) NA_real_ else fit$family$getTheta(TRUE), long_term_edf = unname(summary(fit)$s.table[1, "edf"]), smoothing_parameter = unname(fit$sp[1]), monthly_irr = NA_real_, annual_irr = exp(beta), annual_percent_change = 100 * expm1(beta), annual_ci_low = exp(beta - 1.96 * se), annual_ci_high = exp(beta + 1.96 * se), p_value = if (is.finite(se) && se > 0) 2 * pnorm(-abs(beta / se)) else NA_real_, temporal_p_value = temporal_p, residual_lag1 = lag1, residual_dependence_flag = if (is.finite(lag1)) abs(lag1) > 0.3 else NA, n_model_months = nrow(dat), model_warning = if (length(warnings)) paste(unique(warnings), collapse = "; ") else NA_character_), predictions = predictions)
}

classify_temporal_direction <- function(model_status, fdr_p_value, annual_percent_change, residual_dependence_flag) {
  dplyr::case_when(model_status != "estimated" ~ "Not estimated", is.na(residual_dependence_flag) ~ "Residual diagnostic unavailable: exploratory", residual_dependence_flag ~ "Residual dependence: exploratory", fdr_p_value < 0.05 & annual_percent_change > 0 ~ "Increasing", fdr_p_value < 0.05 & annual_percent_change < 0 ~ "Decreasing", TRUE ~ "No clear net change")
}
