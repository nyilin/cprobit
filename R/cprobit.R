#' Inpernal function: compute geometric mean of a positive variable
#' @param x A numeric vector.
geom_mean <- function(x) {
  exp(mean(log(x)))
}
#' Inpernal function: generate commonly used summary statistics for estimates.
#' @details Vectorised, as long as the length of the input match.
#' @param var Names of variables.
#' @param est Estimated regression coefficients.
#' @param se SE of estimates.
#' @param z_score Z score of estimates, i.e., \code{est / se}.
#' @param pval P-value of estimates.
#' @param value_null Null effects for estimates, either with length 1 or length
#'   of \code{est}. Default is 0.
#' @param ci_lower Lower bound of 95\% CI of estimates.
#' @param ci_upper Upper bound of 95\% CI of estimates.
#' @param prefix Prefix to the column names in the \code{data.frame} returned.
#' @param postfix Postfix to the column names in the \code{data.frame} returned.
compile_est <- function(var, est, se = NULL, z_score = NULL, pval = NULL,
                        value_null = 0, ci_lower = NULL, ci_upper = NULL,
                        prefix = NULL, postfix = NULL) {
  if (is.null(se) & !is.null(z_score)) {
    se <- est / z_score
  }
  if (is.null(ci_lower) & is.null(ci_upper)) {
    ci_lower <- est - 1.96 * se
    ci_upper <- est + 1.96 * se
  }
  if (is.null(pval)) {
    if (!is.null(z_score)) {
      pval <- 2 * pnorm(q = abs(z_score), mean = 0, sd = 1, lower.tail = FALSE)
    } else {
      pval <- 2 * pnorm(q = abs((est - value_null) / se), mean = 0, sd = 1,
                        lower.tail = FALSE)
    }
  }
  ret <- data.frame(var = var, est = est, se = se,
                    ci_lower = ci_lower, ci_upper = ci_upper, pval = pval)
  if (!is.null(prefix)) {
    names(ret) <- paste(prefix, names(ret), sep = "_")
  }
  if (!is.null(postfix)) {
    names(ret) <- paste(names(ret), postfix, sep = "_")
  }
  rownames(ret) <- var
  ret
}
#' Inpernal function: construct design matrix without the intercept term.
#'
#' @param lp Formula for the linear predictor part, as a string.
#' @param dat Data to construct the design matrix from.
#' @param remove_intercept Whether the first column should be removed. Default
#'   is \code{TRUE} (to remove the intercept term).
#'
#' @return Returns a list containing the constructed design matrix and the
#'   original variable names. In the column names of the design matrix returned
#'   , any \code{:} in variable names are replaced with \code{.} to avoid
#'   computational issues when using the design matrix to fit model.
make_design_mat <- function(lp, dat, remove_intercept = TRUE) {
  if (remove_intercept) {
    design_mat <- model.matrix(as.formula(lp), data = dat)[, -1]
  } else {
    design_mat <- model.matrix(as.formula(lp), data = dat)
  }
  if (is.null(dim(design_mat))) {
    # Design matrix is actually a vector. Need to convert to a matrix with 1 col
    design_mat <- matrix(design_mat, ncol = 1)
    colnames(design_mat) <- as.character(as.formula(lp)[[2]])
  }
  nms <- sapply(colnames(design_mat), function(nm) {
    # Need to remove ":" from column names of interaction terms to avoid
    # problems when specifying 1st stage model.
    sub(pattern = ":", replacement = ".", x = nm)
  })
  nms_old <- colnames(design_mat)
  colnames(design_mat) <- nms
  list(design_mat = design_mat, var_names = nms_old)
}
#' Inpernal function: step 1 of the proposed workflow
#'
#' @description Implements the Step 1 of the proposed workflow, where a cprobit
#'   model is applied to analyse whether there is an increase in the outcome
#'   within each subject.
#' @param y_name Name of outcome variable for Step 1.
#' @param x_names Names of covariates for Step 1.
#' @param dat_diff A \code{data.frame} containing the difference data.
#' @param var_names Variable names for the estimates.
#' @return Returns a \code{data.frame} summarising the Step 1 estimates
#' (\code{coef}) and the covariance matrix for the Step 1 estimates (\code{vcov}).
cprobit_step1 <- function(y_name, x_names, dat_diff, var_names = NULL) {
  # Construct the formula for the cprobit model, which does not include the
  # intercept term.
  fml <- paste0(y_name, " ~ -1 + ", paste(x_names, collapse = " + "))
  m <- glm(as.formula(fml), family = binomial(link = "probit"), data = dat_diff)
  summ <- summary(m)$coefficients
  # Check if time-invariant confounders are included in formula
  if (nrow(summ) < length(var_names)) {
    stop(simpleError("Do not include time-invariant confounders in `formula` as they are implicitly controlled for in the cprobit model."))
  }
  if (is.null(var_names)) var_names <- rownames(summ)
  vcov <- vcov(m)
  rownames(vcov) <- var_names
  colnames(vcov) <- var_names
  list(coef = compile_est(var = var_names, est = summ[, 1], se = summ[, 2],
                          pval = summ[, 4]),
       vcov = vcov, x_names = x_names)
}
#' Inpernal function: compute difference in the (transformed) outcome
#'
#' @param y1 Numeric vector of the observed outcome at observation time 1.
#' @param y2 Numeric vector of the observed outcome at observation time 2.
#' @param lambda The Box-Cox transformation parameter. Default is \code{NA},
#'   indicating no need for a transformation. See \code{Details}.
#' @param scaled Whether the difference in the transformed outomes should be
#'   scaled by the Jacobian.
#' @return Returns the difference in the observed outcomes if \code{lambda =
#'   NA}, or the difference in the scaled transformed outcomes with
#'   transformation parameter \code{lambda}.
#' @importFrom car bcPower
get_v <- function(y1, y2, lambda = NA, scaled = TRUE) {
  if (is.na(lambda)) {
    y2 - y1
  } else {
    y_diff <- car::bcPower(U = y2, lambda = lambda, jacobian.adjusted = FALSE) -
      car::bcPower(U = y1, lambda = lambda, jacobian.adjusted = FALSE)
    if (scaled) {
      jacob <- geom_mean(c(y1, y2)) ^ (lambda - 1)
      y_diff / jacob
    } else {
      y_diff
    }
  }
}
#' Inpernal function: estimate the SD of error terms in the difference model
#'
#' @param beta_c Numeric vector of Step 1 estimates.
#' @param design_mat_diff Numeric matrix of the design matrix for difference.
#' @inheritParams get_v
#' @return Returns the estimate for \code{sigma_delta} if \code{lambda = NULL},
#'   or \code{sigma_delta_lambda} on the transformed scale.
estimate_sd_error <- function(beta_c, y1, y2, lambda = NA, design_mat_diff) {
  n <- length(y1)
  v <- get_v(y1 = y1, y2 = y2, lambda = lambda)
  if (ncol(design_mat_diff) == 1) {
    s <- sum(v * design_mat_diff[, 1] * beta_c)
  } else {
    s <- t(v) %*% design_mat_diff %*% beta_c
  }
  as.numeric((sqrt(s ^ 2 + 4 * n * sum(v ^ 2)) - s) / (2 * n))
}
#' Inpernal function: profile log-likelihood of lambda
#'
#' @inheritParams estimate_sd_error
#' @return Returns the profile log likelihood (not the negative value).
profile_llh <- function(lambda, beta_c, y1, y2, design_mat_diff) {
  v_tilde <- get_v(y1 = y1, y2 = y2, lambda = lambda, scaled = TRUE)
  sd_error_tilde <- estimate_sd_error(beta_c = beta_c, y1 = y1, y2 = y2,
                                      lambda = lambda,
                                      design_mat_diff = design_mat_diff)
  if (ncol(design_mat_diff) == 1) {
    e <- v_tilde - sd_error_tilde * design_mat_diff[, 1] * beta_c
  } else {
    e <- v_tilde - sd_error_tilde * design_mat_diff %*% beta_c
  }
  sum(dnorm(x = e, mean = 0, sd = sd_error_tilde, log = TRUE))
}
#' Inpernal function: update Step 1 estimates to obtain linear exposure effect on
#' (transformed) outcome
#'
#' @inheritParams cprobit_step1
#' @param y1_name Name of observed outcome at observation time 1.
#' @param y2_name Name of observed outcome at observation time 2.
#' @param res_step1 Results from Step 1 of the workflow.
#' @param transform Whether the outcome should be transformed. Default is
#'   \code{FALSE}.
#' @return Returns a \code{list}: a \code{data.frame} summarising the estimated
#'   linear exposure effect, the estimated standard deviation of the error terms
#'   from the difference model, the covariance matrix of the estimated exposure
#'   effects, a \code{data.frame} summarising the estimated transforamtion
#'   parameter, and the residuals.
#' @importFrom nortest lillie.test
update_estimate <- function(y1_name, y2_name, var_names = NULL,
                            dat_diff, res_step1, transform = FALSE) {
  y1 <- dat_diff[, y1_name]
  y2 <- dat_diff[, y2_name]
  design_mat_diff <- model.matrix(
    as.formula(paste("~ -1 +", paste(res_step1$x_names, collapse = " + "))),
    dat_diff
  )
  if (transform) {
    obj <- optimize(f = profile_llh, interval = c(-3, 3),
                    beta_c = res_step1$coef$est, y1 = y1, y2 = y2,
                    design_mat_diff = design_mat_diff, maximum = TRUE)
    hessian <- optimHess(par = obj$maximum, fn = profile_llh,
                         beta_c = res_step1$coef$est,
                         y1 = y1, y2 = y2, design_mat_diff = design_mat_diff)
    lambda <- obj$maximum
    lambda_se <- sqrt(solve(-hessian)[1, 1])
    y_tilde <- geom_mean(c(y1, y2))
    jacob <- y_tilde ^ (lambda - 1)
  } else {
    lambda <- NA
    lambda_se <- NA
    jacob <- 1
  }
  sd_error_tilde <- estimate_sd_error(beta_c = res_step1$coef$est,
                                      y1 = y1, y2 = y2, lambda = lambda,
                                      design_mat_diff = design_mat_diff)
  sd_error <- sd_error_tilde * jacob
  if (is.null(var_names)) var_names <- colnames(design_mat_diff)
  coef_mat <- compile_est(var = var_names, est = res_step1$coef$est * sd_error,
                          se = res_step1$coef$se * sd_error)
  v <- get_v(y1 = y1, y2 = y2, lambda = lambda, scaled = FALSE)
  if (ncol(design_mat_diff) == 1) {
    resid <- v - design_mat_diff[, 1] * coef_mat$est
  } else {
    resid <- v - design_mat_diff %*% coef_mat$est
  }
  list(coef = coef_mat, vcov = res_step1$vcov * (sd_error ^ 2),
       sd_error = sd_error,
       transformation = compile_est(var = "lambda", est = lambda,
                                    se = lambda_se, value_null = 1),
       resid = resid, resid_pval = nortest::lillie.test(resid)$p.value)
}
#' Apply the three-step workflow for the analysis of two repeated outcomes from
#' each subject
#'
#' @param formula Formula for the model. Do not convert data type within the
#'   formula (e.g., \code{factor(x)} is not supported in \code{formula}). See
#'   \code{Details}.
#' @param dat A \code{data.frame} in the long format, with each row
#'   corresponding to one measurement from one subject, and two columns
#'   indicating the subject and case ID respecitively. Variable names must not
#'   contain space or special characters.
#' @param index Names of variables indicating subject and case ID. Case ID must
#'   be coded as integers 1 and 2.
#' @param transform Whether a Box-Cox transformation should be applied to the
#'   outcome, taking value \code{NULL} (the default), \code{TRUE} or
#'   \code{FALSE}.
#' @param lambda Value of the Box-Cox transformation parameter to use. Default
#'   is \code{NA}, in which case it will be estimated from data.
#' @param resid_pval_threshold The threshold for the Lilliefors p-value of the
#'   residuals to determine whether a Box-Cox transformation on the outcome is
#'   necessary. Default is 0.05.
#'
#' @details Specify the formula for the repeated measurements instead of the
#'   change in the outcome, but without any time-invariant component that would
#'   have been eliminated after taking the difference. Interaction between two
#'   variables can be specified in the formula using \code{*} or \code{:}, but
#'   users need to create their own variable for interaction involving three or
#'   more variables.
#'
#' If \code{transform = NULL}, the workflow will determine the need for a
#' Box-Cox transforamtion on the outcome (i.e., Step 3) based on the residual
#' diagnostics in Step 2. A Box-Cox transforamtion will be used if the p-value
#' of the Lilliefors test is smaller than \code{resid_pval_threshold} (default
#' is 0.05). If \code{transform = TRUE}, analyses will always be performed on
#' both the observed and Box-Cox transformed outcomes. If \code{transform =
#' FALSE}, analysis will only be performed on the observed outcomes.
#'
#' @return Returns a list.
#'
#' @examples
#' # Apply the three-step workflow to assess the association between the
#' # baseline glucose variability and the change in the glucose variability in
#' # the subsequent two days.
#' # Although age and gender are available, they do not need to be explicitly
#' # adjusted for in the cprobit model.
#' data(bg_variability)
#' head(bg_variability)
#' model <- cprobit(formula = y ~ t + t:sd0, dat = bg_variability,
#'                  index = c("subject_id", "case_id"))
#' summary(model, plot = TRUE)
#' @references
#' \itemize{
#'  \item{GEP Box, DR Cox. An Analysis of Transformations. Journal of the Royal
#'  Statistical Society. Series B (Methodological). 1964;26:211–52.}
#'  \item{DM Hawkins, S Weisberg. Combining the box-cox power and generalised
#'  log transformations to accommodate nonpositive responses in linear and
#'  mixed-effects linear models. South African Stat J. 2017;51:317–28.}
#'  \item{HW Lilliefors. On the Kolmogorov-Smirnov Test for Normality with Mean
#'  and Variance Unknown. J Am Stat Assoc. 1967;62:399.}
#'  \item{Y Ning, NC Støer, PJ Ho, SL Kao, KY Ngiam, EYH Khoo, SC Lee, ES Tai, M
#'  Hartman, M Reilly, CS Tan. Robust estimation of the effect of an exposure on
#'  the change in a continuous outcome. BMC Medical Research Methodology
#'  (in press).}
#' }
#' @export
#' @importFrom car bcPower
cprobit <- function(formula, dat, index, transform = NULL, lambda = NA,
                    resid_pval_threshold = 0.05) {
  dat <- dat[order(dat[, index[1]], dat[, index[2]]), ]
  y_name <- as.character(formula[[2]])
  if ((is.null(transform) || transform) & any(dat[, y_name] <= 0)) {
    # `||` means evaluate the second term only when the first term is FALSE.
    stop(simpleError("Box-Cox transformation cannot be applied to variables with non-positive values."))
  }
  y1 <- dat[dat[, index[2]] == 1, y_name]
  y2 <- dat[dat[, index[2]] == 2, y_name]
  lp <- as.character(formula[[3]])
  if (length(lp) > 1) {
    lp <- paste(lp[-1], collapse = " + ")
  }
  lp <- paste("~", lp)
  design_mat_list <- make_design_mat(lp = lp, dat = dat)
  design_mat <- design_mat_list$design_mat
  x_names <- colnames(design_mat)
  design_mat_diff <- matrix(design_mat[dat[, index[2]] == 2, ] -
                              design_mat[dat[, index[2]] == 1, ],
                            ncol = ncol(design_mat))
  colnames(design_mat_diff) <- x_names
  dat_diff <- as.data.frame(cbind(I_y = y2 > y1, y1 = y1, y2 = y2,
                                  design_mat_diff))
  # Step 1
  res_step1 <- cprobit_step1(y_name = "I_y", x_names = x_names,
                             dat_diff = dat_diff,
                             var_names = design_mat_list$var_names)
  # Step 2
  res_step2 <- update_estimate(y1_name = "y1", y2_name = "y2",
                               var_names = design_mat_list$var_names,
                               dat_diff = dat_diff, res_step1 = res_step1,
                               transform = FALSE)
  res_step2$resid_pval_threshold <- resid_pval_threshold
  res_step3 <- NULL
  if (is.null(transform)) {
    # Follow the workflow: transform only when necessary
    if (res_step2$resid_pval < resid_pval_threshold) {
      res_step3 <- update_estimate(y1_name = "y1", y2_name = "y2",
                                   var_names = design_mat_list$var_names,
                                   dat_diff = dat_diff, res_step1 = res_step1,
                                   transform = TRUE)
    }
  } else if (transform) {
    # Transform anyway
    if (is.na(lambda)) {
      # Estimate the transformation parameter from data
      res_step3 <- update_estimate(y1_name = "y1", y2_name = "y2",
                                   var_names = design_mat_list$var_names,
                                   dat_diff = dat_diff, res_step1 = res_step1,
                                   transform = TRUE)
    } else {
      # With the transformation parameter known
      dat_diff$y1_lambda <- car::bcPower(U = dat_diff$y1, lambda = lambda)
      dat_diff$y2_lambda <- car::bcPower(U = dat_diff$y2, lambda = lambda)
      res_step3 <- update_estimate(y1_name = "y1_lambda", y2_name = "y2_lambda",
                                   var_names = design_mat_list$var_names,
                                   dat_diff = dat_diff, res_step1 = res_step1,
                                   transform = FALSE)
      res_step3$transformation$est[1] <- lambda
    }
  }
  ret <- list(Step1 = res_step1, Step2 = res_step2, Step3 = res_step3)
  class(ret) <- "cprobit"
  ret
}
#' @rdname cprobit
#' @param object Model fitted using \code{cprobit} function.
#' @param plot Wether residual qq-plots should be plotted. Default is \code{FALSE}.
#' @param ... Additional arguments affecting the summary produced (not yet
#'   implemented).
#' @import ggplot2, gridExtra
#' @export
summary.cprobit <- function(object, plot = FALSE, ...) {
  stopifnot(inherits(object, "cprobit"))
  cat("## Results from Step 2:\n\n")
  cat("Lilliefors test p-value for normality assumption\n",
      sprintf("    without transforamtion: %.3f", object$Step2$resid_pval),
      ifelse(object$Step2$resid_pval < object$Step2$resid_pval_threshold,
             "<", ">="),
      sprintf("%.3f", object$Step2$resid_pval_threshold), "\n\n")
  cat("Estimated coefficients for observed outcome:\n\n")
  coef_step2 <- object$Step2$coef
  coef_step2[, -1] <- apply(coef_step2[, -1], 2, function(x) round(x, 3))
  print(coef_step2)
  cat("\n")
  if (!is.null(object$Step3)) {
    cat("## Results from Step 3:\n\n")
    cat("Box-Cox transforamtion on the outcome:\n\n")
    trans_step3 <- object$Step3$transformation
    trans_step3[, -1] <- apply(trans_step3[, -1], 2, function(x) round(x, 3))
    print(trans_step3)
    cat("\n")
    cat("Lilliefors test p-value for normality assumption\n",
        sprintf("    after transformation: %.3f", object$Step3$resid_pval),
        ifelse(object$Step3$resid_pval < object$Step2$resid_pval_threshold,
               "<", ">="),
        sprintf("%.3f", object$Step2$resid_pval_threshold), "\n\n")
    cat("Estimated coefficients for transformed outcome:\n\n")
    coef_step3 <- object$Step3$coef
    coef_step3[, -1] <- apply(coef_step3[, -1], 2, function(x) round(x, 3))
    print(coef_step3)
    if (plot) {
      qqnorm(object$Step3$resid, main = "Normal qq-plot for transformed outcome")
      qqline(object$Step3$resid)
    }
  }
  # Make qq-plots
  if (plot) {
    common_theme <- theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    # Use `aes_string` instead of `aes` to avoid issue of "global variable y not
    # found"
    qq_step2 <- ggplot(data = data.frame(y = object$Step2$resid),
                       aes_string(sample = "y")) +
      stat_qq() +
      stat_qq_line() +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
           title = "Normal qq-plot for observed outcome") +
      common_theme
    if (is.null(object$Step3)) {
      qq_step2
    } else {
      qq_step3 <- ggplot(data = data.frame(y = object$Step3$resid),
                         aes_string(sample = "y")) +
        stat_qq() +
        stat_qq_line() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
             title = "Normal qq-plot for transformed outcome") +
        common_theme
      grid.arrange(qq_step2, qq_step3, ncol = 2)
    }
  }
}
#' @rdname cprobit
#' @param x Model fitted using \code{cprobit} function.
#' @export
print.cprobit <- function(x, ...) summary.cprobit(object = x, plot = FALSE, ...)

#' Inpatient blood glucose data for 1200 patients
#'
#' @description A simulated dataset containing the variability of inpatient
#'   point-of-care blood glucose (BG) measurements from 1200 non-critical care
#'   adult patients in medical ward. BG variability is measured as the standard
#'   deviation of the BG readings within a day. Data was simulated based on real
#'   data.
#'
#' @format A data frame with 1200 rows and 7 variables:
#' \describe{
#'   \item{subject_id}{Subject ID of each patient.}
#'   \item{case_id}{Case ID, with \code{1} and \code{2} referring to the first
#'   and second follow-up respectively.}
#'   \item{y}{BG variability of the first and second follow-up.}
#'   \item{t}{Binary indicator for the second follow-up.}
#'   \item{sd0}{Baseline BG variability.}
#'   \item{age}{Patients' age.}
#'   \item{female}{Binary indicator for being female.}
#' }
"bg_variability"
