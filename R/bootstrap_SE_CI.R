#' Final Project SE and CI
#'
#' @param my_data the data set you're using
#' @param exmaple the name of your example: will be one of either "culcita", "ctsib", "epilepsy", or "tortoise"
#' @param m either the number of repeated observations from n individuals (repeated measures study)
#' or the number of subjects within n groups (grouped effects data)
#'
#' @return a list of the SE and CI of each estimate
#' @export
#'
#' @examples
#' # load tortoise data
#' load("~/Gopher.Fit/data/tortoise.rda")
#'
#' # compute SE and CI
#' SE_CI(tortoise)
SE_CI = function(my_data, example, m){

  # Fit the model
  fit <- run_model(my_data, example)

  # Add fitted values to original data
  my_data <- my_data %>%
    mutate(.fitted = fit$beta[[2]] * prev)

  # generate bootstrap samples
  B = 1050  # generate 1050 bootstrap samples since model may not converge

  bootstrap <- tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = m),  # resample the random effects from N(0, sigmasq)
           eta = fit$beta[[1]] + .fitted + log(Area) + z,  # eta_ij* = beta_0 + beta_1xij1 + log(xi2) + zi*
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(values = suppressMessages(run_model(data = data)), .groups = "drop")

  # check convergence
  Convergence <- bootstrap %>%
    filter(row_number() %% 6 == 5) %>%  # run_model has a list of 6 outputs, 5 to take the convergence
    unnest(cols = values) %>%
    filter(values == -1) %>%  # -1 - algorithm did not converge
    select(B)  # the bootstrap run numbers where convergence didn't occur

  # delete rows with unconverged estimates and select the first 1000 bootstrapped samples
  bootstrap <- bootstrap %>%
    rows_delete(Convergence) %>%
    slice(1:6000)  # to take the first 1000 bootstraps

  # extract bootstrapped estimates
  Intercept <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 6 == 1) %>%  # run_model has a list of 6 outputs, 1 to take the beta values
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 1)  # 1 to take the estimate for the Intercept

  Prev <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 6 == 1) %>%  # run_model has a list of 6 outputs, 1 to take the beta values
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 2)  # 2 to take the estimate for prev

  Sigmasq <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 6 == 2) %>%  # run_model has a list of 6 outputs, 2 to take the sigmasq value
    unnest(cols = values)

  # compute standard errors and confidence intervals
  Intercept <- Intercept %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  Prev <- Prev %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  Sigmasq <- Sigmasq %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  # return values
  return(list(Intercept = Intercept,
              Prev = Prev,
              Sigmasq = Sigmasq))
}
