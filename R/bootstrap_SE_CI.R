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
SE_CI = function(my_data, example, fitted_values, m){

  # Fit the model
  fit <- run_model(my_data, example)

  # Add fitted values to original data
  my_data <- my_data %>%
    mutate(eta_old = fitted_values)

  # generate bootstrap samples
  B = 1050  # generate 1050 bootstrap samples since model may not converge

  bootstrap <- tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/m, mean = 0, sd = sqrt(fit$sigmasq)), each = m),  # resample the random effects from N(0, sigmasq)
           eta_new = eta_old + z,  # eta_ij* = beta_0 + beta_1xij1 + log(xi2) + zi*
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(values = suppressWarnings(suppressMessages(run_model(data = data))), .groups = "drop")

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

  # extract bootstrapped betas
  beta = bootstrap %>%
    select(values) %>%
    filter(row_number() %% 6 == 1) %>%
    unnest_wider(values)

  # extract bootstrapped sigmasquared
  sigmasq <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 6 == 2) %>%  # run_model has a list of 6 outputs, 2 to take the sigmasq value
    unnest(values)

  # rename column
  colnames(sigmasq)[1] = "sigmasq"

  # combine data
  total = beta %>%
    mutate(sigmasq)

  # compute and return SE's and CI's
  tibble(Variable = colnames(total),
       SE = apply(total, 2, sd),
       Lower95 = apply(total, 2, quantile, prob = 0.025),
       Upper95 = apply(total, 2, quantile, prob = 0.975))

}
