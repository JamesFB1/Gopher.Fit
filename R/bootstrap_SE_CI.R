#' Final Project SE and CI
#'
#' @param my_data the data set you're using
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
SE_CI = function(my_data = tortoise){

  my_data = my_data %>%
    rowid_to_column("ID") %>%  # add ID column
    select(-type, -density)    # delete type and density columns

  # fit the model
  fit = run_model(my_data, example = "tortoise")

  # make a column for beta_1 * prev
  my_data = my_data %>%
    mutate(.fitted = fit$beta[[2]] * prev)

  # generate bootstrap samples
  B = 1000

  bootstrap = tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3),
           eta = fit$beta[[1]] + log(Area) + .fitted + z,  # eta_ij* = beta_0 + beta_1xij1 + log(xi2) + zi*
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(values = suppressMessages(run_model(data = data)), .groups = "drop")

  # extract bootstrapped estimates
  Intercept = bootstrap %>%
    select(values) %>%
    filter(row_number() %% 4 == 1) %>%
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 1)

  Prev = bootstrap %>%
    select(values) %>%
    filter(row_number() %% 4 == 1) %>%
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 2)

  Sigmasq = bootstrap %>%
    select(values) %>%
    filter(row_number() %% 4 == 2) %>%
    unnest(cols = values)

  # compute standard errors and confidence intervals
  Intercept = Intercept %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  Prev = Prev %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  Sigmasq = Sigmasq %>%
    summarize(SE = sd(values),
              Lower95 = quantile(values, 0.025),
              Upper95 = quantile(values, 0.975))

  # return values
  return(list(Intercept = Intercept,
              Prev = Prev,
              Sigmasq = Sigmasq))
}
