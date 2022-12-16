#' Final Project Hypothesis Test
#'
#' @param my_data the data set you're using
#'
#' @return the p-value obtained after testing then null hypothesis
#' @export
#'
#' @examples
#' # load tortoise data
#' load("~/Gopher.Fit/data/tortoise.rda")
#'
#' # Compute p-value
#' Hypothesis_Test()
#'
Hypothesis_Test <- function(my_data = tortoise){

  # fit the model
  fit = run_model(my_data, example = "tortoise")

  # extract test statistic for prevalence
  t_value <- fit$test_stat[[2]]

  # generate bootstrap samples under null hypothesis
  B <- 1000

  bootstrap <- tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3), # resample the random effects from N(0, sigmasq)
           eta = fit$beta[[1]] + log(Area) + z,  # eta_ij* = beta_0 + log(xi2) + zi* (no beta_1 this time)
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(values = suppressMessages(run_model(data = data)), .groups = "drop")  # refit

  # Compute p-value

  # first extract bootstrapped test statistics
  p_value <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 4 == 0) %>%  # run_model has a list of 4 outputs, 0 to take the test statistics
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 2)  # 2 to take the test-stat for prev

  # compute p-value
  p_value <- p_value %>%
    summarize(p = mean(abs(values) > abs(t_value)))

  # Return the p-value
  return(p_value$p)

}
