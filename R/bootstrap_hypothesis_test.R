#' Title
#'
#' @param my_data the data set you're using
#'
#' @return
#' @export
#'
#' @examples
Hypothesis_Test <- function(my_data = tortoise){

  my_data <- my_data %>%
    rowid_to_column("ID") %>%  # add ID column
    select(-type, -density)    # delete type and density columns

  # fit the model
  fit = run_model(my_data, example = "tortoise")

  # extract test statistc for prevalence
  t_value <- fit$test_stat[[2]]

  # generate bootstrap samples under null hypothesis
  B <- 1000

  bootstrap <- tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3),
           eta = fit$beta[[1]] + log(Area) + z,  # eta_ij* = beta_0 + log(xi2) + zi* (no beta_1 this time)
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(values = suppressMessages(run_model(data = data)), .groups = "drop")

  # compute p-value

  # first extract bootstrapped test statistics
  p_value <- bootstrap %>%
    select(values) %>%
    filter(row_number() %% 4 == 0) %>%
    unnest(cols = values) %>%
    filter(row_number() %% 4 == 2)

  p_value <- p_value %>%
    summarize(p = mean(abs(values) > abs(t_value)))

  # Return the p-value
  return(p_value$p)

}
