fit.summary = function(my_data, example, fitted_values, m, response, beta_i = 0){
  fit = run_model(my_data, example)
  stand = SE_CI(my_data, example, fitted_values, m, response)
  stand = stand %>%
    mutate(Estimate = c(unname(fit$beta), unname(fit$sigmasq)), .after = Variable)

  return(list("Estimates" = stand,
              "P-Value" = Hypothesis_Test(my_data, example, fitted_values, m, response, beta_i)))
}
