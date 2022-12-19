#' Title
#'
#' @param my_data my_data the data set you're using
#' @param example the name of your example: will be one of either "culcita", "ctsib", "epilepsy", or "tortoise"
#' @param fitted_values fitted values of your linear predictor under the null hypothesis without random effects
#' @param m either the number of repeated observations from n individuals (repeated measures study)
#' or the number of subjects within n groups (grouped effects data)
#' @param response the response variable of your data set in ""
#' @param beta_i the parameter of interest for which we wish to test the null hypothesis for, default value of 0 for the intercept
#'
#' @return a list with a tibble of the estimates, standard errors, and confidence intervals and a p-value for the beta under the null hypothesis
#' @export
#'
#' @examples
#' # load data set
#' data("tortoise")
#'
#' #Fit model to gopher tortoises data
#' tortoise_fit <- run_model(data = tortoise, example = "tortoise")
#'
#' # Obtain the fitted values
#' x1 <- tortoise$prev
#' x2 <- tortoise$Area
#'
#' fitted_values <- tortoise_fit$beta[[1]] + tortoise_fit$beta[[2]] * x1 + log(x2) #fitted values without random effects
#'
#' # Find summary
#' fit_summary(tortoise, "tortoise", fitted_values, 3, "shells", "prev")

fit_summary = function(my_data, example, fitted_values, m, response, beta_i = 0){
  fit = run_model(my_data, example)
  stand = SE_CI(my_data, example, fitted_values, m, response)
  stand = stand %>%
    mutate(Estimate = c(unname(fit$beta), unname(fit$sigmasq)), .after = Variable)

  return(list("Estimates" = stand,
              "P-Value" = Hypothesis_Test(my_data, example, fitted_values, m, response, beta_i)))
}
