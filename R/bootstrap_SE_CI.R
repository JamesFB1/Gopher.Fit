SE_CI = function(my_data = tortoise){

  my_data = my_data %>%
    rowid_to_column("ID") %>%  # add ID column
    select(-type, -density)    # delete type and density columns

  fit = run_model(my_data, example = "tortoise")  # fit the model

  my_data = my_data %>%
    mutate(.fitted = fit$beta[[2]] * prev)  # make a column for beta_1 * prev

  ## Generate bootstrap samples
  B = 1000

  bootstrap = tibble(B = 1:B) %>%
    crossing(my_data) %>%
    mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3),
           eta = fit$beta[[1]] + log(Area) + .fitted + z,  # eta_ij* = beta_0 + beta_1xij1 + log(xi2) + zi*
           shells = rpois(n(), lambda = exp(eta))) %>%  # Yij*~Poisson(lambda) where lambda = mu_ij* = exp(eta_ij*)
    nest_by(B) %>%
    summarize(tidy(suppressMessages(glmer(shells ~ prev + offset(log(Area)) + factor(year) + (1|Site),
                                          family = poisson, data = data))), .groups = "drop")  # fit the model the way Prof. Bonner did

    # summarize(Intercept = run_model(data = data,example = "tortoise")$beta[[1]],  # refit the model, bootstrapped estimate of beta_0
    #           Prev = run_model(data = data,example = "tortoise")$beta[[2]],       # bootstrapped estimate of beta_1
    #           Sigmasq = run_model(data = data,example = "tortoise")$sigmasq) %>%  # bootstrapped estimate of sigma_sqr
    # ungroup()

  ## Compute bootstrap standard errors and confidence intervals for beta_0, beta_1 and sigma_sqr
  Intercept = bootstrap %>%
    filter(term == "(Intercept)") %>%  # filter for the Intercept term (beta_0)
    summarize(SE = sd(estimate),  # estimate column has parameter estimate for Intercept term
              Lower95 = quantile(estimate, 0.025),
              Upper95 = quantile(estimate, 0.975))

  Prev = bootstrap %>%
    filter(term == "prev") %>%  # do the same for beta_1
    summarize(SE = sd(estimate),
              Lower95 = quantile(estimate, 0.025),
              Upper95 = quantile(estimate, 0.975))

  Sigmasq = bootstrap %>%
    filter(term == "sd__(Intercept)") %>%  # not quite sure about this, but checked estimation_functions.R line 103 to try and extract sigmasq from output
    mutate(estimate = ifelse(estimate == 0, 0.100, estimate)) %>%  # followed what was done in estimation_functions.R - if sigmasq == 0, set it to 0.1
    summarize(SE = sd(estimate),
              Lower95 = quantile(estimate, 0.025),
              Upper95 = quantile(estimate, 0.975))


  # Intercept <- bootstrap %>%
  #             summarize(SE = sd(Intercept),
  #             Lower95 = quantile(Intercept, 0.025),
  #             Upper95 = quantile(Intercept, 0.975))
  #
  # Prev <- bootstrap %>%
  #             summarize(SE = sd(Prev),
  #             Lower95 = quantile(Prev, 0.025),
  #             Upper95 = quantile(Prev, 0.975))
  #
  # Sigmasq <- bootstrap %>%
  #             summarize(SE = sd(Sigmasq),
  #             Lower95 = quantile(Sigmasq, 0.025),
  #             Upper95 = quantile(Sigmasq, 0.975))

  return(list(Intercept = Intercept,
              Prev = Prev,
              Sigmasq = Sigmasq))
}
