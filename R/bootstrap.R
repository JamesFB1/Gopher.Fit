tortoise = tortoise %>%
  rowid_to_column("ID") %>%
  select(-type, -density)

fit = run_model(tortoise, example = "tortoise")

tortoise = tortoise %>%
  mutate(.fitted = fit$beta[[2]] * prev)

B = 10000
n = nrow(tortoise)
bootstrap = tibble(B = 1:B) %>%
  crossing(tortoise) %>%
  mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3),
         eta = fit$beta[[1]] + log(Area) + .fitted + z,
         shells = rpois(n(), lambda = exp(eta))) %>%
  nest_by(B)


