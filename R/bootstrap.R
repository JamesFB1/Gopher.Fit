bootsrap = function(my_data = tortoise){
my_data = my_data %>%
  rowid_to_column("ID") %>%
  select(-type, -density)

fit = run_model(my_data, example = "tortoise")

my_data = my_data %>%
  mutate(.fitted = fit$beta[[2]] * prev)

B = 1000

bootstrap = tibble(B = 1:B) %>%
  crossing(my_data) %>%
  mutate(z = rep(rnorm(n()/3, mean = 0, sd = sqrt(fit$sigmasq)), each = 3),
         eta = fit$beta[[1]] + log(Area) + .fitted + z,
         shells = rpois(n(), lambda = exp(eta))) %>%
  nest_by(B) %>%
  summarize(Intercept = run_model(data = data,example = "tortoise")$beta[[1]],
            Prev = run_model(data = data,example = "tortoise")$beta[[2]],
            Sigmasq = run_model(data = data,example = "tortoise")$sigmasq) %>%
  ungroup()

Intercept = bootstrap %>%
            summarize(SE = sd(Intercept),
            Lower95 = quantile(Intercept, 0.025),
            Upper95 = quantile(Intercept, 0.975))

Prev = bootstrap %>%
            summarize(SE = sd(Prev),
            Lower95 = quantile(Prev, 0.025),
            Upper95 = quantile(Prev, 0.975))

Sigmasq = bootstrap %>%
            summarize(SE = sd(Sigmasq),
            Lower95 = quantile(Sigmasq, 0.025),
            Upper95 = quantile(Sigmasq, 0.975))
return(list(Intercept = Intercept,
            Prev = Prev,
            Sigmasq = Sigmasq))
}
