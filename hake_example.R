library(TMB)

library(tidyverse)

library(here)

hake <- read.csv(here("data","hake.csv"))  |> 
  janitor::clean_names() 

compile(here("src","hake.cpp"))               # Compile the C++ file

dyn.load(dynlib(here("src","hake")))

# on.exit(dyn.unload(here("src","hake")))


model <-  MakeADFun(
  data = list(
    catches = hake$catch,
    index = hake$abundance_index,
    years = nrow(hake),
    sigma_proc = 0.025,
    sigma_catch = 0.05
  ),
  parameters = list(
    log_r = log(0.2),
    log_k = log(20 * max(hake$catch)),
    log_index_q = log(0.001),
    log_init_dep = log(1),
    log_sigma_obs = log(0.5),
    effish = rep(-2, nrow(hake)),
    log_pros_error = rep(0, nrow(hake))
  ),
  random = "log_pros_error"
)

fit <- nlminb(model$par, model$fn, model$gr,control=list(iter.max=10000,eval.max=100000))

fit

sd_report <- sdreport(model)

estimates <- summary(sd_report) |> 
  as.data.frame() |> 
  rownames_to_column(var = "param") |> 
  janitor::clean_names()

index_hat <- estimates |> 
  filter(str_detect(param, "log_index_hat")) |> 
  mutate(year = hake$year,
         index_hat = exp(estimate),
         upper_index_hat = exp(estimate + 1.96 * std_error),
         lower_index_hat = exp(estimate - 1.96 * std_error)) |> 
  select(year, index_hat, upper_index_hat, lower_index_hat)
       


fitted_hake <- hake |> 
  left_join(index_hat, by = "year") 

ref_points <- estimates |> 
  filter(str_detect(param, "msy")) |> 
  mutate(param = str_remove_all(param,"\\..*$")) |> 
  select(param, estimate) |> 
  group_by(param) |> 
  mutate(time_step = 1:length(param)) |> 
  ungroup() |> 
  pivot_wider(names_from = param, values_from = estimate)


index_fit_plot <- fitted_hake |>
  ggplot() +
  geom_point(aes(year, abundance_index), size = 4) +
  geom_pointrange(
    aes(year, index_hat, ymax = upper_index_hat, ymin = lower_index_hat),
    color = "red",
    alpha = 0.5
  ) + 
  scale_y_continuous(limits = c(0, NA))


kobe_plot <- ref_points |> 
  ggplot(aes(b_bmsy, u_umsy, fill = time_step)) + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_path() + 
  geom_point(shape = 21, size = 4) + 
  scale_fill_viridis_c() + 
  scale_x_continuous(limits = c(0, NA), expand = expansion(c(0,.05))) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0,0.05)))
