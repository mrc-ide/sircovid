## Create and update rt data set
reference_data_lancelot_rt <- function() {
  load_reference("data/lancelot_rt.rds", {
    p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- lancelot$new(p, 0, np, seed = 1L)

    initial <- lancelot_initial(mod$info(), 10, p)
    mod$update_state(state = initial)
    mod$set_index(integer(0))
    index <- mod$info()$index$S

    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(0, end, by = 1 / p$dt)

    set.seed(1)
    mod$set_index(index)
    y <- mod$simulate(steps)
    rt_1 <- lancelot_Rt(steps, y[, 1, ], p)
    rt_all <- lancelot_Rt_trajectories(steps, y, p)

    list(inputs = list(steps = steps, y = y, p = p),
         outputs = list(rt_1 = rt_1, rt_all = rt_all))
  })
}


calculate_lancelot_rt_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i <- grep("S_", rownames(dat$trajectories$state))
  S <- dat$trajectories$state[i, , ]
  lancelot_Rt_trajectories(dat$trajectories$step, S, p)
}
