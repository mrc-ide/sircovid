## Create and update rt data set
reference_data_carehomes_rt <- function() {
  load_reference("data/carehomes_rt.rds", {
    p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- carehomes$new(p, 0, np, seed = 1L)

    initial <- carehomes_initial(mod$info(), 10, p)
    mod$update_state(state = initial$state, step = initial$step)
    mod$set_index(integer(0))
    index <- mod$info()$index$S

    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(initial$step, end, by = 1 / p$dt)

    set.seed(1)
    mod$set_index(index)
    y <- mod$simulate(steps)
    rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
    rt_all <- carehomes_Rt_trajectories(steps, y, p)

    list(inputs = list(steps = steps, y = y, p = p),
         outputs = list(rt_1 = rt_1, rt_all = rt_all))
  })
}


calculate_carehomes_rt_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i <- grep("S_", rownames(dat$trajectories$state))
  S <- dat$trajectories$state[i, , ]
  carehomes_Rt_trajectories(dat$trajectories$step, S, p)
}
