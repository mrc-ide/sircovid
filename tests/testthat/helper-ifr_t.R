## Create and update rt data set
reference_data_ifr_t <- function() {
  load_reference("data/ifr_t.rds", {
    p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- carehomes$new(p, 0, np, seed = 1L)

    initial <- carehomes_initial(mod$info(), 10, p)
    mod$set_state(initial$state, initial$step)
    mod$set_index(integer(0))
    index <- mod$info()$index$infections_inc

    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(initial$step, end, by = 1 / p$dt)

    set.seed(1)
    y <- dust::dust_iterate(mod, steps, index)
    ifr_t_1 <- carehomes_ifr_t(steps, y[, 1, ], p)
    ifr_t_all <- carehomes_ifr_t_trajectories(steps, y, p)

    list(inputs = list(steps = steps, y = y, p = p),
         outputs = list(ifr_t_1 = ifr_t_1, ifr_t_all = ifr_t_all))
  })
}


calculate_ifr_t_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i <- grep("infections_inc_", rownames(dat$trajectories$state))
  infections_inc <- dat$trajectories$state[i, , ]
  carehomes_ifr_t_trajectories(dat$trajectories$step, infections_inc, p)
}
