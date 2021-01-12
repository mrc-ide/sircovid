## Create and update rt data set
reference_data_rt <- function() {
  load_reference("data/rt.rds", {
    p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- carehomes$new(p, 0, np, seed = 1L)

    initial <- carehomes_initial(mod$info(), 10, p)
    mod$set_state(initial$state, initial$step)
    mod$set_index(integer(0))
    index_S <- mod$info()$index$S
    index_prob_strain <- mod$info()$index$prob_strain
    index <- c(index_S, index_prob_strain)
    
    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(initial$step, end, by = 1 / p$dt)

    set.seed(1)
    x <- dust::dust_iterate(mod, steps)
    
    y <- x[index_S, , ]
    prob_strain <- x[index_prob_strain, , ]
    
    rt_1 <- carehomes_Rt(steps, y[, 1, ], prob_strain[, 1, ], p)
    rt_all <- carehomes_Rt_trajectories(steps, y, prob_strain, p)

    list(inputs = list(steps = steps, y = y, prob_strain = prob_strain, p = p),
         outputs = list(rt_1 = rt_1, rt_all = rt_all))
  })
}


calculate_rt_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i <- grep("S_", rownames(dat$trajectories$state))
  S <- dat$trajectories$state[i, , ]
  carehomes_Rt_trajectories(dat$trajectories$step, S, p)
}
