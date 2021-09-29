## Create and update rt data set
reference_data_carehomes_ifr_t <- function() {
  load_reference("data/ifr_t.rds", {
    p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- carehomes$new(p, 0, np, seed = 1L)

    initial <- carehomes_initial(mod$info(), 10, p)
    mod$set_state(initial$state, initial$step)
    mod$set_index(integer(0))
    index_S <- mod$info()$index$S
    index_I_weighted <- mod$info()$index$I_weighted
    index <- c(index_S, index_I_weighted)

    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(initial$step, end, by = 1 / p$dt)

    set.seed(1)
    mod$set_index(index)
    y <- mod$simulate(steps)
    S <- y[seq_len(length(index_S)), , ]
    I_weighted <- y[-seq_len(length(index_S)), , ]

    ifr_t_1 <- carehomes_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
    ifr_t_all <- carehomes_ifr_t_trajectories(steps, S, I_weighted, p)

    list(inputs = list(steps = steps, S = S, I_weighted = I_weighted, p = p),
         outputs = list(ifr_t_1 = ifr_t_1, ifr_t_all = ifr_t_all))
  })
}


calculate_carehomes_ifr_t_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i_S <- grep("S_", rownames(dat$trajectories$state))
  S <- dat$trajectories$state[i_S, , ]
  i_I_weighted <- grep("I_weighted_", rownames(dat$trajectories$state))
  I_weighted <- dat$trajectories$state[i_I_weighted, , ]
  carehomes_ifr_t_trajectories(dat$trajectories$step, S, I_weighted, p)
}
