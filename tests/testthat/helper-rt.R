## Create and update rt data set
reference_data_rt <- function() {
  path <- "data/rt.rds"
  if (!file.exists(path)) {
    message(sprintf("Reference data '%s' does not exist - generating",
                    path))
    p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
    np <- 3L
    mod <- carehomes$new(p, 0, np, seed = 1L)

    initial <- carehomes_initial(mod$info(), 10, p)
    mod$set_state(initial$state, initial$step)
    mod$set_index(integer(0))
    index <- mod$info()$index$S

    end <- sircovid_date("2020-05-01") / p$dt
    steps <- seq(initial$step, end, by = 1 / p$dt)

    set.seed(1)
    y <- dust::dust_iterate(mod, steps, index)
    rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
    rt_all <- carehomes_Rt_trajectories(steps, y, p)

    data <- list(inputs = list(steps = steps, y = y, p = p),
                 outputs = list(rt_1 = rt_1, rt_all = rt_all))

    dir.create(dirname(path), FALSE, TRUE)
    saveRDS(data, path, version = 2L)
  }
  readRDS(path)
}


calculate_rt_simple <- function(dat) {
  p <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  i <- grep("S_", rownames(dat$trajectories$state))
  S <- dat$trajectories$state[i, , ]
  carehomes_Rt_trajectories(dat$trajectories$step, S, p)
}
