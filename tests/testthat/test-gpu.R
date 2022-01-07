context("gpu")

test_that("can generate GPU interface and pass arguments", {
  skip_if_not_installed("mockery")
  mock_odin_dust <- mockery::mock(lancelot, cycle = TRUE)

  mockery::stub(compile_gpu, "odin.dust::odin_dust_", mock_odin_dust)

  expect_identical(
    compile_gpu(),
    lancelot)
  mockery::expect_called(mock_odin_dust, 1L)
  expect_equal(
    mockery::mock_args(mock_odin_dust)[[1]],
    list(sircovid_file("odin/lancelot.R"), real_type = "float", gpu = TRUE))

  compile_gpu(real_type = "double", rewrite_dims = TRUE,
              gpu = FALSE, gpu_generate = TRUE)
  expect_equal(
    mockery::mock_args(mock_odin_dust)[[2]],
    list(sircovid_file("odin/lancelot.R"), real_type = "double", gpu = FALSE,
         rewrite_dims = TRUE, gpu_generate = TRUE))
})


test_that("can run the gpu model on the cpu", {
  skip_unless_ci()
  skip_if_not_installed("odin")
  skip_if_not_installed("odin.dust")

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = FALSE)
  expect_equal(gen$public_methods$name(), "lancelot")

  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 20)

  mod_cpu <- gen$new(p, 0, 5, seed = 1L)
  mod_gpu <- gen$new(p, 0, 5, seed = 1L, gpu_config = 0)

  end <- sircovid_date("2020-07-31") / p$dt
  info <- mod_cpu$info()

  expect_equal(mod_gpu$info(), mod_cpu$info())

  initial <- lancelot_initial(info, 10, p)
  index <- c(lancelot_index(info)$run,
             deaths_carehomes = info$index[["D_carehomes_tot"]],
             deaths_comm = info$index[["D_comm_tot"]],
             deaths_hosp = info$index[["D_hosp_tot"]],
             admitted = info$index[["cum_admit_conf"]],
             diagnoses = info$index[["cum_new_conf"]],
             sympt_cases = info$index[["cum_sympt_cases"]],
             sympt_cases_over25 = info$index[["cum_sympt_cases_over25"]])

  mod_cpu$update_state(state = initial)
  mod_gpu$update_state(state = initial)
  mod_cpu$set_index(index)
  mod_gpu$set_index(index)

  res_cpu <- mod_cpu$run(end)
  res_gpu <- mod_gpu$run(end)

  expect_equal(res_gpu, res_cpu)
})


test_that("Can run the gpu compare on the cpu", {
  skip_unless_ci()
  skip_if_not_installed("odin")
  skip_if_not_installed("odin.dust")

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = FALSE)

  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england",
                              initial_seed_size = 20)
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               start_date, pars$dt)

  np <- 10
  mod_cpu <- gen$new(pars, 0, np, seed = 1L)
  mod_gpu <- gen$new(pars, 0, np, seed = 1L, gpu_config = 0)

  initial <- lancelot_initial(mod_cpu$info(), np, pars)
  mod_cpu$update_state(state = initial)
  mod_gpu$update_state(state = initial)

  mod_cpu$set_data(dust::dust_data(data, "step_end"))
  mod_gpu$set_data(dust::dust_data(data, "step_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y_cpu <- mod_cpu$run(data$step_end[[i]])
  y_gpu <- mod_gpu$run(data$step_end[[i]])
  expect_equal(mod_cpu$compare_data(), mod_gpu$compare_data())
})
