context("gpu")

test_that("can generate GPU interface and pass arguments", {
  skip_if_not_installed("mockery")
  mock_odin_dust <- mockery::mock(carehomes, cycle = TRUE)

  mockery::stub(compile_gpu, "odin.dust::odin_dust_", mock_odin_dust)

  expect_identical(
    compile_gpu(),
    carehomes)
  mockery::expect_called(mock_odin_dust, 1L)
  expect_equal(
    mockery::mock_args(mock_odin_dust)[[1]],
    list(sircovid_file("odin/carehomes.R"), real_t = "float", gpu = TRUE))

  compile_gpu(real_t = "double", rewrite_dims = TRUE,
              gpu = FALSE, gpu_generate = TRUE)
  expect_equal(
    mockery::mock_args(mock_odin_dust)[[2]],
    list(sircovid_file("odin/carehomes.R"), real_t = "double", gpu = FALSE,
         rewrite_dims = TRUE, gpu_generate = TRUE))
})


test_that("can run the gpu model on the cpu", {
  skip_unless_ci()
  ## TODO; see #286
  skip_on_os("linux")
  skip_if_not_installed("odin")
  skip_if_not_installed("odin.dust")

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = TRUE)
  expect_equal(gen$public_methods$name(), "carehomes")

  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            initial_I = 20)

  mod_cpu <- gen$new(p, 0, 5, seed = 1L)
  mod_gpu <- gen$new(p, 0, 5, seed = 1L, device_id = 0)

  end <- sircovid_date("2020-07-31") / p$dt
  info <- mod_cpu$info()

  expect_equal(mod_gpu$info(), mod_cpu$info())

  initial <- carehomes_initial(info, 10, p)
  index <- c(carehomes_index(info)$run,
             deaths_carehomes = info$index[["D_carehomes_tot"]],
             deaths_comm = info$index[["D_comm_tot"]],
             deaths_hosp = info$index[["D_hosp_tot"]],
             admitted = info$index[["cum_admit_conf"]],
             diagnoses = info$index[["cum_new_conf"]],
             sympt_cases = info$index[["cum_sympt_cases"]],
             sympt_cases_over25 = info$index[["cum_sympt_cases_over25"]])

  mod_cpu$set_state(initial$state, initial$step)
  mod_gpu$set_state(initial$state, initial$step)
  mod_cpu$set_index(index)
  mod_gpu$set_index(index)

  res_cpu <- mod_cpu$run(end)
  res_gpu <- mod_gpu$run(end, device = TRUE)

  expect_equal(res_gpu, res_cpu)
})
