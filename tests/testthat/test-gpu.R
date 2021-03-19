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
