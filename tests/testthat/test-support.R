context("date")

test_that("sircovid_date can convert to dates into 2020", {
  expect_equal(sircovid_date("2020-01-01"), 1)
  expect_equal(sircovid_date("2020-12-31"), 366) # leap year!
  r <- seq(as_date("2020-01-01"), as_date("2020-12-31"), by = 1)
  expect_equal(sircovid_date(r), 1:366)

  expect_equal(sircovid_date_as_date(sircovid_date(r)), r)

  expect_error(sircovid_date("2019-01-01"),
               "Negative dates, sircovid_date likely applied twice")
  expect_error(sircovid_date(c("2020-01-01", "2019-01-01")),
               "Negative dates, sircovid_date likely applied twice")
})


test_that("assert sircovid date throws on non sircovid dates", {
  expect_silent(assert_sircovid_date(1))
  expect_error(assert_sircovid_date(as_date("2020-01-01")),
               "'date' must be numeric - did you forget sircovid_date()?",
               fixed = TRUE)
})


test_that("helper function can convert somewhat helpfully", {
  expect_equal(as_sircovid_date("2020-02-01"), 32)
  expect_equal(as_sircovid_date(as_date("2020-02-01")), 32)
  expect_equal(as_sircovid_date(32), 32)
})


test_that("read population", {
  cache$population <- NULL
  n <- sircovid_population("south_west")
  expect_s3_class(cache$population, "data.frame")
  expected <- c(296855, 325583, 307586, 304574, 339548, 337183, 327525,
                327910, 311774, 374176, 403024, 380109, 338931, 332910,
                333351, 225646, 339312)
  expect_equal(n, expected)
})


test_that("reject invalid regions", {
  expect_error(sircovid_population("oxfordshire"),
               "Population not found for 'oxfordshire': must be one of 'uk',")
  expect_error(sircovid_population(NULL),
               "'region' must not be NULL")
})


test_that("Downcase region name", {
  expect_identical(sircovid_population("UK"), sircovid_population("uk"))
})


test_that("Use constant age bins", {
  expect_equal(
    sircovid_age_bins(),
    list(start = seq(0, 80, by = 5),
         end = c(seq(4, 79, by = 5), 100)))
})
