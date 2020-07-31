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
