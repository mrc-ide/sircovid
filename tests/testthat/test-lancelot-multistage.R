test_that("can run a multistage particle filter", {
  ## This is not really proper test because we need to make sure that
  ## things behave as expected.  But that's for you to write, I think.

  ## First set of parameters is the phase 1; we have one strain
  p1 <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  ## Second set is the multistrain set; we're also seeding in from
  ## shortly after the change.
  p2 <- lancelot_parameters(sircovid_date("2020-01-01"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date = sircovid_date("2020-03-15"),
                            strain_seed_size = 100,
                            strain_seed_pattern = rep(1, 4),
                            cross_immunity = 0,
                            n_real_strains = 2)

  ## Then we put them together as a single multistage parameter set.
  ## The data here is important; that's the changeover date.  So betas
  ## set in p1 before this, or in p2 after it will be discarded.  We'd
  ## want to make sure appropriate things are shared here when
  ## generating the parameters of course.
  p <- mcstate::multistage_parameters(
    p1,
    epochs = list(
      mcstate::multistage_epoch(sircovid_date("2020-03-01"),
                                p2, inflate_state_strains)))

  ## The usual test data set for the filter:
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               0, p1$dt)

  ## Construct the filter with the existing functions
  filter <- helper_lancelot_particle_filter(data, n_particles = 10, seed = 1L)

  ## The first case runs the filter just on the first parameter set
  ## over the data
  expect_silent(filter$run(p1))

  ## This case runs the multistage filter!
  filter <- helper_lancelot_particle_filter(data, n_particles = 10, seed = 1L)
  expect_silent(filter$run(p))
})
