template <typename real_t>
HOSTDEVICE
real_t ll_nbinom(real_t data, real_t model, real_t kappa, real_t exp_noise,
                 dust::rng_state_t<real_t>& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_t mu = model + dust::distr::rexp(rng_state, exp_noise);
  return dust::dnbinom_mu(data, kappa, mu, true);
}

// [[odin.dust::compare_data(icu = real_t, deaths = real_t)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_t compare(const typename T::real_t * state,
                           const typename T::data_t& data,
                           const typename T::internal_t internal,
                           std::shared_ptr<const typename T::shared_t> shared,
                           dust::rng_state_t<typename T::real_t>& rng_state) {
  typedef typename T::real_t real_t;
  real_t ll_icu = ll_nbinom(data.icu, odin(phi_ICU) * odin(I_ICU_tot),
                            odin(kappa_ICU), odin(exp_noise),
                            rng_state);
  real_t ll_deaths = ll_nbinom(data.deaths, odin(phi_death) * odin(D_inc),
                               odin(kappa_death), odin(exp_noise),
                               rng_state);
  return ll_icu + ll_deaths;
}
