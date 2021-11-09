template <typename real_type, typename rng_state_type>
__host__ __device__
real_type ll_nbinom(real_type data, real_type model, real_type kappa,
                    real_type exp_noise,
                    rng_state_type& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_type mu = model +
    dust::random::exponential<real_type>(rng_state, exp_noise);
  return dust::density::negative_binomial_mu(data, kappa, mu, true);
}

// [[odin.dust::compare_data(icu = real_type, deaths = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
compare(const typename T::real_type * state,
        const typename T::data_type& data,
        const typename T::internal_type internal,
        std::shared_ptr<const typename T::shared_type> shared,
        typename T::rng_state_type& rng_state) {
  typedef typename T::real_type real_type;
  real_type ll_icu = ll_nbinom(data.icu, odin(phi_ICU) * odin(I_ICU_tot),
                            odin(kappa_ICU), odin(exp_noise),
                            rng_state);
  real_type ll_deaths = ll_nbinom(data.deaths, odin(phi_death) * odin(D_inc),
                               odin(kappa_death), odin(exp_noise),
                               rng_state);
  return ll_icu + ll_deaths;
}
