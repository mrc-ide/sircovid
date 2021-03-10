#include <dust/dust.hpp>
#include <dust/interface.hpp>

template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum1(const container x, size_t from, size_t to);
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum3(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int dim_x_1, int dim_x_12);
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum4(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int from_l, int to_l, int dim_x_1, int dim_x_12, int dim_x_123);
template <typename real_t, typename T, typename U>
HOSTDEVICE real_t fmodr(T x, U y) {
  real_t tmp = std::fmod(static_cast<real_t>(x), static_cast<real_t>(y));
  if (tmp * y < 0) {
    tmp += y;
  }
  return tmp;
}

// These exist to support the model on the gpu, as in C++14 std::min
// and std::max are constexpr and error without --expt-relaxed-constexpr
template <typename T>
HOSTDEVICE T odin_min(T x, T y) {
  return x < y ? x : y;
}

template <typename T>
HOSTDEVICE T odin_max(T x, T y) {
  return x > y ? x : y;
}
template <typename real_t>
real_t ll_nbinom(real_t data, real_t model, real_t kappa, real_t exp_noise,
                 dust::rng_state_t<real_t>& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_t mu = model + dust::distr::rexp(rng_state, exp_noise);
  return dust::dnbinom(data, kappa, mu, true);
}

template <typename real_t>
real_t ll_binom(real_t data_x, real_t data_size, real_t model_prob) {
  if (std::isnan(data_x) || std::isnan(data_size)) {
    return 0;
  }
  return dust::dbinom(data_x, data_size, model_prob, true);
}

template <typename real_t>
real_t ll_betabinom(real_t data_x, real_t data_size, real_t model_prob,
                    real_t rho) {
  if (std::isnan(data_x) || std::isnan(data_size)) {
    return 0;
  }
  return dust::dbetabinom(data_x, data_size, model_prob, rho, true);
}

template <typename real_t>
real_t test_prob_pos(real_t pos, real_t neg, real_t sensitivity,
                     real_t specificity, real_t exp_noise,
                     dust::rng_state_t<real_t>& rng_state) {
  // We add some exponential noise to the number of positives and negatives
  // to help ensure prob_pos is not 0 or 1. If e.g. prob_pos were 0 and there
  // were individuals who tested positive, this would result in a weight of 0
  // for a particle. If all particles have weights of 0, the particle filter
  // breaks. The exponential noise produces small non-zero weights in these
  // circumstances to prevent the particle filter from breaking.
  pos += dust::distr::rexp(rng_state, exp_noise);
  neg += dust::distr::rexp(rng_state, exp_noise);
  return (sensitivity * pos + (1 - specificity) * neg) / (pos + neg);
}

// Lots of data items here. Even though some of these are counts (in
// fact, most of them are counts) we'll pull them in as floating point
// numbers (double or float) because that will help with the NA
// handling. Once this is working reliably we could cast these
// explicitly as ints, but there is not likely to be any major gain
// here.
//
// [[odin.dust::compare_data(icu = real_t)]]
// [[odin.dust::compare_data(general = real_t)]]
// [[odin.dust::compare_data(hosp = real_t)]]
// [[odin.dust::compare_data(deaths_hosp = real_t)]]
// [[odin.dust::compare_data(deaths_comm = real_t)]]
// [[odin.dust::compare_data(deaths_carehomes = real_t)]]
// [[odin.dust::compare_data(deaths = real_t)]]
// [[odin.dust::compare_data(deaths_non_hosp = real_t)]]
// [[odin.dust::compare_data(admitted = real_t)]]
// [[odin.dust::compare_data(diagnoses = real_t)]]
// [[odin.dust::compare_data(all_admission = real_t)]]
// [[odin.dust::compare_data(npos_15_64 = real_t)]]
// [[odin.dust::compare_data(ntot_15_64 = real_t)]]
// [[odin.dust::compare_data(pillar2_pos = real_t)]]
// [[odin.dust::compare_data(pillar2_tot = real_t)]]
// [[odin.dust::compare_data(pillar2_cases = real_t)]]
// [[odin.dust::compare_data(pillar2_over25_pos = real_t)]]
// [[odin.dust::compare_data(pillar2_over25_tot = real_t)]]
// [[odin.dust::compare_data(pillar2_over25_cases = real_t)]]
// [[odin.dust::compare_data(react_pos = real_t)]]
// [[odin.dust::compare_data(react_tot = real_t)]]
// [[odin.dust::compare_data(strain_non_variant = real_t)]]
// [[odin.dust::compare_data(strain_tot = real_t)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_t compare(const typename T::real_t * state,
                           const typename T::data_t& data,
                           const typename T::internal_t internal,
                           std::shared_ptr<const typename T::shared_t> shared,
                           dust::rng_state_t<typename T::real_t>& rng_state) {
  typedef typename T::real_t real_t;

  // State variables; these largely correspond to the quantities in data
  const real_t model_icu = state[9];
  const real_t model_general = state[10];
  const real_t model_hosp = model_icu + model_general;
  const real_t model_deaths_carehomes = state[16];
  const real_t model_deaths_comm = state[14];
  const real_t model_deaths_hosp = state[17];
  const real_t model_admitted = state[1];
  const real_t model_diagnoses = state[2];
  const real_t model_all_admission = model_admitted + model_diagnoses;
  const real_t model_sero_pos = state[19];
  const real_t model_sympt_cases = state[23];
  const real_t model_sympt_cases_over25 = state[24];
  const real_t model_sympt_cases_non_variant_over25 =
    state[25];
  const real_t model_react_pos = state[26];

  // This is used over and over
  const real_t exp_noise = shared->exp_noise;

  const real_t pillar2_negs =
    shared->p_NC * (shared->N_tot_all - model_sympt_cases);
  const real_t model_pillar2_prob_pos =
    test_prob_pos(model_sympt_cases,
                  pillar2_negs,
                  shared->pillar2_sensitivity,
                  shared->pillar2_specificity,
                  exp_noise,
                  rng_state);

  const real_t pillar2_over25_negs =
    shared->p_NC * (shared->N_tot_over25 - model_sympt_cases_over25);
  const real_t model_pillar2_over25_prob_pos =
    test_prob_pos(model_sympt_cases_over25,
                  pillar2_over25_negs,
                  shared->pillar2_sensitivity,
                  shared->pillar2_specificity,
                  exp_noise,
                  rng_state);

  const real_t N_tot_react = shared->N_tot_react; // sum(pars$N_tot[2:18])
  const real_t model_react_prob_pos =
    test_prob_pos(model_react_pos,
                  N_tot_react - model_react_pos,
                  shared->react_sensitivity,
                  shared->react_specificity,
                  exp_noise,
                  rng_state);

  // serology
  const real_t model_sero_prob_pos =
    test_prob_pos(model_sero_pos,
                  shared->N_tot_15_64 - model_sero_pos,
                  shared->sero_sensitivity,
                  shared->sero_specificity,
                  exp_noise,
                  rng_state);

  // Strain
  real_t strain_sensitivity = 1.0;
  real_t strain_specificity = 1.0;
  const real_t model_strain_over25_prob_pos =
    test_prob_pos(model_sympt_cases_non_variant_over25,
                  model_sympt_cases_over25 -
                  model_sympt_cases_non_variant_over25,
                  strain_sensitivity,
                  strain_specificity,
                  exp_noise,
                  rng_state);

  // Note that in ll_nbinom, the purpose of exp_noise is to allow a
  // non-zero probability when the model value is 0 and the observed
  // value is non-zero (i.e. there is overreporting)
  const real_t ll_icu =
    ll_nbinom(data.icu, shared->phi_ICU * model_icu,
              shared->kappa_ICU, exp_noise, rng_state);
  const real_t ll_general =
    ll_nbinom(data.general, shared->phi_general * model_general,
              shared->kappa_general, exp_noise, rng_state);
  const real_t ll_hosp =
    ll_nbinom(data.hosp, shared->phi_hosp * model_hosp,
              shared->kappa_hosp, exp_noise, rng_state);

  // We will compute one of the following:
  // 1. ll_deaths_hosp, ll_deaths_carehomes and ll_deaths_comm
  // 2. ll_deaths_hosp and ll_deaths_non_hosp
  // 3. ll_deaths
  const real_t ll_deaths_hosp =
    ll_nbinom(data.deaths_hosp, shared->phi_death_hosp * model_deaths_hosp,
              shared->kappa_death_hosp, exp_noise, rng_state);
  const real_t ll_deaths_carehomes =
    ll_nbinom(data.deaths_carehomes,
              shared->phi_death_carehomes * model_deaths_carehomes,
              shared->kappa_death_carehomes, exp_noise, rng_state);
  const real_t ll_deaths_comm =
    ll_nbinom(data.deaths_comm, shared->phi_death_comm * model_deaths_comm,
              shared->kappa_death_comm, exp_noise, rng_state);
  const real_t ll_deaths_non_hosp =
    ll_nbinom(data.deaths_non_hosp,
              shared->phi_death_carehomes * model_deaths_carehomes +
                shared->phi_death_comm * model_deaths_comm,
              shared->kappa_death_non_hosp, exp_noise, rng_state);
  const real_t ll_deaths =
    ll_nbinom(data.deaths,
              shared->phi_death_hosp * model_deaths_hosp +
                shared->phi_death_carehomes * model_deaths_carehomes +
                shared->phi_death_comm * model_deaths_comm,
              shared->kappa_death, exp_noise, rng_state);

  const real_t ll_admitted =
    ll_nbinom(data.admitted, shared->phi_admitted * model_admitted,
              shared->kappa_admitted, exp_noise, rng_state);
  const real_t ll_diagnoses =
    ll_nbinom(data.diagnoses, shared->phi_diagnoses * model_diagnoses,
              shared->kappa_diagnoses, exp_noise, rng_state);
  const real_t ll_all_admission =
    ll_nbinom(data.all_admission, shared->phi_all_admission * model_all_admission,
              shared->kappa_all_admission, exp_noise, rng_state);

  const real_t ll_serology =
    ll_binom(data.npos_15_64, data.ntot_15_64, model_sero_prob_pos);

  const real_t ll_pillar2_tests =
    ll_betabinom(data.pillar2_pos, data.pillar2_tot,
                 model_pillar2_prob_pos, shared->rho_pillar2_tests);
  const real_t ll_pillar2_cases =
    ll_nbinom(data.pillar2_cases,
              shared->phi_pillar2_cases * model_sympt_cases,
              shared->kappa_pillar2_cases, exp_noise, rng_state);

  const real_t ll_pillar2_over25_tests =
    ll_betabinom(data.pillar2_over25_pos, data.pillar2_over25_tot,
                 model_pillar2_over25_prob_pos, shared->rho_pillar2_tests);
  const real_t ll_pillar2_over25_cases =
    ll_nbinom(data.pillar2_over25_cases,
              shared->phi_pillar2_cases * model_sympt_cases_over25,
              shared->kappa_pillar2_cases, exp_noise, rng_state);

  const real_t ll_react =
    ll_binom(data.react_pos, data.react_tot,
             model_react_prob_pos);
  const real_t ll_strain_over25 =
    ll_binom(data.strain_non_variant, data.strain_tot,
             model_strain_over25_prob_pos);

  return ll_icu + ll_general + ll_hosp + ll_deaths_hosp + ll_deaths_carehomes +
    ll_deaths_comm + ll_deaths_non_hosp + ll_deaths + ll_admitted +
    ll_diagnoses + ll_all_admission + ll_serology +
    ll_pillar2_tests + ll_pillar2_cases + ll_pillar2_over25_tests +
    ll_pillar2_over25_cases + ll_react + ll_strain_over25;
}
// [[dust::class(carehomes)]]
// [[dust::param(G_D_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(ICU_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_A_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_C_1_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_C_2_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_P_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_tot_15_64, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_tot_all, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_tot_over25, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_tot_react, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(exp_noise, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(hosp_transmission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(index_dose, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_A, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_C_1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_C_2, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_E, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_G_D, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_H_D, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_H_R, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_ICU_D, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_ICU_W_D, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_ICU_W_R, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_ICU_pre, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_P, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_PCR_pos, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_PCR_pre, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_W_D, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_W_R, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(k_sero_pos, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_ICU, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_admitted, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_all_admission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_death, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_death_carehomes, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_death_comm, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_death_hosp, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_death_non_hosp, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_diagnoses, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_general, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_hosp, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kappa_pillar2_cases, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(m, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_C, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_G_D_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_H_D_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_H_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_ICU_D_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_ICU_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_NC, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_W_D_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_sero_pos, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_star_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_ICU, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_admitted, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_all_admission, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_death_carehomes, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_death_comm, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_death_hosp, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_diagnoses, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_general, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_hosp, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi_pillar2_cases, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(pillar2_sensitivity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(pillar2_specificity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_G_D, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_H, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_H_D, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_ICU, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_ICU_D, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_W_D, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(psi_star, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(react_sensitivity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(react_specificity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rel_infectivity, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rel_p_hosp_if_sympt, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rel_p_sympt, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rel_susceptibility, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rho_pillar2_tests, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(sero_sensitivity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(sero_specificity, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(steps_per_day, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(strain_seed_step, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(strain_transmission, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vaccine_dose_step, has_default = FALSE, default_value = NULL, rank = 3, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vaccine_progression_rate_base, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(waning_rate, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_A, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_C_1, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_C_2, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_E, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_G_D, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_H_D, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_H_R, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_ICU_D, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_ICU_W_D, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_ICU_W_R, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_ICU_pre, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_P, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_PCR_pos, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_PCR_pre, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_U, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_W_D, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_W_R, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_sero_pos, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_sero_pre_1, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_sero_pre_2, has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(model_pcr_and_serology_user, has_default = TRUE, default_value = 1L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_sero_pre_1, has_default = TRUE, default_value = 0.5, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class carehomes {
public:
  typedef float real_t;
  struct data_t {
    real_t icu;
    real_t general;
    real_t hosp;
    real_t deaths_hosp;
    real_t deaths_comm;
    real_t deaths_carehomes;
    real_t deaths;
    real_t deaths_non_hosp;
    real_t admitted;
    real_t diagnoses;
    real_t all_admission;
    real_t npos_15_64;
    real_t ntot_15_64;
    real_t pillar2_pos;
    real_t pillar2_tot;
    real_t pillar2_cases;
    real_t pillar2_over25_pos;
    real_t pillar2_over25_tot;
    real_t pillar2_over25_cases;
    real_t react_pos;
    real_t react_tot;
    real_t strain_non_variant;
    real_t strain_tot;
  };
  struct shared_t {
    real_t G_D_transmission;
    real_t ICU_transmission;
    real_t I_A_transmission;
    real_t I_C_1_transmission;
    real_t I_C_2_transmission;
    real_t I_P_transmission;
    real_t N_tot_15_64;
    real_t N_tot_all;
    real_t N_tot_over25;
    real_t N_tot_react;
    std::vector<real_t> beta_step;
    int dim_E;
    int dim_E_12;
    int dim_E_123;
    int dim_G_D;
    int dim_G_D_123;
    int dim_H_D_unconf;
    int dim_H_D_unconf_123;
    int dim_H_R_unconf;
    int dim_H_R_unconf_123;
    int dim_ICU_D_unconf;
    int dim_ICU_D_unconf_123;
    int dim_ICU_W_D_unconf;
    int dim_ICU_W_D_unconf_123;
    int dim_ICU_W_R_unconf;
    int dim_ICU_W_R_unconf_123;
    int dim_ICU_pre_unconf;
    int dim_ICU_pre_unconf_123;
    int dim_I_A;
    int dim_I_A_123;
    int dim_I_C_1;
    int dim_I_C_1_123;
    int dim_I_C_2;
    int dim_I_C_2_123;
    int dim_I_P;
    int dim_I_P_123;
    int dim_R;
    int dim_T_PCR_pos;
    int dim_T_PCR_pos_123;
    int dim_T_PCR_pre;
    int dim_T_PCR_pre_123;
    int dim_T_sero_pos;
    int dim_T_sero_pos_123;
    int dim_T_sero_pre;
    int dim_T_sero_pre_123;
    int dim_W_D_unconf;
    int dim_W_D_unconf_123;
    int dim_W_R_unconf;
    int dim_W_R_unconf_123;
    int dim_beta_step;
    int dim_cum_n_S_vaccinated;
    int dim_p_G_D_step;
    int dim_p_H_D_step;
    int dim_p_H_step;
    int dim_p_ICU_D_step;
    int dim_p_ICU_step;
    int dim_p_W_D_step;
    int dim_p_star_step;
    int dim_rel_susceptibility;
    int dim_rel_susceptibility_1;
    int dim_rel_susceptibility_2;
    int dim_s_ij;
    int dim_strain_seed_step;
    int dim_strain_transmission;
    int dim_vaccine_dose_step;
    int dim_vaccine_dose_step_1;
    int dim_vaccine_dose_step_12;
    int dim_vaccine_dose_step_2;
    int dim_vaccine_dose_step_3;
    real_t dt;
    real_t exp_noise;
    real_t gamma_A;
    real_t gamma_C_1;
    real_t gamma_C_2;
    real_t gamma_E;
    real_t gamma_G_D;
    real_t gamma_H_D;
    real_t gamma_H_R;
    real_t gamma_ICU_D;
    real_t gamma_ICU_W_D;
    real_t gamma_ICU_W_R;
    real_t gamma_ICU_pre;
    real_t gamma_P;
    real_t gamma_PCR_pos;
    real_t gamma_PCR_pre;
    real_t gamma_U;
    real_t gamma_W_D;
    real_t gamma_W_R;
    real_t gamma_sero_pos;
    std::vector<real_t> gamma_sero_pre;
    real_t gamma_sero_pre_1;
    real_t gamma_sero_pre_2;
    real_t hosp_transmission;
    std::vector<int> index_dose;
    real_t initial_D_carehomes_inc;
    real_t initial_D_carehomes_tot;
    real_t initial_D_comm_inc;
    real_t initial_D_comm_tot;
    std::vector<real_t> initial_D_hosp;
    real_t initial_D_hosp_inc;
    real_t initial_D_hosp_tot;
    std::vector<real_t> initial_D_non_hosp;
    real_t initial_D_tot;
    std::vector<real_t> initial_E;
    std::vector<real_t> initial_G_D;
    std::vector<real_t> initial_H_D_conf;
    std::vector<real_t> initial_H_D_unconf;
    std::vector<real_t> initial_H_R_conf;
    std::vector<real_t> initial_H_R_unconf;
    std::vector<real_t> initial_ICU_D_conf;
    std::vector<real_t> initial_ICU_D_unconf;
    std::vector<real_t> initial_ICU_W_D_conf;
    std::vector<real_t> initial_ICU_W_D_unconf;
    std::vector<real_t> initial_ICU_W_R_conf;
    std::vector<real_t> initial_ICU_W_R_unconf;
    std::vector<real_t> initial_ICU_pre_conf;
    std::vector<real_t> initial_ICU_pre_unconf;
    real_t initial_ICU_tot;
    std::vector<real_t> initial_I_A;
    std::vector<real_t> initial_I_C_1;
    std::vector<real_t> initial_I_C_2;
    std::vector<real_t> initial_I_P;
    std::vector<real_t> initial_I_weighted;
    std::vector<real_t> initial_N_tot;
    real_t initial_N_tot2;
    real_t initial_N_tot3;
    std::vector<real_t> initial_R;
    std::vector<real_t> initial_S;
    std::vector<real_t> initial_T_PCR_neg;
    std::vector<real_t> initial_T_PCR_pos;
    std::vector<real_t> initial_T_PCR_pre;
    std::vector<real_t> initial_T_sero_neg;
    std::vector<real_t> initial_T_sero_pos;
    std::vector<real_t> initial_T_sero_pre;
    std::vector<real_t> initial_W_D_conf;
    std::vector<real_t> initial_W_D_unconf;
    std::vector<real_t> initial_W_R_conf;
    std::vector<real_t> initial_W_R_unconf;
    real_t initial_admit_conf_inc;
    real_t initial_beta_out;
    std::vector<real_t> initial_cum_admit_by_age;
    real_t initial_cum_admit_conf;
    real_t initial_cum_infections;
    std::vector<real_t> initial_cum_infections_per_strain;
    std::vector<real_t> initial_cum_n_E_vaccinated;
    std::vector<real_t> initial_cum_n_I_A_vaccinated;
    std::vector<real_t> initial_cum_n_I_P_vaccinated;
    std::vector<real_t> initial_cum_n_R_vaccinated;
    std::vector<real_t> initial_cum_n_S_vaccinated;
    std::vector<real_t> initial_cum_n_vaccinated;
    real_t initial_cum_new_conf;
    real_t initial_cum_sympt_cases;
    real_t initial_cum_sympt_cases_non_variant_over25;
    real_t initial_cum_sympt_cases_over25;
    real_t initial_general_tot;
    real_t initial_hosp_tot;
    real_t initial_new_conf_inc;
    std::vector<real_t> initial_prob_strain;
    real_t initial_react_pos;
    real_t initial_sero_pos;
    real_t initial_sympt_cases_inc;
    real_t initial_sympt_cases_non_variant_over25_inc;
    real_t initial_sympt_cases_over25_inc;
    real_t initial_time;
    std::vector<real_t> initial_tmp_vaccine_n_candidates;
    std::vector<real_t> initial_tmp_vaccine_probability;
    int k_A;
    int k_C_1;
    int k_C_2;
    int k_E;
    int k_G_D;
    int k_H_D;
    int k_H_R;
    int k_ICU_D;
    int k_ICU_W_D;
    int k_ICU_W_R;
    int k_ICU_pre;
    int k_P;
    int k_PCR_pos;
    int k_PCR_pre;
    int k_W_D;
    int k_W_R;
    int k_sero_pos;
    real_t kappa_ICU;
    real_t kappa_admitted;
    real_t kappa_all_admission;
    real_t kappa_death;
    real_t kappa_death_carehomes;
    real_t kappa_death_comm;
    real_t kappa_death_hosp;
    real_t kappa_death_non_hosp;
    real_t kappa_diagnoses;
    real_t kappa_general;
    real_t kappa_hosp;
    real_t kappa_pillar2_cases;
    std::vector<real_t> m;
    real_t model_pcr_and_serology;
    real_t model_pcr_and_serology_user;
    int n_age_groups;
    int n_doses;
    int n_groups;
    int n_strains;
    int n_vacc_classes;
    int offset_variable_E;
    int offset_variable_G_D;
    int offset_variable_H_D_conf;
    int offset_variable_H_D_unconf;
    int offset_variable_H_R_conf;
    int offset_variable_H_R_unconf;
    int offset_variable_ICU_D_conf;
    int offset_variable_ICU_D_unconf;
    int offset_variable_ICU_W_D_conf;
    int offset_variable_ICU_W_D_unconf;
    int offset_variable_ICU_W_R_conf;
    int offset_variable_ICU_W_R_unconf;
    int offset_variable_ICU_pre_conf;
    int offset_variable_ICU_pre_unconf;
    int offset_variable_I_A;
    int offset_variable_I_C_1;
    int offset_variable_I_C_2;
    int offset_variable_I_P;
    int offset_variable_I_weighted;
    int offset_variable_R;
    int offset_variable_S;
    int offset_variable_T_PCR_neg;
    int offset_variable_T_PCR_pos;
    int offset_variable_T_PCR_pre;
    int offset_variable_T_sero_neg;
    int offset_variable_T_sero_pos;
    int offset_variable_T_sero_pre;
    int offset_variable_W_D_conf;
    int offset_variable_W_D_unconf;
    int offset_variable_W_R_conf;
    int offset_variable_W_R_unconf;
    int offset_variable_cum_n_E_vaccinated;
    int offset_variable_cum_n_I_A_vaccinated;
    int offset_variable_cum_n_I_P_vaccinated;
    int offset_variable_cum_n_R_vaccinated;
    int offset_variable_cum_n_S_vaccinated;
    int offset_variable_cum_n_vaccinated;
    int offset_variable_prob_strain;
    int offset_variable_tmp_vaccine_probability;
    std::vector<real_t> p_C;
    real_t p_E_progress;
    real_t p_G_D_progress;
    std::vector<real_t> p_G_D_step;
    real_t p_H_D_progress;
    std::vector<real_t> p_H_D_step;
    real_t p_H_R_progress;
    std::vector<real_t> p_H_step;
    real_t p_ICU_D_progress;
    std::vector<real_t> p_ICU_D_step;
    real_t p_ICU_W_D_progress;
    real_t p_ICU_W_R_progress;
    real_t p_ICU_pre_progress;
    std::vector<real_t> p_ICU_step;
    real_t p_I_A_progress;
    real_t p_I_C_1_progress;
    real_t p_I_C_2_progress;
    real_t p_I_P_progress;
    real_t p_NC;
    std::vector<real_t> p_RS;
    real_t p_T_PCR_pos_progress;
    real_t p_T_PCR_pre_progress;
    real_t p_T_sero_pos_progress;
    std::vector<real_t> p_T_sero_pre_progress;
    real_t p_W_D_progress;
    std::vector<real_t> p_W_D_step;
    real_t p_W_R_progress;
    std::vector<real_t> p_sero_pos;
    real_t p_sero_pre_1;
    std::vector<real_t> p_star_step;
    real_t p_test;
    real_t phi_ICU;
    real_t phi_admitted;
    real_t phi_all_admission;
    real_t phi_death_carehomes;
    real_t phi_death_comm;
    real_t phi_death_hosp;
    real_t phi_diagnoses;
    real_t phi_general;
    real_t phi_hosp;
    real_t phi_pillar2_cases;
    real_t pillar2_sensitivity;
    real_t pillar2_specificity;
    std::vector<real_t> psi_G_D;
    std::vector<real_t> psi_H;
    std::vector<real_t> psi_H_D;
    std::vector<real_t> psi_ICU;
    std::vector<real_t> psi_ICU_D;
    std::vector<real_t> psi_W_D;
    std::vector<real_t> psi_star;
    real_t react_sensitivity;
    real_t react_specificity;
    std::vector<real_t> rel_infectivity;
    std::vector<real_t> rel_p_hosp_if_sympt;
    std::vector<real_t> rel_p_sympt;
    std::vector<real_t> rel_susceptibility;
    real_t rho_pillar2_tests;
    real_t sero_sensitivity;
    real_t sero_specificity;
    int steps_per_day;
    std::vector<real_t> strain_seed_step;
    std::vector<real_t> strain_transmission;
    std::vector<real_t> vaccine_dose_step;
    std::vector<real_t> vaccine_progression_rate_base;
    std::vector<real_t> waning_rate;
  };
  struct internal_t {
    std::vector<real_t> I_weighted_strain;
    std::vector<real_t> I_with_diff_trans;
    std::vector<real_t> aux_E;
    std::vector<real_t> aux_G_D;
    std::vector<real_t> aux_H_D_conf;
    std::vector<real_t> aux_H_D_unconf;
    std::vector<real_t> aux_H_R_conf;
    std::vector<real_t> aux_H_R_unconf;
    std::vector<real_t> aux_ICU_D_conf;
    std::vector<real_t> aux_ICU_D_unconf;
    std::vector<real_t> aux_ICU_W_D_conf;
    std::vector<real_t> aux_ICU_W_D_unconf;
    std::vector<real_t> aux_ICU_W_R_conf;
    std::vector<real_t> aux_ICU_W_R_unconf;
    std::vector<real_t> aux_ICU_pre_conf;
    std::vector<real_t> aux_ICU_pre_unconf;
    std::vector<real_t> aux_I_A;
    std::vector<real_t> aux_I_C_1;
    std::vector<real_t> aux_I_C_2;
    std::vector<real_t> aux_I_P;
    std::vector<real_t> aux_W_D_conf;
    std::vector<real_t> aux_W_D_unconf;
    std::vector<real_t> aux_W_R_conf;
    std::vector<real_t> aux_W_R_unconf;
    std::vector<real_t> delta_D_hosp;
    std::vector<real_t> delta_D_non_hosp;
    std::vector<real_t> lambda;
    std::vector<real_t> n_EE;
    std::vector<real_t> n_EE_next_vacc_class;
    std::vector<real_t> n_EI_A;
    std::vector<real_t> n_EI_A_next_vacc_class;
    std::vector<real_t> n_EI_P;
    std::vector<real_t> n_EI_P_next_vacc_class;
    std::vector<real_t> n_E_next_vacc_class;
    std::vector<real_t> n_E_progress;
    std::vector<real_t> n_G_D_progress;
    std::vector<real_t> n_H_D_conf_progress;
    std::vector<real_t> n_H_D_unconf_progress;
    std::vector<real_t> n_H_D_unconf_to_conf;
    std::vector<real_t> n_H_R_conf_progress;
    std::vector<real_t> n_H_R_unconf_progress;
    std::vector<real_t> n_H_R_unconf_to_conf;
    std::vector<real_t> n_ICU_D_conf_progress;
    std::vector<real_t> n_ICU_D_unconf_progress;
    std::vector<real_t> n_ICU_D_unconf_to_conf;
    std::vector<real_t> n_ICU_W_D_conf_progress;
    std::vector<real_t> n_ICU_W_D_unconf_progress;
    std::vector<real_t> n_ICU_W_D_unconf_to_conf;
    std::vector<real_t> n_ICU_W_R_conf_progress;
    std::vector<real_t> n_ICU_W_R_unconf_progress;
    std::vector<real_t> n_ICU_W_R_unconf_to_conf;
    std::vector<real_t> n_ICU_pre_conf_progress;
    std::vector<real_t> n_ICU_pre_conf_to_ICU_D_conf;
    std::vector<real_t> n_ICU_pre_conf_to_ICU_W_D_conf;
    std::vector<real_t> n_ICU_pre_conf_to_ICU_W_R_conf;
    std::vector<real_t> n_ICU_pre_unconf_progress;
    std::vector<real_t> n_ICU_pre_unconf_to_ICU_D_unconf;
    std::vector<real_t> n_ICU_pre_unconf_to_ICU_W_D_unconf;
    std::vector<real_t> n_ICU_pre_unconf_to_ICU_W_R_unconf;
    std::vector<real_t> n_ICU_pre_unconf_to_conf;
    std::vector<real_t> n_II_A;
    std::vector<real_t> n_II_A_next_vacc_class;
    std::vector<real_t> n_II_P;
    std::vector<real_t> n_II_P_next_vacc_class;
    std::vector<real_t> n_I_A_next_vacc_class;
    std::vector<real_t> n_I_A_progress;
    std::vector<real_t> n_I_C_1_progress;
    std::vector<real_t> n_I_C_2_progress;
    std::vector<real_t> n_I_C_2_to_G_D;
    std::vector<real_t> n_I_C_2_to_H_D;
    std::vector<real_t> n_I_C_2_to_H_D_conf;
    std::vector<real_t> n_I_C_2_to_H_R;
    std::vector<real_t> n_I_C_2_to_H_R_conf;
    std::vector<real_t> n_I_C_2_to_ICU_pre;
    std::vector<real_t> n_I_C_2_to_ICU_pre_conf;
    std::vector<real_t> n_I_C_2_to_R;
    std::vector<real_t> n_I_C_2_to_hosp;
    std::vector<real_t> n_I_P_next_vacc_class;
    std::vector<real_t> n_I_P_progress;
    std::vector<real_t> n_RS;
    std::vector<real_t> n_RS_next_vacc_class;
    std::vector<real_t> n_R_next_vacc_class;
    std::vector<real_t> n_R_next_vacc_class_capped;
    std::vector<real_t> n_R_next_vacc_class_tmp;
    std::vector<real_t> n_R_progress;
    std::vector<real_t> n_R_progress_capped;
    std::vector<real_t> n_R_progress_tmp;
    std::vector<real_t> n_SE;
    std::vector<real_t> n_SE_next_vacc_class;
    std::vector<real_t> n_S_next_vacc_class;
    std::vector<real_t> n_S_progress;
    std::vector<real_t> n_S_progress_tot;
    std::vector<real_t> n_T_PCR_pos_progress;
    std::vector<real_t> n_T_PCR_pre_progress;
    std::vector<real_t> n_T_sero_pos_progress;
    std::vector<real_t> n_T_sero_pre_progress;
    std::vector<real_t> n_T_sero_pre_to_T_sero_pos;
    std::vector<real_t> n_W_D_conf_progress;
    std::vector<real_t> n_W_D_unconf_progress;
    std::vector<real_t> n_W_D_unconf_to_conf;
    std::vector<real_t> n_W_R_conf_progress;
    std::vector<real_t> n_W_R_unconf_progress;
    std::vector<real_t> n_W_R_unconf_to_conf;
    std::vector<real_t> n_com_to_T_sero_pre;
    std::vector<real_t> n_hosp_non_ICU;
    std::vector<real_t> new_E;
    std::vector<real_t> new_G_D;
    std::vector<real_t> new_H_D_conf;
    std::vector<real_t> new_H_D_unconf;
    std::vector<real_t> new_H_R_conf;
    std::vector<real_t> new_H_R_unconf;
    std::vector<real_t> new_ICU_D_conf;
    std::vector<real_t> new_ICU_D_unconf;
    std::vector<real_t> new_ICU_W_D_conf;
    std::vector<real_t> new_ICU_W_D_unconf;
    std::vector<real_t> new_ICU_W_R_conf;
    std::vector<real_t> new_ICU_W_R_unconf;
    std::vector<real_t> new_ICU_pre_conf;
    std::vector<real_t> new_ICU_pre_unconf;
    std::vector<real_t> new_I_A;
    std::vector<real_t> new_I_C_1;
    std::vector<real_t> new_I_C_2;
    std::vector<real_t> new_I_P;
    std::vector<real_t> new_R;
    std::vector<real_t> new_S;
    std::vector<real_t> new_T_PCR_neg;
    std::vector<real_t> new_T_PCR_pos;
    std::vector<real_t> new_T_PCR_pre;
    std::vector<real_t> new_T_sero_neg;
    std::vector<real_t> new_T_sero_pos;
    std::vector<real_t> new_T_sero_pre;
    std::vector<real_t> new_W_D_conf;
    std::vector<real_t> new_W_D_unconf;
    std::vector<real_t> new_W_R_conf;
    std::vector<real_t> new_W_R_unconf;
    std::vector<real_t> p_E_next_vacc_class;
    std::vector<real_t> p_G_D_by_age;
    std::vector<real_t> p_H_D_by_age;
    std::vector<real_t> p_H_by_age;
    std::vector<real_t> p_ICU_D_by_age;
    std::vector<real_t> p_ICU_by_age;
    std::vector<real_t> p_I_A_next_vacc_class;
    std::vector<real_t> p_I_P_next_vacc_class;
    std::vector<real_t> p_R_next_vacc_class;
    std::vector<real_t> p_SE;
    std::vector<real_t> p_S_next_vacc_class;
    std::vector<real_t> p_W_D_by_age;
    std::vector<real_t> p_star_by_age;
    std::vector<real_t> s_ij;
    std::vector<real_t> vaccine_n_candidates;
    std::vector<real_t> vaccine_probability;
    std::vector<real_t> vaccine_probability_doses;
  };
  carehomes(const dust::pars_t<carehomes>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() {
    return shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_T_sero_pre + shared->dim_T_sero_pos + shared->dim_T_PCR_pre + shared->dim_T_PCR_pos + 141;
  }
  std::vector<real_t> initial(size_t step) {
    std::vector<real_t> state(shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_T_sero_pre + shared->dim_T_sero_pos + shared->dim_T_PCR_pre + shared->dim_T_PCR_pos + 141);
    state[0] = shared->initial_time;
    state[1] = shared->initial_admit_conf_inc;
    state[2] = shared->initial_new_conf_inc;
    state[3] = shared->initial_cum_infections;
    state[4] = shared->initial_cum_admit_conf;
    state[5] = shared->initial_cum_new_conf;
    state[6] = shared->initial_beta_out;
    state[7] = shared->initial_N_tot2;
    state[8] = shared->initial_N_tot3;
    state[9] = shared->initial_ICU_tot;
    state[10] = shared->initial_general_tot;
    state[11] = shared->initial_hosp_tot;
    state[12] = shared->initial_D_hosp_tot;
    state[13] = shared->initial_D_comm_tot;
    state[14] = shared->initial_D_comm_inc;
    state[15] = shared->initial_D_carehomes_tot;
    state[16] = shared->initial_D_carehomes_inc;
    state[17] = shared->initial_D_hosp_inc;
    state[18] = shared->initial_D_tot;
    state[19] = shared->initial_sero_pos;
    state[20] = shared->initial_cum_sympt_cases;
    state[21] = shared->initial_cum_sympt_cases_over25;
    state[22] = shared->initial_cum_sympt_cases_non_variant_over25;
    state[23] = shared->initial_sympt_cases_inc;
    state[24] = shared->initial_sympt_cases_over25_inc;
    state[25] = shared->initial_sympt_cases_non_variant_over25_inc;
    state[26] = shared->initial_react_pos;
    std::copy(shared->initial_D_hosp.begin(), shared->initial_D_hosp.end(), state.begin() + 27);
    std::copy(shared->initial_D_non_hosp.begin(), shared->initial_D_non_hosp.end(), state.begin() + 46);
    std::copy(shared->initial_cum_admit_by_age.begin(), shared->initial_cum_admit_by_age.end(), state.begin() + 65);
    std::copy(shared->initial_N_tot.begin(), shared->initial_N_tot.end(), state.begin() + 84);
    std::copy(shared->initial_tmp_vaccine_n_candidates.begin(), shared->initial_tmp_vaccine_n_candidates.end(), state.begin() + 103);
    std::copy(shared->initial_cum_infections_per_strain.begin(), shared->initial_cum_infections_per_strain.end(), state.begin() + 141);
    std::copy(shared->initial_cum_n_S_vaccinated.begin(), shared->initial_cum_n_S_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_S_vaccinated);
    std::copy(shared->initial_cum_n_E_vaccinated.begin(), shared->initial_cum_n_E_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_E_vaccinated);
    std::copy(shared->initial_cum_n_I_A_vaccinated.begin(), shared->initial_cum_n_I_A_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_I_A_vaccinated);
    std::copy(shared->initial_cum_n_I_P_vaccinated.begin(), shared->initial_cum_n_I_P_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_I_P_vaccinated);
    std::copy(shared->initial_cum_n_R_vaccinated.begin(), shared->initial_cum_n_R_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_R_vaccinated);
    std::copy(shared->initial_cum_n_vaccinated.begin(), shared->initial_cum_n_vaccinated.end(), state.begin() + shared->offset_variable_cum_n_vaccinated);
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + shared->offset_variable_S);
    std::copy(shared->initial_prob_strain.begin(), shared->initial_prob_strain.end(), state.begin() + shared->offset_variable_prob_strain);
    std::copy(shared->initial_I_weighted.begin(), shared->initial_I_weighted.end(), state.begin() + shared->offset_variable_I_weighted);
    std::copy(shared->initial_tmp_vaccine_probability.begin(), shared->initial_tmp_vaccine_probability.end(), state.begin() + shared->offset_variable_tmp_vaccine_probability);
    std::copy(shared->initial_T_sero_neg.begin(), shared->initial_T_sero_neg.end(), state.begin() + shared->offset_variable_T_sero_neg);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_T_PCR_neg.begin(), shared->initial_T_PCR_neg.end(), state.begin() + shared->offset_variable_T_PCR_neg);
    std::copy(shared->initial_E.begin(), shared->initial_E.end(), state.begin() + shared->offset_variable_E);
    std::copy(shared->initial_I_A.begin(), shared->initial_I_A.end(), state.begin() + shared->offset_variable_I_A);
    std::copy(shared->initial_I_P.begin(), shared->initial_I_P.end(), state.begin() + shared->offset_variable_I_P);
    std::copy(shared->initial_I_C_1.begin(), shared->initial_I_C_1.end(), state.begin() + shared->offset_variable_I_C_1);
    std::copy(shared->initial_I_C_2.begin(), shared->initial_I_C_2.end(), state.begin() + shared->offset_variable_I_C_2);
    std::copy(shared->initial_G_D.begin(), shared->initial_G_D.end(), state.begin() + shared->offset_variable_G_D);
    std::copy(shared->initial_ICU_pre_unconf.begin(), shared->initial_ICU_pre_unconf.end(), state.begin() + shared->offset_variable_ICU_pre_unconf);
    std::copy(shared->initial_ICU_pre_conf.begin(), shared->initial_ICU_pre_conf.end(), state.begin() + shared->offset_variable_ICU_pre_conf);
    std::copy(shared->initial_H_R_unconf.begin(), shared->initial_H_R_unconf.end(), state.begin() + shared->offset_variable_H_R_unconf);
    std::copy(shared->initial_H_R_conf.begin(), shared->initial_H_R_conf.end(), state.begin() + shared->offset_variable_H_R_conf);
    std::copy(shared->initial_H_D_unconf.begin(), shared->initial_H_D_unconf.end(), state.begin() + shared->offset_variable_H_D_unconf);
    std::copy(shared->initial_H_D_conf.begin(), shared->initial_H_D_conf.end(), state.begin() + shared->offset_variable_H_D_conf);
    std::copy(shared->initial_ICU_W_R_unconf.begin(), shared->initial_ICU_W_R_unconf.end(), state.begin() + shared->offset_variable_ICU_W_R_unconf);
    std::copy(shared->initial_ICU_W_R_conf.begin(), shared->initial_ICU_W_R_conf.end(), state.begin() + shared->offset_variable_ICU_W_R_conf);
    std::copy(shared->initial_ICU_W_D_unconf.begin(), shared->initial_ICU_W_D_unconf.end(), state.begin() + shared->offset_variable_ICU_W_D_unconf);
    std::copy(shared->initial_ICU_W_D_conf.begin(), shared->initial_ICU_W_D_conf.end(), state.begin() + shared->offset_variable_ICU_W_D_conf);
    std::copy(shared->initial_ICU_D_unconf.begin(), shared->initial_ICU_D_unconf.end(), state.begin() + shared->offset_variable_ICU_D_unconf);
    std::copy(shared->initial_ICU_D_conf.begin(), shared->initial_ICU_D_conf.end(), state.begin() + shared->offset_variable_ICU_D_conf);
    std::copy(shared->initial_W_R_unconf.begin(), shared->initial_W_R_unconf.end(), state.begin() + shared->offset_variable_W_R_unconf);
    std::copy(shared->initial_W_R_conf.begin(), shared->initial_W_R_conf.end(), state.begin() + shared->offset_variable_W_R_conf);
    std::copy(shared->initial_W_D_unconf.begin(), shared->initial_W_D_unconf.end(), state.begin() + shared->offset_variable_W_D_unconf);
    std::copy(shared->initial_W_D_conf.begin(), shared->initial_W_D_conf.end(), state.begin() + shared->offset_variable_W_D_conf);
    std::copy(shared->initial_T_sero_pre.begin(), shared->initial_T_sero_pre.end(), state.begin() + shared->offset_variable_T_sero_pre);
    std::copy(shared->initial_T_sero_pos.begin(), shared->initial_T_sero_pos.end(), state.begin() + shared->offset_variable_T_sero_pos);
    std::copy(shared->initial_T_PCR_pre.begin(), shared->initial_T_PCR_pre.end(), state.begin() + shared->offset_variable_T_PCR_pre);
    std::copy(shared->initial_T_PCR_pos.begin(), shared->initial_T_PCR_pos.end(), state.begin() + shared->offset_variable_T_PCR_pos);
    return state;
  }
  HOST void update(size_t step, const real_t * state, dust::rng_state_t<real_t>& rng_state, real_t * state_next) {
    const real_t * cum_n_S_vaccinated = state + shared->offset_variable_cum_n_S_vaccinated;
    const real_t * cum_n_E_vaccinated = state + shared->offset_variable_cum_n_E_vaccinated;
    const real_t * cum_n_I_A_vaccinated = state + shared->offset_variable_cum_n_I_A_vaccinated;
    const real_t * cum_n_I_P_vaccinated = state + shared->offset_variable_cum_n_I_P_vaccinated;
    const real_t * cum_n_R_vaccinated = state + shared->offset_variable_cum_n_R_vaccinated;
    const real_t * S = state + shared->offset_variable_S;
    const real_t * E = state + shared->offset_variable_E;
    const real_t * I_A = state + shared->offset_variable_I_A;
    const real_t * I_P = state + shared->offset_variable_I_P;
    const real_t * I_C_1 = state + shared->offset_variable_I_C_1;
    const real_t * I_C_2 = state + shared->offset_variable_I_C_2;
    const real_t * G_D = state + shared->offset_variable_G_D;
    const real_t * ICU_pre_unconf = state + shared->offset_variable_ICU_pre_unconf;
    const real_t * ICU_pre_conf = state + shared->offset_variable_ICU_pre_conf;
    const real_t * H_R_unconf = state + shared->offset_variable_H_R_unconf;
    const real_t * H_R_conf = state + shared->offset_variable_H_R_conf;
    const real_t * H_D_unconf = state + shared->offset_variable_H_D_unconf;
    const real_t * H_D_conf = state + shared->offset_variable_H_D_conf;
    const real_t * ICU_W_R_unconf = state + shared->offset_variable_ICU_W_R_unconf;
    const real_t * ICU_W_R_conf = state + shared->offset_variable_ICU_W_R_conf;
    const real_t * ICU_W_D_unconf = state + shared->offset_variable_ICU_W_D_unconf;
    const real_t * ICU_W_D_conf = state + shared->offset_variable_ICU_W_D_conf;
    const real_t * ICU_D_unconf = state + shared->offset_variable_ICU_D_unconf;
    const real_t * ICU_D_conf = state + shared->offset_variable_ICU_D_conf;
    const real_t * W_R_unconf = state + shared->offset_variable_W_R_unconf;
    const real_t * W_R_conf = state + shared->offset_variable_W_R_conf;
    const real_t * W_D_unconf = state + shared->offset_variable_W_D_unconf;
    const real_t * W_D_conf = state + shared->offset_variable_W_D_conf;
    const real_t * T_sero_pre = state + shared->offset_variable_T_sero_pre;
    const real_t * T_sero_pos = state + shared->offset_variable_T_sero_pos;
    const real_t * T_sero_neg = state + shared->offset_variable_T_sero_neg;
    const real_t * R = state + shared->offset_variable_R;
    const real_t * D_hosp = state + 27;
    const real_t * D_non_hosp = state + 46;
    const real_t * T_PCR_pre = state + shared->offset_variable_T_PCR_pre;
    const real_t * T_PCR_pos = state + shared->offset_variable_T_PCR_pos;
    const real_t * T_PCR_neg = state + shared->offset_variable_T_PCR_neg;
    const real_t cum_admit_conf = state[4];
    const real_t cum_new_conf = state[5];
    const real_t admit_conf_inc = state[1];
    const real_t new_conf_inc = state[2];
    const real_t * cum_admit_by_age = state + 65;
    const real_t cum_infections = state[3];
    const real_t * cum_infections_per_strain = state + 141;
    const real_t D_hosp_tot = state[12];
    const real_t D_comm_tot = state[13];
    const real_t D_comm_inc = state[14];
    const real_t D_carehomes_tot = state[15];
    const real_t D_carehomes_inc = state[16];
    const real_t D_hosp_inc = state[17];
    const real_t D_tot = state[18];
    const real_t cum_sympt_cases = state[20];
    const real_t cum_sympt_cases_over25 = state[21];
    const real_t cum_sympt_cases_non_variant_over25 = state[22];
    const real_t sympt_cases_inc = state[23];
    const real_t sympt_cases_over25_inc = state[24];
    const real_t sympt_cases_non_variant_over25_inc = state[25];
    state_next[7] = odin_sum1<real_t>(S, 0, shared->dim_cum_n_S_vaccinated) + odin_sum1<real_t>(T_sero_pre, 0, shared->dim_T_sero_pre) + odin_sum1<real_t>(T_sero_pos, 0, shared->dim_T_sero_pos) + odin_sum1<real_t>(T_sero_neg, 0, shared->dim_R) + odin_sum1<real_t>(E, 0, shared->dim_E);
    state_next[8] = odin_sum1<real_t>(S, 0, shared->dim_cum_n_S_vaccinated) + odin_sum1<real_t>(T_PCR_pre, 0, shared->dim_T_PCR_pre) + odin_sum1<real_t>(T_PCR_pos, 0, shared->dim_T_PCR_pos) + odin_sum1<real_t>(T_PCR_neg, 0, shared->dim_R);
    real_t beta = (static_cast<int>(step) >= shared->dim_beta_step ? shared->beta_step[shared->dim_beta_step - 1] : shared->beta_step[step + 1 - 1]);
    real_t p_G_D = (static_cast<int>(step) >= shared->dim_p_G_D_step ? shared->p_G_D_step[shared->dim_p_G_D_step - 1] : shared->p_G_D_step[step + 1 - 1]);
    real_t p_H = (static_cast<int>(step) >= shared->dim_p_H_step ? shared->p_H_step[shared->dim_p_H_step - 1] : shared->p_H_step[step + 1 - 1]);
    real_t p_H_D = (static_cast<int>(step) >= shared->dim_p_H_D_step ? shared->p_H_D_step[shared->dim_p_H_D_step - 1] : shared->p_H_D_step[step + 1 - 1]);
    real_t p_ICU = (static_cast<int>(step) >= shared->dim_p_ICU_step ? shared->p_ICU_step[shared->dim_p_ICU_step - 1] : shared->p_ICU_step[step + 1 - 1]);
    real_t p_ICU_D = (static_cast<int>(step) >= shared->dim_p_ICU_D_step ? shared->p_ICU_D_step[shared->dim_p_ICU_D_step - 1] : shared->p_ICU_D_step[step + 1 - 1]);
    real_t p_W_D = (static_cast<int>(step) >= shared->dim_p_W_D_step ? shared->p_W_D_step[shared->dim_p_W_D_step - 1] : shared->p_W_D_step[step + 1 - 1]);
    real_t p_star = (static_cast<int>(step) >= shared->dim_p_star_step ? shared->p_star_step[shared->dim_p_star_step - 1] : shared->p_star_step[step + 1 - 1]);
    real_t strain_seed = ((static_cast<int>(step) >= shared->dim_strain_seed_step ? shared->strain_seed_step[shared->dim_strain_seed_step - 1] : shared->strain_seed_step[step + 1 - 1]));
    state_next[0] = (step + 1) * shared->dt;
    for (int i = 1; i <= 19; ++i) {
      internal.p_G_D_by_age[i - 1] = p_G_D * shared->psi_G_D[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_H_D_by_age[i - 1] = p_H_D * shared->psi_H_D[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_H_by_age[i - 1] = p_H * shared->psi_H[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_ICU_D_by_age[i - 1] = p_ICU_D * shared->psi_ICU_D[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_ICU_by_age[i - 1] = p_ICU * shared->psi_ICU[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_W_D_by_age[i - 1] = p_W_D * shared->psi_W_D[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      internal.p_star_by_age[i - 1] = p_star * shared->psi_star[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      state_next[84 + i - 1] = odin_sum2<real_t>(S, i - 1, i, 0, shared->dim_rel_susceptibility_2, 19) + odin_sum3<real_t>(R, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12) + D_hosp[i - 1] + odin_sum4<real_t>(E, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_E, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_E_123) + odin_sum4<real_t>(I_A, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_A, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_I_A_123) + odin_sum4<real_t>(I_P, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_P, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_I_P_123) + odin_sum4<real_t>(I_C_1, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_C_1, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_I_C_1_123) + odin_sum4<real_t>(I_C_2, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_C_2, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_I_C_2_123) + odin_sum4<real_t>(ICU_pre_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_pre, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(ICU_pre_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_pre, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(H_R_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_H_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(H_R_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_H_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(H_D_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_H_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(H_D_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_H_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(ICU_W_R_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_W_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_R_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_W_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_D_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_W_D_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_D_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(ICU_D_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_ICU_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(W_R_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_W_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_R_unconf_123) + odin_sum4<real_t>(W_R_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_W_R, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_R_unconf_123) + odin_sum4<real_t>(W_D_conf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_D_unconf_123) + odin_sum4<real_t>(W_D_unconf, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_D_unconf_123) + odin_sum4<real_t>(G_D, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_G_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_G_D_123) + D_non_hosp[i - 1];
    }
    state_next[6] = beta;
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_S_vaccinated[19 * (j - 1) + i - 1] + cum_n_E_vaccinated[19 * (j - 1) + i - 1] + cum_n_I_A_vaccinated[19 * (j - 1) + i - 1] + cum_n_I_P_vaccinated[19 * (j - 1) + i - 1] + cum_n_R_vaccinated[19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= 2; ++j) {
        internal.vaccine_n_candidates[i - 1 + 19 * (j - 1)] = S[19 * (shared->index_dose[j - 1] - 1) + i - 1] + odin_sum4<real_t>(E, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_E, shared->index_dose[j - 1] - 1, shared->index_dose[j - 1], 19, shared->dim_E_12, shared->dim_E_123) + odin_sum4<real_t>(I_A, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_A, shared->index_dose[j - 1] - 1, shared->index_dose[j - 1], 19, shared->dim_E_12, shared->dim_I_A_123) + odin_sum4<real_t>(I_P, i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_P, shared->index_dose[j - 1] - 1, shared->index_dose[j - 1], 19, shared->dim_E_12, shared->dim_I_P_123) + odin_sum3<real_t>(R, i - 1, i, 0, shared->dim_strain_transmission, shared->index_dose[j - 1] - 1, shared->index_dose[j - 1], 19, shared->dim_E_12);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.I_with_diff_trans[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = shared->rel_infectivity[19 * (k - 1) + i - 1] * shared->strain_transmission[j - 1] * (shared->I_A_transmission * odin_sum4<real_t>(I_A, i - 1, i, j - 1, j, 0, shared->k_A, k - 1, k, 19, shared->dim_E_12, shared->dim_I_A_123) + shared->I_P_transmission * odin_sum4<real_t>(I_P, i - 1, i, j - 1, j, 0, shared->k_P, k - 1, k, 19, shared->dim_E_12, shared->dim_I_P_123) + shared->I_C_1_transmission * odin_sum4<real_t>(I_C_1, i - 1, i, j - 1, j, 0, shared->k_C_1, k - 1, k, 19, shared->dim_E_12, shared->dim_I_C_1_123) + shared->I_C_2_transmission * odin_sum4<real_t>(I_C_2, i - 1, i, j - 1, j, 0, shared->k_C_2, k - 1, k, 19, shared->dim_E_12, shared->dim_I_C_2_123) + shared->hosp_transmission * (odin_sum4<real_t>(ICU_pre_unconf, i - 1, i, j - 1, j, 0, shared->k_ICU_pre, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(ICU_pre_conf, i - 1, i, j - 1, j, 0, shared->k_ICU_pre, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(H_R_unconf, i - 1, i, j - 1, j, 0, shared->k_H_R, k - 1, k, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(H_R_conf, i - 1, i, j - 1, j, 0, shared->k_H_R, k - 1, k, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(H_D_unconf, i - 1, i, j - 1, j, 0, shared->k_H_D, k - 1, k, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(H_D_conf, i - 1, i, j - 1, j, 0, shared->k_H_D, k - 1, k, 19, shared->dim_E_12, shared->dim_H_D_unconf_123)) + shared->ICU_transmission * (odin_sum4<real_t>(ICU_W_R_unconf, i - 1, i, j - 1, j, 0, shared->k_ICU_W_R, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_R_conf, i - 1, i, j - 1, j, 0, shared->k_ICU_W_R, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_D_unconf, i - 1, i, j - 1, j, 0, shared->k_ICU_W_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_W_D_conf, i - 1, i, j - 1, j, 0, shared->k_ICU_W_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_D_unconf, i - 1, i, j - 1, j, 0, shared->k_ICU_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(ICU_D_conf, i - 1, i, j - 1, j, 0, shared->k_ICU_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123)) + shared->G_D_transmission * odin_sum4<real_t>(G_D, i - 1, i, j - 1, j, 0, shared->k_G_D, k - 1, k, 19, shared->dim_E_12, shared->dim_G_D_123));
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_E_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_E_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_G_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_G_D_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(G_D[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_G_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_D_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_H_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_D_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_H_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_R_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_H_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_R_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_H_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_D_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_D_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_D_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_W_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_D_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_W_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_R_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_W_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_R_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_W_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_pre_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_pre_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_pre_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_ICU_pre_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_A_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_I_A_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_1; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_C_1_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_C_1[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_I_C_1_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_C_2_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_C_2[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_I_C_2_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_P_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_I_P_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_progress_tmp[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_RS[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_T_PCR_pos_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_PCR_pos[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_T_PCR_pos_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_T_PCR_pre_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_PCR_pre[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_T_PCR_pre_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_sero_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_T_sero_pos_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_sero_pos[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_T_sero_pos_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_D_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_W_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_D_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_W_D_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_R_conf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_W_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_R_unconf_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_W_R_progress);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= 2; ++j) {
        state_next[103 + i - 1 + 19 * (j - 1)] = internal.vaccine_n_candidates[19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= 2; ++j) {
        internal.vaccine_probability_doses[i - 1 + 19 * (j - 1)] = ((static_cast<int>(step) >= shared->dim_vaccine_dose_step_3 || internal.vaccine_n_candidates[19 * (j - 1) + i - 1] == 0 ? 0 : shared->vaccine_dose_step[shared->dim_vaccine_dose_step_12 * (step + 1 - 1) + shared->dim_vaccine_dose_step_1 * (j - 1) + i - 1] / (real_t) internal.vaccine_n_candidates[19 * (j - 1) + i - 1]));
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_D_conf_progress[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_D_conf_progress[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_D_unconf_progress[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_D_unconf_progress[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_R_conf_progress[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_R_conf_progress[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_R_unconf_progress[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_R_unconf_progress[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_conf_progress[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_conf_progress[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_unconf_progress[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_unconf_progress[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_I_C_2[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = internal.n_I_C_1_progress[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (shared->k_C_1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_C_2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_C_2[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = internal.n_I_C_2_progress[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_C_2[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = internal.aux_I_C_2[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_progress[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_D_conf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (shared->k_ICU_W_D - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_D_conf_progress[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_D_conf_progress[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_D_unconf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (shared->k_ICU_W_D - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_D_unconf_progress[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_D_unconf_progress[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_R_conf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (shared->k_ICU_W_R - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_R_conf_progress[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_R_conf_progress[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_R_unconf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (shared->k_ICU_W_R - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_R_unconf_progress[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_R_unconf_progress[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      internal.delta_D_hosp[i - 1] = odin_sum4<real_t>(internal.n_H_D_unconf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_H_D - 1, shared->k_H_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(internal.n_H_D_conf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_H_D - 1, shared->k_H_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(internal.n_ICU_D_unconf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_ICU_D - 1, shared->k_ICU_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(internal.n_ICU_D_conf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_ICU_D - 1, shared->k_ICU_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(internal.n_W_D_unconf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_W_D - 1, shared->k_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_D_unconf_123) + odin_sum4<real_t>(internal.n_W_D_conf_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_W_D - 1, shared->k_W_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_W_D_unconf_123);
    }
    for (int i = 1; i <= 19; ++i) {
      internal.delta_D_non_hosp[i - 1] = odin_sum4<real_t>(internal.n_G_D_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->k_G_D - 1, shared->k_G_D, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_G_D_123);
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_conf_to_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_ICU_pre_conf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1]), internal.p_ICU_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_unconf_to_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_ICU_pre_unconf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1]), internal.p_ICU_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_progress[shared->dim_I_C_2_123 * (k - 1) + shared->dim_E_12 * (shared->k_C_2 - 1) + 19 * (j - 1) + i - 1]), 1 - internal.p_H_by_age[i - 1] * shared->rel_p_hosp_if_sympt[19 * (k - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_progress_capped[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = std::min(internal.n_R_progress_tmp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], std::min(T_sero_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], T_PCR_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]));
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= 2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_T_sero_pre_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_sero_pre[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_T_sero_pre_progress[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_PCR_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = T_PCR_pos[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_T_PCR_pos_progress[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_T_PCR_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = internal.new_T_PCR_pos[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_T_PCR_pre_progress[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (shared->k_PCR_pre - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_PCR_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_PCR_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = internal.new_T_PCR_pos[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_T_PCR_pos_progress[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= 19; ++j) {
        for (int k = 1; k <= shared->dim_strain_transmission; ++k) {
          internal.s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = shared->m[19 * (j - 1) + i - 1] * odin_sum3<real_t>(internal.I_with_diff_trans.data(), j - 1, j, k - 1, k, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12);
        }
      }
    }
    for (int i = 1; i <= shared->n_age_groups; ++i) {
      for (int j = 1; j <= shared->n_groups; ++j) {
        for (int k = 1; k <= shared->dim_strain_transmission; ++k) {
          internal.s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = beta * internal.s_ij[361 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = (shared->n_age_groups + 1); i <= shared->n_groups; ++i) {
      for (int j = 1; j <= shared->n_age_groups; ++j) {
        for (int k = 1; k <= shared->dim_strain_transmission; ++k) {
          internal.s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = beta * internal.s_ij[361 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.vaccine_probability[i - 1 + 19 * (j - 1)] = 1 - std::exp(- shared->vaccine_progression_rate_base[19 * (j - 1) + i - 1] * shared->dt);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      int j = shared->index_dose[0];
      internal.vaccine_probability[i - 1 + 19 * (j - 1)] = internal.vaccine_probability_doses[19 * 0 + i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      int j = shared->index_dose[1];
      internal.vaccine_probability[i - 1 + 19 * (j - 1)] = internal.vaccine_probability_doses[19 * 1 + i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_conf_to_ICU_D_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_D_conf_progress[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_D_conf_progress[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_unconf_to_ICU_D_unconf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_D_unconf_progress[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_D_unconf_progress[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    real_t delta_D_carehomes_tot = internal.delta_D_non_hosp[18];
    real_t delta_D_comm_tot = odin_sum1<real_t>(internal.delta_D_non_hosp.data(), 0, 18);
    real_t delta_D_hosp_tot = odin_sum1<real_t>(internal.delta_D_hosp.data(), 0, 19);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        internal.lambda[i - 1 + 19 * (j - 1)] = odin_sum3<real_t>(internal.s_ij.data(), i - 1, i, 0, 19, j - 1, j, 19, 361);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_D_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_H_R_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_conf_to_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_ICU_pre_conf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_conf_to_ICU_D_conf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_W_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_unconf_to_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_ICU_pre_unconf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_unconf_to_ICU_D_unconf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_W_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_pre_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_progress[shared->dim_I_C_2_123 * (k - 1) + shared->dim_E_12 * (shared->k_C_2 - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_G_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = (shared->model_pcr_and_serology == 1 ? internal.n_R_progress_capped[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] : internal.n_R_progress_tmp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_T_sero_pre_to_T_sero_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(odin_sum4<real_t>(internal.n_T_sero_pre_progress.data(), i - 1, i, j - 1, j, 0, 2, k - 1, k, 19, shared->dim_E_12, shared->dim_T_sero_pre_123)), shared->p_sero_pos[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_D_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_W_R_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_I_C_2[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = I_C_2[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_I_C_2[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.p_E_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.vaccine_probability[19 * (l - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.p_I_A_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.vaccine_probability[19 * (l - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.p_I_P_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.vaccine_probability[19 * (l - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.p_R_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.vaccine_probability[19 * (k - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.p_S_next_vacc_class[i - 1 + 19 * (j - 1)] = internal.vaccine_probability[19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      state_next[27 + i - 1] = D_hosp[i - 1] + internal.delta_D_hosp[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      state_next[46 + i - 1] = D_non_hosp[i - 1] + internal.delta_D_non_hosp[i - 1];
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_T_PCR_pos + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = internal.new_T_PCR_pos[shared->dim_T_PCR_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    state_next[26] = odin_sum4<real_t>(internal.new_T_PCR_pos.data(), 1, 18, 0, shared->dim_strain_transmission, 0, shared->k_PCR_pos, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_T_PCR_pos_123);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_tmp_vaccine_probability + i - 1 + 19 * (j - 1)] = internal.vaccine_probability[19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = internal.n_I_C_2_to_G_D[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_G_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = internal.n_G_D_progress[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_G_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = internal.aux_G_D[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_G_D_progress[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_conf_to_ICU_W_D_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_D_conf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_D_conf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_unconf_to_ICU_W_D_unconf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_D_unconf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_D_unconf_progress[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_EE_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_E_progress[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_E_next_vacc_class[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_E_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_E_progress[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_E_next_vacc_class[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_D_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_conf_to_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_ICU_pre_conf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_conf_to_ICU_D_conf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_conf_to_ICU_W_D_conf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_ICU_pre_unconf_to_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_ICU_pre_unconf_progress[shared->dim_ICU_pre_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_unconf_to_ICU_D_unconf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_unconf_to_ICU_W_D_unconf[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_II_A_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_A_progress[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_I_A_next_vacc_class[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_II_P_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_P_progress[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_I_P_next_vacc_class[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_A_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_A_progress[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_I_A_next_vacc_class[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_hosp[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_I_C_2_progress[shared->dim_I_C_2_123 * (k - 1) + shared->dim_E_12 * (shared->k_C_2 - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_G_D[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_I_P_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_P_progress[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_I_P_next_vacc_class[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_RS_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_R_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_next_vacc_class_tmp[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_R_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_sero_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_sero_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = T_sero_pos[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_T_sero_pos_progress[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_T_sero_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = internal.new_T_sero_pos[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_T_sero_pre_to_T_sero_pos[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_sero_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_sero_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = internal.new_T_sero_pos[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_T_sero_pos_progress[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_D_unconf_to_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.aux_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_D_unconf_to_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_W_R_unconf_to_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.aux_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_W_R_unconf_to_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.p_SE[i - 1 + 19 * (j - 1)] = 1 - std::exp(- odin_sum2<real_t>(internal.lambda.data(), i - 1, i, 0, shared->dim_strain_transmission, 19) * shared->rel_susceptibility[shared->dim_rel_susceptibility_1 * (j - 1) + i - 1] * shared->dt);
      }
    }
    state_next[16] = (fmodr<real_t>(step, shared->steps_per_day) == 0 ? delta_D_carehomes_tot : D_carehomes_inc + delta_D_carehomes_tot);
    state_next[15] = D_carehomes_tot + delta_D_carehomes_tot;
    state_next[14] = (fmodr<real_t>(step, shared->steps_per_day) == 0 ? delta_D_comm_tot : D_comm_inc + delta_D_comm_tot);
    state_next[13] = D_comm_tot + delta_D_comm_tot;
    state_next[17] = (fmodr<real_t>(step, shared->steps_per_day) == 0 ? delta_D_hosp_tot : D_hosp_inc + delta_D_hosp_tot);
    state_next[12] = D_hosp_tot + delta_D_hosp_tot;
    state_next[18] = D_tot + delta_D_hosp_tot + delta_D_comm_tot + delta_D_carehomes_tot;
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_I_C_2 + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = internal.new_I_C_2[shared->dim_I_C_2_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        state_next[shared->offset_variable_prob_strain + i - 1 + 19 * (j - 1)] = internal.lambda[19 * (j - 1) + i - 1] / (real_t) odin_sum2<real_t>(internal.lambda.data(), i - 1, i, 0, shared->dim_strain_transmission, 19);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_conf_to_ICU_W_R_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_R_conf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_R_conf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_unconf_to_ICU_W_R_unconf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_R_unconf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_R_unconf_progress[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_EE[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.n_E_progress[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_EE_next_vacc_class[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_EI_A_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_EE_next_vacc_class[shared->dim_E_123 * (k - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1]), 1 - shared->p_C[i - 1] * shared->rel_p_sympt[19 * (k - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_D_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_II_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.n_I_A_progress[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_II_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.n_I_P_progress[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_ICU_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_to_hosp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_ICU_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_RS[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_RS_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_next_vacc_class_capped[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = std::min(internal.n_R_next_vacc_class_tmp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], std::min(T_sero_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], T_PCR_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]));
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.n_S_progress_tot[i - 1 + 19 * (j - 1)] = dust::distr::rbinom(rng_state, std::round(S[19 * (j - 1) + i - 1]), internal.p_SE[19 * (j - 1) + i - 1]);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_G_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = G_D[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_G_D[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_D_unconf_to_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.aux_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_D_unconf_to_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_sero_pos; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_T_sero_pos + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = internal.new_T_sero_pos[shared->dim_T_sero_pos_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_W_D_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.new_W_D_conf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_W_D_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = internal.new_W_D_unconf[shared->dim_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_W_R_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.new_W_R_conf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_W_R_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = internal.new_W_R_unconf[shared->dim_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      state_next[65 + i - 1] = cum_admit_by_age[i - 1] + odin_sum3<real_t>(internal.n_I_C_2_to_hosp.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12);
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_E_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_E_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(internal.n_E_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_E, j - 1, j, 19, shared->dim_E_12, shared->dim_E_123) + odin_sum4<real_t>(internal.n_EE_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_E, j - 1, j, 19, shared->dim_E_12, shared->dim_E_123);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_I_A_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_I_A_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(internal.n_I_A_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_A, j - 1, j, 19, shared->dim_E_12, shared->dim_I_A_123) + odin_sum4<real_t>(internal.n_II_A_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_A, j - 1, j, 19, shared->dim_E_12, shared->dim_I_A_123);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_I_P_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_I_P_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(internal.n_I_P_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_P, j - 1, j, 19, shared->dim_E_12, shared->dim_I_P_123) + odin_sum4<real_t>(internal.n_II_P_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, 0, shared->k_P, j - 1, j, 19, shared->dim_E_12, shared->dim_I_P_123);
      }
    }
    state_next[19] = odin_sum4<real_t>(internal.new_T_sero_pos.data(), 3, 13, 0, shared->dim_strain_transmission, 0, shared->k_sero_pos, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12, shared->dim_T_sero_pos_123);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        int l = 1;
        internal.aux_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = internal.n_II_P[shared->dim_I_P_123 * 0 + shared->dim_E_12 * (shared->k_P - 1) + 19 * (j - 1) + i - 1] + internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (shared->k_P - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.aux_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = internal.n_II_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (shared->k_P - 1) + 19 * (j - 1) + i - 1] + internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (l - 1 - 1) + shared->dim_E_12 * (shared->k_P - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_C_1; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = internal.n_I_C_1_progress[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_1; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = internal.aux_I_C_1[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_1_progress[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_EI_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_EE[shared->dim_E_123 * (k - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1]), 1 - shared->p_C[i - 1] * shared->rel_p_sympt[19 * (k - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_EI_P_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_EE_next_vacc_class[shared->dim_E_123 * (k - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] - internal.n_EI_A_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.n_ICU_W_R_unconf_to_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.aux_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), shared->p_test);
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_to_ICU_pre[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_star_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_R_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = (shared->model_pcr_and_serology == 1 ? internal.n_R_next_vacc_class_capped[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] : internal.n_R_next_vacc_class_tmp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      int j = 1;
      for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
        internal.n_S_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_S_progress_tot[19 * (k - 1) + i - 1]), internal.lambda[19 * 0 + i - 1] / (real_t) odin_sum2<real_t>(internal.lambda.data(), i - 1, i, 0, shared->dim_strain_transmission, 19));
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 2; j <= shared->n_strains; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_S_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_S_progress_tot[19 * (k - 1) + i - 1] - odin_sum3<real_t>(internal.n_S_progress.data(), i - 1, i, 0, j - 1, k - 1, k, 19, shared->dim_E_12)), internal.lambda[19 * (j - 1) + i - 1] / (real_t) odin_sum2<real_t>(internal.lambda.data(), i - 1, i, j - 1, shared->n_strains, 19));
        }
      }
    }
    {
       int i = 4;
       for (int j = 2; j <= shared->n_strains; ++j) {
         int k = 1;
         internal.n_S_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = std::min(internal.n_S_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + strain_seed, internal.n_S_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + S[19 * (k - 1) + i - 1] - odin_sum3<real_t>(internal.n_S_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, k - 1, k, 19, shared->dim_E_12));
       }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        int l = 1;
        internal.n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_EE[shared->dim_E_123 * 0 + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1]), shared->p_sero_pre_1);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_EE[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (l - 1 - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1]), shared->p_sero_pre_1);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 2;
        int l = 1;
        internal.n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = internal.n_EE[shared->dim_E_123 * 0 + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] - internal.n_com_to_T_sero_pre[shared->dim_T_sero_pre_123 * 0 + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 2;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = internal.n_EE[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (l - 1 - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] - internal.n_com_to_T_sero_pre[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_hosp_non_ICU[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_I_C_2_to_hosp[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_ICU_pre[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_D_unconf_to_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.aux_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_D_unconf_to_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_G_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_G_D + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = internal.new_G_D[shared->dim_G_D_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_D_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.new_ICU_D_conf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_D_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = internal.new_ICU_D_unconf[shared->dim_ICU_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.n_EI_A[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.n_II_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_A_next_vacc_class[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          int l = 1;
          internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_I_A_next_vacc_class[shared->dim_I_A_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_I_A_next_vacc_class[shared->dim_I_A_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        int l = 1;
        internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EI_A_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EI_A_next_vacc_class[shared->dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_A; ++k) {
          int l = 1;
          internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_A; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    real_t delta_new_conf = odin_sum1<real_t>(internal.n_H_D_unconf_to_conf.data(), 0, shared->dim_H_D_unconf) + odin_sum1<real_t>(internal.n_H_R_unconf_to_conf.data(), 0, shared->dim_H_R_unconf) + odin_sum1<real_t>(internal.n_ICU_pre_unconf_to_conf.data(), 0, shared->dim_ICU_pre_unconf) + odin_sum1<real_t>(internal.n_ICU_D_unconf_to_conf.data(), 0, shared->dim_ICU_D_unconf) + odin_sum1<real_t>(internal.n_ICU_W_R_unconf_to_conf.data(), 0, shared->dim_ICU_W_R_unconf) + odin_sum1<real_t>(internal.n_ICU_W_D_unconf_to_conf.data(), 0, shared->dim_ICU_W_D_unconf) + odin_sum1<real_t>(internal.n_W_R_unconf_to_conf.data(), 0, shared->dim_W_R_unconf) + odin_sum1<real_t>(internal.n_W_D_unconf_to_conf.data(), 0, shared->dim_W_D_unconf);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_EI_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_EE[shared->dim_E_123 * (k - 1) + shared->dim_E_12 * (shared->k_E - 1) + 19 * (j - 1) + i - 1] - internal.n_EI_A[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_H_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_hosp_non_ICU[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_H_D_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_SE_next_vacc_class[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_S_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_S_next_vacc_class[19 * (k - 1) + i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.n_S_next_vacc_class[i - 1 + 19 * (j - 1)] = dust::distr::rbinom(rng_state, std::round(S[19 * (j - 1) + i - 1] - odin_sum3<real_t>(internal.n_S_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12)), internal.p_S_next_vacc_class[19 * (j - 1) + i - 1]);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_W_R_unconf_to_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.aux_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_W_R_unconf_to_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_ICU_pre_unconf_to_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.new_ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_ICU_pre_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.aux_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_ICU_pre_unconf_to_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.new_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_ICU_pre[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_ICU_pre_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_1; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = I_C_1[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_I_C_1[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.new_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_A[shared->dim_I_A_123 * (k - 1) + shared->dim_E_12 * (shared->k_A - 1) + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_R_conf_progress[shared->dim_H_R_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_H_R - 1) + 19 * (j - 1) + i - 1] + internal.n_H_R_unconf_progress[shared->dim_H_R_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_H_R - 1) + 19 * (j - 1) + i - 1] + internal.n_W_R_conf_progress[shared->dim_W_R_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_W_R - 1) + 19 * (j - 1) + i - 1] + internal.n_W_R_unconf_progress[shared->dim_W_R_unconf_123 * (k - 1) + shared->dim_E_12 * (shared->k_W_R - 1) + 19 * (j - 1) + i - 1] - internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        internal.new_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_R[shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (shared->k_A - 1) + 19 * (j - 1) + i - 1] + internal.n_R_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->n_vacc_classes; ++k) {
          internal.new_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_A_next_vacc_class[shared->dim_I_A_123 * (k - 1 - 1) + shared->dim_E_12 * (shared->k_A - 1) + 19 * (j - 1) + i - 1] + internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.new_T_PCR_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = T_PCR_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_T_PCR_pos_progress[shared->dim_T_PCR_pos_123 * (k - 1) + shared->dim_E_12 * (shared->k_PCR_pos - 1) + 19 * (j - 1) + i - 1] - shared->model_pcr_and_serology * internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        internal.new_T_PCR_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_PCR_neg[shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->n_vacc_classes; ++k) {
          internal.new_T_PCR_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_PCR_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_PCR_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = T_PCR_pre[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_T_PCR_pre_progress[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_T_PCR_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = internal.new_T_PCR_pre[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_S_progress[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_PCR_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_PCR_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = internal.new_T_PCR_pre[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_T_PCR_pre_progress[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.new_T_sero_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = T_sero_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + odin_sum4<real_t>(internal.n_T_sero_pre_progress.data(), i - 1, i, j - 1, j, 0, 2, k - 1, k, 19, shared->dim_E_12, shared->dim_T_sero_pre_123) - internal.n_T_sero_pre_to_T_sero_pos[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_T_sero_pos_progress[shared->dim_T_sero_pos_123 * (k - 1) + shared->dim_E_12 * (shared->k_sero_pos - 1) + 19 * (j - 1) + i - 1] - shared->model_pcr_and_serology * internal.n_R_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        internal.new_T_sero_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_sero_neg[shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->n_vacc_classes; ++k) {
          internal.new_T_sero_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_sero_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + shared->model_pcr_and_serology * internal.n_R_next_vacc_class[shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= 2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = T_sero_pre[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_com_to_T_sero_pre[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_T_sero_pre_progress[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_W_D_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.new_ICU_W_D_conf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_W_D_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = internal.new_ICU_W_D_unconf[shared->dim_ICU_W_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    state_next[3] = cum_infections + odin_sum1<real_t>(internal.n_S_progress.data(), 0, shared->dim_R);
    for (int i = 1; i <= shared->dim_strain_transmission; ++i) {
      state_next[141 + i - 1] = cum_infections_per_strain[i - 1] + odin_sum3<real_t>(internal.n_S_progress.data(), 0, 19, i - 1, i, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12);
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_R_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_R_vaccinated[19 * (j - 1) + i - 1] + odin_sum3<real_t>(internal.n_R_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12) + odin_sum3<real_t>(internal.n_RS_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12);
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.n_EI_P[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.n_II_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_P_next_vacc_class[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          int l = 1;
          internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_I_P_next_vacc_class[shared->dim_I_P_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_I_P_next_vacc_class[shared->dim_I_P_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        int l = 1;
        internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EI_P_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EI_P_next_vacc_class[shared->dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_P; ++k) {
          int l = 1;
          internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_P; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_II_P_next_vacc_class[shared->dim_I_P_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_to_H_D[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_star_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_H_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_hosp_non_ICU[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_H_D[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->n_vacc_classes; ++k) {
          internal.n_SE[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.n_S_progress[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_SE_next_vacc_class[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    real_t new_ICU_tot = odin_sum1<real_t>(internal.new_ICU_W_R_conf.data(), 0, shared->dim_ICU_W_R_unconf) + odin_sum1<real_t>(internal.new_ICU_W_D_conf.data(), 0, shared->dim_ICU_W_D_unconf) + odin_sum1<real_t>(internal.new_ICU_D_conf.data(), 0, shared->dim_ICU_D_unconf);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        internal.new_S[i - 1 + 19 * (j - 1)] = S[19 * (j - 1) + i - 1] + odin_sum3<real_t>(internal.n_RS.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12) - odin_sum3<real_t>(internal.n_S_progress.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12) - internal.n_S_next_vacc_class[19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      int j = 1;
      internal.new_S[i - 1 + 19 * (j - 1)] = internal.new_S[19 * 0 + i - 1] + internal.n_S_next_vacc_class[19 * (shared->n_vacc_classes - 1) + i - 1] + odin_sum3<real_t>(internal.n_RS_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, shared->n_vacc_classes - 1, shared->n_vacc_classes, 19, shared->dim_E_12);
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 2; j <= shared->n_vacc_classes; ++j) {
        internal.new_S[i - 1 + 19 * (j - 1)] = internal.new_S[19 * (j - 1) + i - 1] + internal.n_S_next_vacc_class[19 * (j - 1 - 1) + i - 1] + odin_sum3<real_t>(internal.n_RS_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1 - 1, j - 1, 19, shared->dim_E_12);
      }
    }
    real_t new_sympt_cases = odin_sum1<real_t>(internal.n_EI_P.data(), 0, shared->dim_R) + odin_sum1<real_t>(internal.n_EI_P_next_vacc_class.data(), 0, shared->dim_R);
    real_t new_sympt_cases_non_variant_over25 = odin_sum3<real_t>(internal.n_EI_P.data(), 5, shared->n_groups, 0, 1, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12) + odin_sum3<real_t>(internal.n_EI_P_next_vacc_class.data(), 5, shared->n_groups, 0, 1, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12);
    real_t new_sympt_cases_over25 = odin_sum3<real_t>(internal.n_EI_P.data(), 5, shared->n_groups, 0, shared->dim_strain_transmission, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12) + odin_sum3<real_t>(internal.n_EI_P_next_vacc_class.data(), 5, shared->n_groups, 0, shared->dim_strain_transmission, 0, shared->dim_rel_susceptibility_2, 19, shared->dim_E_12);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_W_R_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.new_ICU_W_R_conf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_W_R_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = internal.new_ICU_W_R_unconf[shared->dim_ICU_W_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_pre_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.new_ICU_pre_conf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_ICU_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_ICU_pre_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = internal.new_ICU_pre_unconf[shared->dim_ICU_pre_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_C_1; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_I_C_1 + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = internal.new_I_C_1[shared->dim_I_C_1_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          state_next[shared->offset_variable_R + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          state_next[shared->offset_variable_T_PCR_neg + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_PCR_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_PCR_pre; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_T_PCR_pre + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = internal.new_T_PCR_pre[shared->dim_T_PCR_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          state_next[shared->offset_variable_T_sero_neg + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = internal.new_T_sero_neg[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= 2; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_T_sero_pre + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = internal.new_T_sero_pre[shared->dim_T_sero_pre_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_cum_n_S_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_S_vaccinated[19 * (j - 1) + i - 1] + internal.n_S_next_vacc_class[19 * (j - 1) + i - 1] + odin_sum3<real_t>(internal.n_SE_next_vacc_class.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12);
      }
    }
    state_next[5] = cum_new_conf + delta_new_conf;
    state_next[2] = (fmodr<real_t>(step, shared->steps_per_day) == 0 ? delta_new_conf : new_conf_inc + delta_new_conf);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.n_SE[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.n_EE[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_EE[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_EE_next_vacc_class[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_E_next_vacc_class[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          int l = 1;
          internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * 0 + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_E_next_vacc_class[shared->dim_E_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_E_next_vacc_class[shared->dim_E_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        int l = 1;
        internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * 0 + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_SE_next_vacc_class[shared->dim_E_12 * (shared->n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 2; l <= shared->n_vacc_classes; ++l) {
          internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_SE_next_vacc_class[shared->dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_E; ++k) {
          int l = 1;
          internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (shared->n_vacc_classes - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 2; k <= shared->k_E; ++k) {
          for (int l = 2; l <= shared->n_vacc_classes; ++l) {
            internal.aux_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_EE_next_vacc_class[shared->dim_E_123 * (l - 1 - 1) + shared->dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.n_I_C_2_to_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(internal.n_I_C_2_to_H_R[shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), internal.p_star_by_age[i - 1]);
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_D_unconf_to_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.new_H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_H_D_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.aux_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_D_unconf_to_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.new_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_H_D[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_H_D_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    state_next[9] = new_ICU_tot;
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_A; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_I_A + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = internal.new_I_A[shared->dim_I_A_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_S + i - 1 + 19 * (j - 1)] = internal.new_S[19 * (j - 1) + i - 1];
      }
    }
    state_next[20] = cum_sympt_cases + new_sympt_cases;
    state_next[22] = cum_sympt_cases_non_variant_over25 + new_sympt_cases_non_variant_over25;
    state_next[21] = cum_sympt_cases_over25 + new_sympt_cases_over25;
    state_next[23] = ((fmodr<real_t>(step, shared->steps_per_day) == 0 ? new_sympt_cases : sympt_cases_inc + new_sympt_cases));
    state_next[25] = ((fmodr<real_t>(step, shared->steps_per_day) == 0 ? new_sympt_cases_non_variant_over25 : sympt_cases_non_variant_over25_inc + new_sympt_cases_non_variant_over25));
    state_next[24] = ((fmodr<real_t>(step, shared->steps_per_day) == 0 ? new_sympt_cases_over25 : sympt_cases_over25_inc + new_sympt_cases_over25));
    real_t delta_admit_conf = odin_sum1<real_t>(internal.n_I_C_2_to_H_D_conf.data(), 0, shared->dim_R) + odin_sum1<real_t>(internal.n_I_C_2_to_H_R_conf.data(), 0, shared->dim_R) + odin_sum1<real_t>(internal.n_I_C_2_to_ICU_pre_conf.data(), 0, shared->dim_R);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.aux_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + internal.n_H_R_unconf_to_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.new_H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_H_R_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            internal.new_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.aux_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - internal.n_H_R_unconf_to_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        int k = 1;
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          internal.new_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.new_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * 0 + 19 * (j - 1) + i - 1] + internal.n_I_C_2_to_H_R[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - internal.n_I_C_2_to_H_R_conf[shared->dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_H_D_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.new_H_D_conf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_D; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_H_D_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = internal.new_H_D_unconf[shared->dim_H_D_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_P; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_I_P + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = internal.new_I_P[shared->dim_I_P_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
          internal.I_weighted_strain[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = shared->strain_transmission[j - 1] * (shared->I_A_transmission * odin_sum4<real_t>(internal.new_I_A.data(), i - 1, i, j - 1, j, 0, shared->k_A, k - 1, k, 19, shared->dim_E_12, shared->dim_I_A_123) + shared->I_P_transmission * odin_sum4<real_t>(internal.new_I_P.data(), i - 1, i, j - 1, j, 0, shared->k_P, k - 1, k, 19, shared->dim_E_12, shared->dim_I_P_123) + shared->I_C_1_transmission * odin_sum4<real_t>(internal.new_I_C_1.data(), i - 1, i, j - 1, j, 0, shared->k_C_1, k - 1, k, 19, shared->dim_E_12, shared->dim_I_C_1_123) + shared->I_C_2_transmission * odin_sum4<real_t>(internal.new_I_C_2.data(), i - 1, i, j - 1, j, 0, shared->k_C_2, k - 1, k, 19, shared->dim_E_12, shared->dim_I_C_2_123) + shared->hosp_transmission * (odin_sum4<real_t>(internal.new_ICU_pre_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_pre, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(internal.new_ICU_pre_conf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_pre, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_pre_unconf_123) + odin_sum4<real_t>(internal.new_H_R_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_H_R, k - 1, k, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(internal.new_H_R_conf.data(), i - 1, i, j - 1, j, 0, shared->k_H_R, k - 1, k, 19, shared->dim_E_12, shared->dim_H_R_unconf_123) + odin_sum4<real_t>(internal.new_H_D_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_H_D, k - 1, k, 19, shared->dim_E_12, shared->dim_H_D_unconf_123) + odin_sum4<real_t>(internal.new_H_D_conf.data(), i - 1, i, j - 1, j, 0, shared->k_H_D, k - 1, k, 19, shared->dim_E_12, shared->dim_H_D_unconf_123)) + shared->ICU_transmission * (odin_sum4<real_t>(internal.new_ICU_W_R_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_W_R, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(internal.new_ICU_W_R_conf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_W_R, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(internal.new_ICU_W_D_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_W_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(internal.new_ICU_W_D_conf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_W_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(internal.new_ICU_D_unconf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123) + odin_sum4<real_t>(internal.new_ICU_D_conf.data(), i - 1, i, j - 1, j, 0, shared->k_ICU_D, k - 1, k, 19, shared->dim_E_12, shared->dim_ICU_D_unconf_123)) + shared->G_D_transmission * odin_sum4<real_t>(internal.new_G_D.data(), i - 1, i, j - 1, j, 0, shared->k_G_D, k - 1, k, 19, shared->dim_E_12, shared->dim_G_D_123));
        }
      }
    }
    real_t new_general_tot = odin_sum1<real_t>(internal.new_ICU_pre_conf.data(), 0, shared->dim_ICU_pre_unconf) + odin_sum1<real_t>(internal.new_H_R_conf.data(), 0, shared->dim_H_R_unconf) + odin_sum1<real_t>(internal.new_H_D_conf.data(), 0, shared->dim_H_D_unconf) + odin_sum1<real_t>(internal.new_W_R_conf.data(), 0, shared->dim_W_R_unconf) + odin_sum1<real_t>(internal.new_W_D_conf.data(), 0, shared->dim_W_D_unconf);
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_E; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_E + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = internal.new_E[shared->dim_E_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_H_R_conf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.new_H_R_conf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
        for (int k = 1; k <= shared->k_H_R; ++k) {
          for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
            state_next[shared->offset_variable_H_R_unconf + i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = internal.new_H_R_unconf[shared->dim_H_R_unconf_123 * (l - 1) + shared->dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
          }
        }
      }
    }
    state_next[1] = (fmodr<real_t>(step, shared->steps_per_day) == 0 ? delta_admit_conf : admit_conf_inc + delta_admit_conf);
    state_next[4] = cum_admit_conf + delta_admit_conf;
    for (int i = 1; i <= 19; ++i) {
      for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
        state_next[shared->offset_variable_I_weighted + i - 1 + 19 * (j - 1)] = odin_sum3<real_t>(internal.I_weighted_strain.data(), i - 1, i, 0, shared->dim_strain_transmission, j - 1, j, 19, shared->dim_E_12);
      }
    }
    state_next[10] = new_general_tot;
    state_next[11] = new_ICU_tot + new_general_tot;
  }
  real_t compare_data(const real_t * state, const data_t& data, dust::rng_state_t<real_t>& rng_state) {
    return compare<carehomes>(state, data, internal, shared, rng_state);
  }
private:
  std::shared_ptr<const shared_t> shared;
  internal_t internal;
};
namespace dust {
template <>
struct has_gpu_support<carehomes> : std::true_type {};
}
template <>
size_t dust::device_shared_size_int<carehomes>(dust::shared_ptr<carehomes> shared) {
  return 119 + 2;
}
template <>
size_t dust::device_shared_size_real<carehomes>(dust::shared_ptr<carehomes> shared) {
  return 28 + shared->dim_beta_step + shared->dim_p_G_D_step + shared->dim_p_H_step + shared->dim_p_H_D_step + shared->dim_p_ICU_step + shared->dim_p_ICU_D_step + shared->dim_p_W_D_step + shared->dim_p_star_step + shared->dim_strain_seed_step + 19 + 19 + 19 + 19 + 19 + 19 + 19 + shared->dim_strain_transmission + 19 + 19 + 19 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 361 + shared->dim_cum_n_S_vaccinated + shared->dim_rel_susceptibility + shared->dim_cum_n_S_vaccinated + shared->dim_vaccine_dose_step + shared->dim_T_sero_pre;
}
template <>
size_t dust::device_internal_size_int<carehomes>(dust::shared_ptr<carehomes> shared) {
  return 0;
}
template <>
size_t dust::device_internal_size_real<carehomes>(dust::shared_ptr<carehomes> shared) {
  return 0 + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_G_D + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_I_A + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_I_P + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + 19 + 19 + shared->dim_E_12 + shared->dim_E + shared->dim_E + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_E + shared->dim_G_D + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_pre_unconf + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_ICU_pre_unconf + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_ICU_pre_unconf + shared->dim_I_A + shared->dim_I_A + shared->dim_I_P + shared->dim_I_P + shared->dim_I_A + shared->dim_I_A + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_I_P + shared->dim_I_P + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_cum_n_S_vaccinated + shared->dim_T_PCR_pos + shared->dim_T_PCR_pre + shared->dim_T_sero_pos + shared->dim_T_sero_pre + shared->dim_R + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_T_sero_pre + shared->dim_R + shared->dim_E + shared->dim_G_D + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_I_A + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_I_P + shared->dim_R + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_T_PCR_pos + shared->dim_T_PCR_pre + shared->dim_R + shared->dim_T_sero_pos + shared->dim_T_sero_pre + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_E + 19 + 19 + 19 + 19 + 19 + shared->dim_I_A + shared->dim_I_P + shared->dim_R + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 19 + 19 + shared->dim_s_ij + 38 + shared->dim_cum_n_S_vaccinated + 38;
}
template <>
void dust::device_shared_copy<carehomes>(dust::shared_ptr<carehomes> shared, int * shared_int, carehomes::real_t * shared_real) {
  shared_int = dust::shared_copy(shared_int, shared->dim_beta_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_G_D_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_H_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_H_D_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_ICU_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_ICU_D_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_W_D_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_p_star_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_strain_seed_step);
  shared_int = dust::shared_copy(shared_int, shared->dim_cum_n_S_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->dim_R);
  shared_int = dust::shared_copy(shared_int, shared->dim_E);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_A);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_P);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_C_1);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_C_2);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_pre_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_H_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_H_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_W_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_W_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_W_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_W_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->dim_G_D);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_PCR_pos);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_PCR_pre);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_sero_pos);
  shared_int = dust::shared_copy(shared_int, shared->dim_vaccine_dose_step);
  shared_int = dust::shared_copy(shared_int, shared->k_H_D);
  shared_int = dust::shared_copy(shared_int, shared->k_H_R);
  shared_int = dust::shared_copy(shared_int, shared->k_ICU_pre);
  shared_int = dust::shared_copy(shared_int, shared->k_C_1);
  shared_int = dust::shared_copy(shared_int, shared->k_C_2);
  shared_int = dust::shared_copy(shared_int, shared->k_ICU_W_D);
  shared_int = dust::shared_copy(shared_int, shared->k_W_D);
  shared_int = dust::shared_copy(shared_int, shared->k_ICU_W_R);
  shared_int = dust::shared_copy(shared_int, shared->k_W_R);
  shared_int = dust::shared_copy(shared_int, shared->k_ICU_D);
  shared_int = dust::shared_copy(shared_int, shared->k_G_D);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_sero_pre);
  shared_int = dust::shared_copy(shared_int, shared->k_PCR_pre);
  shared_int = dust::shared_copy(shared_int, shared->k_PCR_pos);
  shared_int = dust::shared_copy(shared_int, shared->dim_s_ij);
  shared_int = dust::shared_copy(shared_int, shared->n_age_groups);
  shared_int = dust::shared_copy(shared_int, shared->n_groups);
  shared_int = dust::shared_copy(shared_int, shared->dim_E_12);
  shared_int = dust::shared_copy(shared_int, shared->k_sero_pos);
  shared_int = dust::shared_copy(shared_int, shared->steps_per_day);
  shared_int = dust::shared_copy(shared_int, shared->k_E);
  shared_int = dust::shared_copy(shared_int, shared->k_P);
  shared_int = dust::shared_copy(shared_int, shared->n_vacc_classes);
  shared_int = dust::shared_copy(shared_int, shared->n_strains);
  shared_int = dust::shared_copy(shared_int, shared->k_A);
  shared_int = dust::shared_copy(shared_int, shared->dim_strain_transmission);
  shared_int = dust::shared_copy(shared_int, shared->dim_E_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_G_D_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_H_D_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_H_R_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_D_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_W_D_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_W_R_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_ICU_pre_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_A_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_C_1_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_C_2_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_I_P_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_PCR_pos_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_PCR_pre_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_sero_pos_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_T_sero_pre_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_W_D_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_W_R_unconf_123);
  shared_int = dust::shared_copy(shared_int, shared->dim_rel_susceptibility);
  shared_int = dust::shared_copy(shared_int, shared->dim_rel_susceptibility_1);
  shared_int = dust::shared_copy(shared_int, shared->dim_rel_susceptibility_2);
  shared_int = dust::shared_copy(shared_int, shared->dim_vaccine_dose_step_1);
  shared_int = dust::shared_copy(shared_int, shared->dim_vaccine_dose_step_12);
  shared_int = dust::shared_copy(shared_int, shared->dim_vaccine_dose_step_2);
  shared_int = dust::shared_copy(shared_int, shared->dim_vaccine_dose_step_3);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_E);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_G_D);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_H_D_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_H_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_H_R_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_H_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_D_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_W_D_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_W_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_W_R_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_W_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_pre_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_ICU_pre_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_I_A);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_I_C_1);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_I_C_2);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_I_P);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_I_weighted);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_R);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_S);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_PCR_neg);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_PCR_pos);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_PCR_pre);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_sero_neg);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_sero_pos);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_T_sero_pre);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_W_D_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_W_D_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_W_R_conf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_W_R_unconf);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_E_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_I_A_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_I_P_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_R_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_S_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_cum_n_vaccinated);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_prob_strain);
  shared_int = dust::shared_copy(shared_int, shared->offset_variable_tmp_vaccine_probability);
  shared_int = dust::shared_copy(shared_int, shared->index_dose);
  shared_real = dust::shared_copy(shared_real, shared->dt);
  shared_real = dust::shared_copy(shared_real, shared->I_A_transmission);
  shared_real = dust::shared_copy(shared_real, shared->I_P_transmission);
  shared_real = dust::shared_copy(shared_real, shared->I_C_1_transmission);
  shared_real = dust::shared_copy(shared_real, shared->I_C_2_transmission);
  shared_real = dust::shared_copy(shared_real, shared->hosp_transmission);
  shared_real = dust::shared_copy(shared_real, shared->ICU_transmission);
  shared_real = dust::shared_copy(shared_real, shared->G_D_transmission);
  shared_real = dust::shared_copy(shared_real, shared->p_E_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_G_D_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_H_D_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_H_R_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_D_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_W_D_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_W_R_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_pre_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_I_A_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_I_C_1_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_I_C_2_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_I_P_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_T_PCR_pos_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_T_PCR_pre_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_T_sero_pos_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_W_D_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_W_R_progress);
  shared_real = dust::shared_copy(shared_real, shared->p_test);
  shared_real = dust::shared_copy(shared_real, shared->model_pcr_and_serology);
  shared_real = dust::shared_copy(shared_real, shared->p_sero_pre_1);
  shared_real = dust::shared_copy(shared_real, shared->beta_step);
  shared_real = dust::shared_copy(shared_real, shared->p_G_D_step);
  shared_real = dust::shared_copy(shared_real, shared->p_H_step);
  shared_real = dust::shared_copy(shared_real, shared->p_H_D_step);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_step);
  shared_real = dust::shared_copy(shared_real, shared->p_ICU_D_step);
  shared_real = dust::shared_copy(shared_real, shared->p_W_D_step);
  shared_real = dust::shared_copy(shared_real, shared->p_star_step);
  shared_real = dust::shared_copy(shared_real, shared->strain_seed_step);
  shared_real = dust::shared_copy(shared_real, shared->psi_G_D);
  shared_real = dust::shared_copy(shared_real, shared->psi_H_D);
  shared_real = dust::shared_copy(shared_real, shared->psi_H);
  shared_real = dust::shared_copy(shared_real, shared->psi_ICU_D);
  shared_real = dust::shared_copy(shared_real, shared->psi_ICU);
  shared_real = dust::shared_copy(shared_real, shared->psi_W_D);
  shared_real = dust::shared_copy(shared_real, shared->psi_star);
  shared_real = dust::shared_copy(shared_real, shared->strain_transmission);
  shared_real = dust::shared_copy(shared_real, shared->p_RS);
  shared_real = dust::shared_copy(shared_real, shared->p_sero_pos);
  shared_real = dust::shared_copy(shared_real, shared->p_C);
  shared_real = dust::shared_copy(shared_real, shared->rel_infectivity);
  shared_real = dust::shared_copy(shared_real, shared->rel_p_hosp_if_sympt);
  shared_real = dust::shared_copy(shared_real, shared->m);
  shared_real = dust::shared_copy(shared_real, shared->vaccine_progression_rate_base);
  shared_real = dust::shared_copy(shared_real, shared->rel_susceptibility);
  shared_real = dust::shared_copy(shared_real, shared->rel_p_sympt);
  shared_real = dust::shared_copy(shared_real, shared->vaccine_dose_step);
  shared_real = dust::shared_copy(shared_real, shared->p_T_sero_pre_progress);
}
template<>
DEVICE void update_device<carehomes>(size_t step, const dust::interleaved<carehomes::real_t> state, dust::interleaved<int> internal_int, dust::interleaved<carehomes::real_t> internal_real, const int * shared_int, const carehomes::real_t * shared_real, dust::rng_state_t<carehomes::real_t>& rng_state, dust::interleaved<carehomes::real_t> state_next) {
  typedef carehomes::real_t real_t;
  int dim_beta_step = shared_int[0];
  int dim_p_G_D_step = shared_int[1];
  int dim_p_H_step = shared_int[2];
  int dim_p_H_D_step = shared_int[3];
  int dim_p_ICU_step = shared_int[4];
  int dim_p_ICU_D_step = shared_int[5];
  int dim_p_W_D_step = shared_int[6];
  int dim_p_star_step = shared_int[7];
  int dim_strain_seed_step = shared_int[8];
  int dim_cum_n_S_vaccinated = shared_int[9];
  int dim_R = shared_int[10];
  int dim_E = shared_int[11];
  int dim_I_A = shared_int[12];
  int dim_I_P = shared_int[13];
  int dim_I_C_1 = shared_int[14];
  int dim_I_C_2 = shared_int[15];
  int dim_ICU_pre_unconf = shared_int[16];
  int dim_H_R_unconf = shared_int[17];
  int dim_H_D_unconf = shared_int[18];
  int dim_ICU_W_R_unconf = shared_int[19];
  int dim_ICU_W_D_unconf = shared_int[20];
  int dim_ICU_D_unconf = shared_int[21];
  int dim_W_R_unconf = shared_int[22];
  int dim_W_D_unconf = shared_int[23];
  int dim_G_D = shared_int[24];
  int dim_T_PCR_pos = shared_int[25];
  int dim_T_PCR_pre = shared_int[26];
  int dim_T_sero_pos = shared_int[27];
  int dim_vaccine_dose_step = shared_int[28];
  int k_H_D = shared_int[29];
  int k_H_R = shared_int[30];
  int k_ICU_pre = shared_int[31];
  int k_C_1 = shared_int[32];
  int k_C_2 = shared_int[33];
  int k_ICU_W_D = shared_int[34];
  int k_W_D = shared_int[35];
  int k_ICU_W_R = shared_int[36];
  int k_W_R = shared_int[37];
  int k_ICU_D = shared_int[38];
  int k_G_D = shared_int[39];
  int dim_T_sero_pre = shared_int[40];
  int k_PCR_pre = shared_int[41];
  int k_PCR_pos = shared_int[42];
  int dim_s_ij = shared_int[43];
  int n_age_groups = shared_int[44];
  int n_groups = shared_int[45];
  int dim_E_12 = shared_int[46];
  int k_sero_pos = shared_int[47];
  int steps_per_day = shared_int[48];
  int k_E = shared_int[49];
  int k_P = shared_int[50];
  int n_vacc_classes = shared_int[51];
  int n_strains = shared_int[52];
  int k_A = shared_int[53];
  int dim_strain_transmission = shared_int[54];
  int dim_E_123 = shared_int[55];
  int dim_G_D_123 = shared_int[56];
  int dim_H_D_unconf_123 = shared_int[57];
  int dim_H_R_unconf_123 = shared_int[58];
  int dim_ICU_D_unconf_123 = shared_int[59];
  int dim_ICU_W_D_unconf_123 = shared_int[60];
  int dim_ICU_W_R_unconf_123 = shared_int[61];
  int dim_ICU_pre_unconf_123 = shared_int[62];
  int dim_I_A_123 = shared_int[63];
  int dim_I_C_1_123 = shared_int[64];
  int dim_I_C_2_123 = shared_int[65];
  int dim_I_P_123 = shared_int[66];
  int dim_T_PCR_pos_123 = shared_int[67];
  int dim_T_PCR_pre_123 = shared_int[68];
  int dim_T_sero_pos_123 = shared_int[69];
  int dim_T_sero_pre_123 = shared_int[70];
  int dim_W_D_unconf_123 = shared_int[71];
  int dim_W_R_unconf_123 = shared_int[72];
  int dim_rel_susceptibility = shared_int[73];
  int dim_rel_susceptibility_1 = shared_int[74];
  int dim_rel_susceptibility_2 = shared_int[75];
  int dim_vaccine_dose_step_1 = shared_int[76];
  int dim_vaccine_dose_step_12 = shared_int[77];
  int dim_vaccine_dose_step_2 = shared_int[78];
  int dim_vaccine_dose_step_3 = shared_int[79];
  int offset_variable_E = shared_int[80];
  int offset_variable_G_D = shared_int[81];
  int offset_variable_H_D_conf = shared_int[82];
  int offset_variable_H_D_unconf = shared_int[83];
  int offset_variable_H_R_conf = shared_int[84];
  int offset_variable_H_R_unconf = shared_int[85];
  int offset_variable_ICU_D_conf = shared_int[86];
  int offset_variable_ICU_D_unconf = shared_int[87];
  int offset_variable_ICU_W_D_conf = shared_int[88];
  int offset_variable_ICU_W_D_unconf = shared_int[89];
  int offset_variable_ICU_W_R_conf = shared_int[90];
  int offset_variable_ICU_W_R_unconf = shared_int[91];
  int offset_variable_ICU_pre_conf = shared_int[92];
  int offset_variable_ICU_pre_unconf = shared_int[93];
  int offset_variable_I_A = shared_int[94];
  int offset_variable_I_C_1 = shared_int[95];
  int offset_variable_I_C_2 = shared_int[96];
  int offset_variable_I_P = shared_int[97];
  int offset_variable_I_weighted = shared_int[98];
  int offset_variable_R = shared_int[99];
  int offset_variable_S = shared_int[100];
  int offset_variable_T_PCR_neg = shared_int[101];
  int offset_variable_T_PCR_pos = shared_int[102];
  int offset_variable_T_PCR_pre = shared_int[103];
  int offset_variable_T_sero_neg = shared_int[104];
  int offset_variable_T_sero_pos = shared_int[105];
  int offset_variable_T_sero_pre = shared_int[106];
  int offset_variable_W_D_conf = shared_int[107];
  int offset_variable_W_D_unconf = shared_int[108];
  int offset_variable_W_R_conf = shared_int[109];
  int offset_variable_W_R_unconf = shared_int[110];
  int offset_variable_cum_n_E_vaccinated = shared_int[111];
  int offset_variable_cum_n_I_A_vaccinated = shared_int[112];
  int offset_variable_cum_n_I_P_vaccinated = shared_int[113];
  int offset_variable_cum_n_R_vaccinated = shared_int[114];
  int offset_variable_cum_n_S_vaccinated = shared_int[115];
  int offset_variable_cum_n_vaccinated = shared_int[116];
  int offset_variable_prob_strain = shared_int[117];
  int offset_variable_tmp_vaccine_probability = shared_int[118];
  const int * index_dose = shared_int + 119;
  real_t dt = shared_real[0];
  real_t I_A_transmission = shared_real[1];
  real_t I_P_transmission = shared_real[2];
  real_t I_C_1_transmission = shared_real[3];
  real_t I_C_2_transmission = shared_real[4];
  real_t hosp_transmission = shared_real[5];
  real_t ICU_transmission = shared_real[6];
  real_t G_D_transmission = shared_real[7];
  real_t p_E_progress = shared_real[8];
  real_t p_G_D_progress = shared_real[9];
  real_t p_H_D_progress = shared_real[10];
  real_t p_H_R_progress = shared_real[11];
  real_t p_ICU_D_progress = shared_real[12];
  real_t p_ICU_W_D_progress = shared_real[13];
  real_t p_ICU_W_R_progress = shared_real[14];
  real_t p_ICU_pre_progress = shared_real[15];
  real_t p_I_A_progress = shared_real[16];
  real_t p_I_C_1_progress = shared_real[17];
  real_t p_I_C_2_progress = shared_real[18];
  real_t p_I_P_progress = shared_real[19];
  real_t p_T_PCR_pos_progress = shared_real[20];
  real_t p_T_PCR_pre_progress = shared_real[21];
  real_t p_T_sero_pos_progress = shared_real[22];
  real_t p_W_D_progress = shared_real[23];
  real_t p_W_R_progress = shared_real[24];
  real_t p_test = shared_real[25];
  real_t model_pcr_and_serology = shared_real[26];
  real_t p_sero_pre_1 = shared_real[27];
  const real_t * beta_step = shared_real + 28;
  const real_t * p_G_D_step = beta_step + dim_beta_step;
  const real_t * p_H_step = p_G_D_step + dim_p_G_D_step;
  const real_t * p_H_D_step = p_H_step + dim_p_H_step;
  const real_t * p_ICU_step = p_H_D_step + dim_p_H_D_step;
  const real_t * p_ICU_D_step = p_ICU_step + dim_p_ICU_step;
  const real_t * p_W_D_step = p_ICU_D_step + dim_p_ICU_D_step;
  const real_t * p_star_step = p_W_D_step + dim_p_W_D_step;
  const real_t * strain_seed_step = p_star_step + dim_p_star_step;
  const real_t * psi_G_D = strain_seed_step + dim_strain_seed_step;
  const real_t * psi_H_D = psi_G_D + 19;
  const real_t * psi_H = psi_H_D + 19;
  const real_t * psi_ICU_D = psi_H + 19;
  const real_t * psi_ICU = psi_ICU_D + 19;
  const real_t * psi_W_D = psi_ICU + 19;
  const real_t * psi_star = psi_W_D + 19;
  const real_t * strain_transmission = psi_star + 19;
  const real_t * p_RS = strain_transmission + dim_strain_transmission;
  const real_t * p_sero_pos = p_RS + 19;
  const real_t * p_C = p_sero_pos + 19;
  const real_t * rel_infectivity = p_C + 19;
  const real_t * rel_p_hosp_if_sympt = rel_infectivity + dim_cum_n_S_vaccinated;
  const real_t * m = rel_p_hosp_if_sympt + dim_cum_n_S_vaccinated;
  const real_t * vaccine_progression_rate_base = m + 361;
  const real_t * rel_susceptibility = vaccine_progression_rate_base + dim_cum_n_S_vaccinated;
  const real_t * rel_p_sympt = rel_susceptibility + dim_rel_susceptibility;
  const real_t * vaccine_dose_step = rel_p_sympt + dim_cum_n_S_vaccinated;
  const real_t * p_T_sero_pre_progress = vaccine_dose_step + dim_vaccine_dose_step;
  dust::interleaved<carehomes::real_t> I_weighted_strain = internal_real + 0;
  dust::interleaved<carehomes::real_t> I_with_diff_trans = I_weighted_strain + dim_R;
  dust::interleaved<carehomes::real_t> aux_E = I_with_diff_trans + dim_R;
  dust::interleaved<carehomes::real_t> aux_G_D = aux_E + dim_E;
  dust::interleaved<carehomes::real_t> aux_H_D_conf = aux_G_D + dim_G_D;
  dust::interleaved<carehomes::real_t> aux_H_D_unconf = aux_H_D_conf + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> aux_H_R_conf = aux_H_D_unconf + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> aux_H_R_unconf = aux_H_R_conf + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_D_conf = aux_H_R_unconf + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_D_unconf = aux_ICU_D_conf + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_W_D_conf = aux_ICU_D_unconf + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_W_D_unconf = aux_ICU_W_D_conf + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_W_R_conf = aux_ICU_W_D_unconf + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_W_R_unconf = aux_ICU_W_R_conf + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_pre_conf = aux_ICU_W_R_unconf + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> aux_ICU_pre_unconf = aux_ICU_pre_conf + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> aux_I_A = aux_ICU_pre_unconf + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> aux_I_C_1 = aux_I_A + dim_I_A;
  dust::interleaved<carehomes::real_t> aux_I_C_2 = aux_I_C_1 + dim_I_C_1;
  dust::interleaved<carehomes::real_t> aux_I_P = aux_I_C_2 + dim_I_C_2;
  dust::interleaved<carehomes::real_t> aux_W_D_conf = aux_I_P + dim_I_P;
  dust::interleaved<carehomes::real_t> aux_W_D_unconf = aux_W_D_conf + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> aux_W_R_conf = aux_W_D_unconf + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> aux_W_R_unconf = aux_W_R_conf + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> delta_D_hosp = aux_W_R_unconf + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> delta_D_non_hosp = delta_D_hosp + 19;
  dust::interleaved<carehomes::real_t> lambda = delta_D_non_hosp + 19;
  dust::interleaved<carehomes::real_t> n_EE = lambda + dim_E_12;
  dust::interleaved<carehomes::real_t> n_EE_next_vacc_class = n_EE + dim_E;
  dust::interleaved<carehomes::real_t> n_EI_A = n_EE_next_vacc_class + dim_E;
  dust::interleaved<carehomes::real_t> n_EI_A_next_vacc_class = n_EI_A + dim_R;
  dust::interleaved<carehomes::real_t> n_EI_P = n_EI_A_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> n_EI_P_next_vacc_class = n_EI_P + dim_R;
  dust::interleaved<carehomes::real_t> n_E_next_vacc_class = n_EI_P_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> n_E_progress = n_E_next_vacc_class + dim_E;
  dust::interleaved<carehomes::real_t> n_G_D_progress = n_E_progress + dim_E;
  dust::interleaved<carehomes::real_t> n_H_D_conf_progress = n_G_D_progress + dim_G_D;
  dust::interleaved<carehomes::real_t> n_H_D_unconf_progress = n_H_D_conf_progress + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> n_H_D_unconf_to_conf = n_H_D_unconf_progress + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> n_H_R_conf_progress = n_H_D_unconf_to_conf + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> n_H_R_unconf_progress = n_H_R_conf_progress + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> n_H_R_unconf_to_conf = n_H_R_unconf_progress + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_D_conf_progress = n_H_R_unconf_to_conf + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_D_unconf_progress = n_ICU_D_conf_progress + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_D_unconf_to_conf = n_ICU_D_unconf_progress + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_D_conf_progress = n_ICU_D_unconf_to_conf + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_D_unconf_progress = n_ICU_W_D_conf_progress + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_D_unconf_to_conf = n_ICU_W_D_unconf_progress + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_R_conf_progress = n_ICU_W_D_unconf_to_conf + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_R_unconf_progress = n_ICU_W_R_conf_progress + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_W_R_unconf_to_conf = n_ICU_W_R_unconf_progress + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_pre_conf_progress = n_ICU_W_R_unconf_to_conf + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_pre_conf_to_ICU_D_conf = n_ICU_pre_conf_progress + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_pre_conf_to_ICU_W_D_conf = n_ICU_pre_conf_to_ICU_D_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_ICU_pre_conf_to_ICU_W_R_conf = n_ICU_pre_conf_to_ICU_W_D_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_ICU_pre_unconf_progress = n_ICU_pre_conf_to_ICU_W_R_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_ICU_pre_unconf_to_ICU_D_unconf = n_ICU_pre_unconf_progress + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> n_ICU_pre_unconf_to_ICU_W_D_unconf = n_ICU_pre_unconf_to_ICU_D_unconf + dim_R;
  dust::interleaved<carehomes::real_t> n_ICU_pre_unconf_to_ICU_W_R_unconf = n_ICU_pre_unconf_to_ICU_W_D_unconf + dim_R;
  dust::interleaved<carehomes::real_t> n_ICU_pre_unconf_to_conf = n_ICU_pre_unconf_to_ICU_W_R_unconf + dim_R;
  dust::interleaved<carehomes::real_t> n_II_A = n_ICU_pre_unconf_to_conf + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> n_II_A_next_vacc_class = n_II_A + dim_I_A;
  dust::interleaved<carehomes::real_t> n_II_P = n_II_A_next_vacc_class + dim_I_A;
  dust::interleaved<carehomes::real_t> n_II_P_next_vacc_class = n_II_P + dim_I_P;
  dust::interleaved<carehomes::real_t> n_I_A_next_vacc_class = n_II_P_next_vacc_class + dim_I_P;
  dust::interleaved<carehomes::real_t> n_I_A_progress = n_I_A_next_vacc_class + dim_I_A;
  dust::interleaved<carehomes::real_t> n_I_C_1_progress = n_I_A_progress + dim_I_A;
  dust::interleaved<carehomes::real_t> n_I_C_2_progress = n_I_C_1_progress + dim_I_C_1;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_G_D = n_I_C_2_progress + dim_I_C_2;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_H_D = n_I_C_2_to_G_D + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_H_D_conf = n_I_C_2_to_H_D + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_H_R = n_I_C_2_to_H_D_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_H_R_conf = n_I_C_2_to_H_R + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_ICU_pre = n_I_C_2_to_H_R_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_ICU_pre_conf = n_I_C_2_to_ICU_pre + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_R = n_I_C_2_to_ICU_pre_conf + dim_R;
  dust::interleaved<carehomes::real_t> n_I_C_2_to_hosp = n_I_C_2_to_R + dim_R;
  dust::interleaved<carehomes::real_t> n_I_P_next_vacc_class = n_I_C_2_to_hosp + dim_R;
  dust::interleaved<carehomes::real_t> n_I_P_progress = n_I_P_next_vacc_class + dim_I_P;
  dust::interleaved<carehomes::real_t> n_RS = n_I_P_progress + dim_I_P;
  dust::interleaved<carehomes::real_t> n_RS_next_vacc_class = n_RS + dim_R;
  dust::interleaved<carehomes::real_t> n_R_next_vacc_class = n_RS_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> n_R_next_vacc_class_capped = n_R_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> n_R_next_vacc_class_tmp = n_R_next_vacc_class_capped + dim_R;
  dust::interleaved<carehomes::real_t> n_R_progress = n_R_next_vacc_class_tmp + dim_R;
  dust::interleaved<carehomes::real_t> n_R_progress_capped = n_R_progress + dim_R;
  dust::interleaved<carehomes::real_t> n_R_progress_tmp = n_R_progress_capped + dim_R;
  dust::interleaved<carehomes::real_t> n_SE = n_R_progress_tmp + dim_R;
  dust::interleaved<carehomes::real_t> n_SE_next_vacc_class = n_SE + dim_R;
  dust::interleaved<carehomes::real_t> n_S_next_vacc_class = n_SE_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> n_S_progress = n_S_next_vacc_class + dim_cum_n_S_vaccinated;
  dust::interleaved<carehomes::real_t> n_S_progress_tot = n_S_progress + dim_R;
  dust::interleaved<carehomes::real_t> n_T_PCR_pos_progress = n_S_progress_tot + dim_cum_n_S_vaccinated;
  dust::interleaved<carehomes::real_t> n_T_PCR_pre_progress = n_T_PCR_pos_progress + dim_T_PCR_pos;
  dust::interleaved<carehomes::real_t> n_T_sero_pos_progress = n_T_PCR_pre_progress + dim_T_PCR_pre;
  dust::interleaved<carehomes::real_t> n_T_sero_pre_progress = n_T_sero_pos_progress + dim_T_sero_pos;
  dust::interleaved<carehomes::real_t> n_T_sero_pre_to_T_sero_pos = n_T_sero_pre_progress + dim_T_sero_pre;
  dust::interleaved<carehomes::real_t> n_W_D_conf_progress = n_T_sero_pre_to_T_sero_pos + dim_R;
  dust::interleaved<carehomes::real_t> n_W_D_unconf_progress = n_W_D_conf_progress + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_W_D_unconf_to_conf = n_W_D_unconf_progress + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_W_R_conf_progress = n_W_D_unconf_to_conf + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> n_W_R_unconf_progress = n_W_R_conf_progress + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_W_R_unconf_to_conf = n_W_R_unconf_progress + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_com_to_T_sero_pre = n_W_R_unconf_to_conf + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> n_hosp_non_ICU = n_com_to_T_sero_pre + dim_T_sero_pre;
  dust::interleaved<carehomes::real_t> new_E = n_hosp_non_ICU + dim_R;
  dust::interleaved<carehomes::real_t> new_G_D = new_E + dim_E;
  dust::interleaved<carehomes::real_t> new_H_D_conf = new_G_D + dim_G_D;
  dust::interleaved<carehomes::real_t> new_H_D_unconf = new_H_D_conf + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> new_H_R_conf = new_H_D_unconf + dim_H_D_unconf;
  dust::interleaved<carehomes::real_t> new_H_R_unconf = new_H_R_conf + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_D_conf = new_H_R_unconf + dim_H_R_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_D_unconf = new_ICU_D_conf + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_W_D_conf = new_ICU_D_unconf + dim_ICU_D_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_W_D_unconf = new_ICU_W_D_conf + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_W_R_conf = new_ICU_W_D_unconf + dim_ICU_W_D_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_W_R_unconf = new_ICU_W_R_conf + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_pre_conf = new_ICU_W_R_unconf + dim_ICU_W_R_unconf;
  dust::interleaved<carehomes::real_t> new_ICU_pre_unconf = new_ICU_pre_conf + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> new_I_A = new_ICU_pre_unconf + dim_ICU_pre_unconf;
  dust::interleaved<carehomes::real_t> new_I_C_1 = new_I_A + dim_I_A;
  dust::interleaved<carehomes::real_t> new_I_C_2 = new_I_C_1 + dim_I_C_1;
  dust::interleaved<carehomes::real_t> new_I_P = new_I_C_2 + dim_I_C_2;
  dust::interleaved<carehomes::real_t> new_R = new_I_P + dim_I_P;
  dust::interleaved<carehomes::real_t> new_S = new_R + dim_R;
  dust::interleaved<carehomes::real_t> new_T_PCR_neg = new_S + dim_cum_n_S_vaccinated;
  dust::interleaved<carehomes::real_t> new_T_PCR_pos = new_T_PCR_neg + dim_R;
  dust::interleaved<carehomes::real_t> new_T_PCR_pre = new_T_PCR_pos + dim_T_PCR_pos;
  dust::interleaved<carehomes::real_t> new_T_sero_neg = new_T_PCR_pre + dim_T_PCR_pre;
  dust::interleaved<carehomes::real_t> new_T_sero_pos = new_T_sero_neg + dim_R;
  dust::interleaved<carehomes::real_t> new_T_sero_pre = new_T_sero_pos + dim_T_sero_pos;
  dust::interleaved<carehomes::real_t> new_W_D_conf = new_T_sero_pre + dim_T_sero_pre;
  dust::interleaved<carehomes::real_t> new_W_D_unconf = new_W_D_conf + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> new_W_R_conf = new_W_D_unconf + dim_W_D_unconf;
  dust::interleaved<carehomes::real_t> new_W_R_unconf = new_W_R_conf + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> p_E_next_vacc_class = new_W_R_unconf + dim_W_R_unconf;
  dust::interleaved<carehomes::real_t> p_G_D_by_age = p_E_next_vacc_class + dim_E;
  dust::interleaved<carehomes::real_t> p_H_D_by_age = p_G_D_by_age + 19;
  dust::interleaved<carehomes::real_t> p_H_by_age = p_H_D_by_age + 19;
  dust::interleaved<carehomes::real_t> p_ICU_D_by_age = p_H_by_age + 19;
  dust::interleaved<carehomes::real_t> p_ICU_by_age = p_ICU_D_by_age + 19;
  dust::interleaved<carehomes::real_t> p_I_A_next_vacc_class = p_ICU_by_age + 19;
  dust::interleaved<carehomes::real_t> p_I_P_next_vacc_class = p_I_A_next_vacc_class + dim_I_A;
  dust::interleaved<carehomes::real_t> p_R_next_vacc_class = p_I_P_next_vacc_class + dim_I_P;
  dust::interleaved<carehomes::real_t> p_SE = p_R_next_vacc_class + dim_R;
  dust::interleaved<carehomes::real_t> p_S_next_vacc_class = p_SE + dim_cum_n_S_vaccinated;
  dust::interleaved<carehomes::real_t> p_W_D_by_age = p_S_next_vacc_class + dim_cum_n_S_vaccinated;
  dust::interleaved<carehomes::real_t> p_star_by_age = p_W_D_by_age + 19;
  dust::interleaved<carehomes::real_t> s_ij = p_star_by_age + 19;
  dust::interleaved<carehomes::real_t> vaccine_n_candidates = s_ij + dim_s_ij;
  dust::interleaved<carehomes::real_t> vaccine_probability = vaccine_n_candidates + 38;
  dust::interleaved<carehomes::real_t> vaccine_probability_doses = vaccine_probability + dim_cum_n_S_vaccinated;
  const dust::interleaved<real_t> cum_n_S_vaccinated = state + offset_variable_cum_n_S_vaccinated;
  const dust::interleaved<real_t> cum_n_E_vaccinated = state + offset_variable_cum_n_E_vaccinated;
  const dust::interleaved<real_t> cum_n_I_A_vaccinated = state + offset_variable_cum_n_I_A_vaccinated;
  const dust::interleaved<real_t> cum_n_I_P_vaccinated = state + offset_variable_cum_n_I_P_vaccinated;
  const dust::interleaved<real_t> cum_n_R_vaccinated = state + offset_variable_cum_n_R_vaccinated;
  const dust::interleaved<real_t> S = state + offset_variable_S;
  const dust::interleaved<real_t> E = state + offset_variable_E;
  const dust::interleaved<real_t> I_A = state + offset_variable_I_A;
  const dust::interleaved<real_t> I_P = state + offset_variable_I_P;
  const dust::interleaved<real_t> I_C_1 = state + offset_variable_I_C_1;
  const dust::interleaved<real_t> I_C_2 = state + offset_variable_I_C_2;
  const dust::interleaved<real_t> G_D = state + offset_variable_G_D;
  const dust::interleaved<real_t> ICU_pre_unconf = state + offset_variable_ICU_pre_unconf;
  const dust::interleaved<real_t> ICU_pre_conf = state + offset_variable_ICU_pre_conf;
  const dust::interleaved<real_t> H_R_unconf = state + offset_variable_H_R_unconf;
  const dust::interleaved<real_t> H_R_conf = state + offset_variable_H_R_conf;
  const dust::interleaved<real_t> H_D_unconf = state + offset_variable_H_D_unconf;
  const dust::interleaved<real_t> H_D_conf = state + offset_variable_H_D_conf;
  const dust::interleaved<real_t> ICU_W_R_unconf = state + offset_variable_ICU_W_R_unconf;
  const dust::interleaved<real_t> ICU_W_R_conf = state + offset_variable_ICU_W_R_conf;
  const dust::interleaved<real_t> ICU_W_D_unconf = state + offset_variable_ICU_W_D_unconf;
  const dust::interleaved<real_t> ICU_W_D_conf = state + offset_variable_ICU_W_D_conf;
  const dust::interleaved<real_t> ICU_D_unconf = state + offset_variable_ICU_D_unconf;
  const dust::interleaved<real_t> ICU_D_conf = state + offset_variable_ICU_D_conf;
  const dust::interleaved<real_t> W_R_unconf = state + offset_variable_W_R_unconf;
  const dust::interleaved<real_t> W_R_conf = state + offset_variable_W_R_conf;
  const dust::interleaved<real_t> W_D_unconf = state + offset_variable_W_D_unconf;
  const dust::interleaved<real_t> W_D_conf = state + offset_variable_W_D_conf;
  const dust::interleaved<real_t> T_sero_pre = state + offset_variable_T_sero_pre;
  const dust::interleaved<real_t> T_sero_pos = state + offset_variable_T_sero_pos;
  const dust::interleaved<real_t> T_sero_neg = state + offset_variable_T_sero_neg;
  const dust::interleaved<real_t> R = state + offset_variable_R;
  const dust::interleaved<real_t> D_hosp = state + 27;
  const dust::interleaved<real_t> D_non_hosp = state + 46;
  const dust::interleaved<real_t> T_PCR_pre = state + offset_variable_T_PCR_pre;
  const dust::interleaved<real_t> T_PCR_pos = state + offset_variable_T_PCR_pos;
  const dust::interleaved<real_t> T_PCR_neg = state + offset_variable_T_PCR_neg;
  const real_t cum_admit_conf = state[4];
  const real_t cum_new_conf = state[5];
  const real_t admit_conf_inc = state[1];
  const real_t new_conf_inc = state[2];
  const dust::interleaved<real_t> cum_admit_by_age = state + 65;
  const real_t cum_infections = state[3];
  const dust::interleaved<real_t> cum_infections_per_strain = state + 141;
  const real_t D_hosp_tot = state[12];
  const real_t D_comm_tot = state[13];
  const real_t D_comm_inc = state[14];
  const real_t D_carehomes_tot = state[15];
  const real_t D_carehomes_inc = state[16];
  const real_t D_hosp_inc = state[17];
  const real_t D_tot = state[18];
  const real_t cum_sympt_cases = state[20];
  const real_t cum_sympt_cases_over25 = state[21];
  const real_t cum_sympt_cases_non_variant_over25 = state[22];
  const real_t sympt_cases_inc = state[23];
  const real_t sympt_cases_over25_inc = state[24];
  const real_t sympt_cases_non_variant_over25_inc = state[25];
  state_next[7] = odin_sum1<real_t>(S, 0, dim_cum_n_S_vaccinated) + odin_sum1<real_t>(T_sero_pre, 0, dim_T_sero_pre) + odin_sum1<real_t>(T_sero_pos, 0, dim_T_sero_pos) + odin_sum1<real_t>(T_sero_neg, 0, dim_R) + odin_sum1<real_t>(E, 0, dim_E);
  state_next[8] = odin_sum1<real_t>(S, 0, dim_cum_n_S_vaccinated) + odin_sum1<real_t>(T_PCR_pre, 0, dim_T_PCR_pre) + odin_sum1<real_t>(T_PCR_pos, 0, dim_T_PCR_pos) + odin_sum1<real_t>(T_PCR_neg, 0, dim_R);
  real_t beta = (static_cast<int>(step) >= dim_beta_step ? beta_step[dim_beta_step - 1] : beta_step[step + 1 - 1]);
  real_t p_G_D = (static_cast<int>(step) >= dim_p_G_D_step ? p_G_D_step[dim_p_G_D_step - 1] : p_G_D_step[step + 1 - 1]);
  real_t p_H = (static_cast<int>(step) >= dim_p_H_step ? p_H_step[dim_p_H_step - 1] : p_H_step[step + 1 - 1]);
  real_t p_H_D = (static_cast<int>(step) >= dim_p_H_D_step ? p_H_D_step[dim_p_H_D_step - 1] : p_H_D_step[step + 1 - 1]);
  real_t p_ICU = (static_cast<int>(step) >= dim_p_ICU_step ? p_ICU_step[dim_p_ICU_step - 1] : p_ICU_step[step + 1 - 1]);
  real_t p_ICU_D = (static_cast<int>(step) >= dim_p_ICU_D_step ? p_ICU_D_step[dim_p_ICU_D_step - 1] : p_ICU_D_step[step + 1 - 1]);
  real_t p_W_D = (static_cast<int>(step) >= dim_p_W_D_step ? p_W_D_step[dim_p_W_D_step - 1] : p_W_D_step[step + 1 - 1]);
  real_t p_star = (static_cast<int>(step) >= dim_p_star_step ? p_star_step[dim_p_star_step - 1] : p_star_step[step + 1 - 1]);
  real_t strain_seed = ((static_cast<int>(step) >= dim_strain_seed_step ? strain_seed_step[dim_strain_seed_step - 1] : strain_seed_step[step + 1 - 1]));
  state_next[0] = (step + 1) * dt;
  for (int i = 1; i <= 19; ++i) {
    p_G_D_by_age[i - 1] = p_G_D * psi_G_D[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_H_D_by_age[i - 1] = p_H_D * psi_H_D[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_H_by_age[i - 1] = p_H * psi_H[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_ICU_D_by_age[i - 1] = p_ICU_D * psi_ICU_D[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_ICU_by_age[i - 1] = p_ICU * psi_ICU[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_W_D_by_age[i - 1] = p_W_D * psi_W_D[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    p_star_by_age[i - 1] = p_star * psi_star[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    state_next[84 + i - 1] = odin_sum2<real_t>(S, i - 1, i, 0, dim_rel_susceptibility_2, 19) + odin_sum3<real_t>(R, i - 1, i, 0, dim_strain_transmission, 0, dim_rel_susceptibility_2, 19, dim_E_12) + D_hosp[i - 1] + odin_sum4<real_t>(E, i - 1, i, 0, dim_strain_transmission, 0, k_E, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_E_123) + odin_sum4<real_t>(I_A, i - 1, i, 0, dim_strain_transmission, 0, k_A, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_I_A_123) + odin_sum4<real_t>(I_P, i - 1, i, 0, dim_strain_transmission, 0, k_P, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_I_P_123) + odin_sum4<real_t>(I_C_1, i - 1, i, 0, dim_strain_transmission, 0, k_C_1, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_I_C_1_123) + odin_sum4<real_t>(I_C_2, i - 1, i, 0, dim_strain_transmission, 0, k_C_2, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_I_C_2_123) + odin_sum4<real_t>(ICU_pre_conf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_pre, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(ICU_pre_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_pre, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(H_R_conf, i - 1, i, 0, dim_strain_transmission, 0, k_H_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(H_R_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_H_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(H_D_conf, i - 1, i, 0, dim_strain_transmission, 0, k_H_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(H_D_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_H_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(ICU_W_R_conf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_W_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_R_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_W_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_D_conf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_W_D_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_D_conf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(ICU_D_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_ICU_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(W_R_conf, i - 1, i, 0, dim_strain_transmission, 0, k_W_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_R_unconf_123) + odin_sum4<real_t>(W_R_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_W_R, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_R_unconf_123) + odin_sum4<real_t>(W_D_conf, i - 1, i, 0, dim_strain_transmission, 0, k_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_D_unconf_123) + odin_sum4<real_t>(W_D_unconf, i - 1, i, 0, dim_strain_transmission, 0, k_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_D_unconf_123) + odin_sum4<real_t>(G_D, i - 1, i, 0, dim_strain_transmission, 0, k_G_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_G_D_123) + D_non_hosp[i - 1];
  }
  state_next[6] = beta;
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_S_vaccinated[19 * (j - 1) + i - 1] + cum_n_E_vaccinated[19 * (j - 1) + i - 1] + cum_n_I_A_vaccinated[19 * (j - 1) + i - 1] + cum_n_I_P_vaccinated[19 * (j - 1) + i - 1] + cum_n_R_vaccinated[19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= 2; ++j) {
      vaccine_n_candidates[i - 1 + 19 * (j - 1)] = S[19 * (index_dose[j - 1] - 1) + i - 1] + odin_sum4<real_t>(E, i - 1, i, 0, dim_strain_transmission, 0, k_E, index_dose[j - 1] - 1, index_dose[j - 1], 19, dim_E_12, dim_E_123) + odin_sum4<real_t>(I_A, i - 1, i, 0, dim_strain_transmission, 0, k_A, index_dose[j - 1] - 1, index_dose[j - 1], 19, dim_E_12, dim_I_A_123) + odin_sum4<real_t>(I_P, i - 1, i, 0, dim_strain_transmission, 0, k_P, index_dose[j - 1] - 1, index_dose[j - 1], 19, dim_E_12, dim_I_P_123) + odin_sum3<real_t>(R, i - 1, i, 0, dim_strain_transmission, index_dose[j - 1] - 1, index_dose[j - 1], 19, dim_E_12);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        I_with_diff_trans[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = rel_infectivity[19 * (k - 1) + i - 1] * strain_transmission[j - 1] * (I_A_transmission * odin_sum4<real_t>(I_A, i - 1, i, j - 1, j, 0, k_A, k - 1, k, 19, dim_E_12, dim_I_A_123) + I_P_transmission * odin_sum4<real_t>(I_P, i - 1, i, j - 1, j, 0, k_P, k - 1, k, 19, dim_E_12, dim_I_P_123) + I_C_1_transmission * odin_sum4<real_t>(I_C_1, i - 1, i, j - 1, j, 0, k_C_1, k - 1, k, 19, dim_E_12, dim_I_C_1_123) + I_C_2_transmission * odin_sum4<real_t>(I_C_2, i - 1, i, j - 1, j, 0, k_C_2, k - 1, k, 19, dim_E_12, dim_I_C_2_123) + hosp_transmission * (odin_sum4<real_t>(ICU_pre_unconf, i - 1, i, j - 1, j, 0, k_ICU_pre, k - 1, k, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(ICU_pre_conf, i - 1, i, j - 1, j, 0, k_ICU_pre, k - 1, k, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(H_R_unconf, i - 1, i, j - 1, j, 0, k_H_R, k - 1, k, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(H_R_conf, i - 1, i, j - 1, j, 0, k_H_R, k - 1, k, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(H_D_unconf, i - 1, i, j - 1, j, 0, k_H_D, k - 1, k, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(H_D_conf, i - 1, i, j - 1, j, 0, k_H_D, k - 1, k, 19, dim_E_12, dim_H_D_unconf_123)) + ICU_transmission * (odin_sum4<real_t>(ICU_W_R_unconf, i - 1, i, j - 1, j, 0, k_ICU_W_R, k - 1, k, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_R_conf, i - 1, i, j - 1, j, 0, k_ICU_W_R, k - 1, k, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(ICU_W_D_unconf, i - 1, i, j - 1, j, 0, k_ICU_W_D, k - 1, k, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_W_D_conf, i - 1, i, j - 1, j, 0, k_ICU_W_D, k - 1, k, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(ICU_D_unconf, i - 1, i, j - 1, j, 0, k_ICU_D, k - 1, k, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(ICU_D_conf, i - 1, i, j - 1, j, 0, k_ICU_D, k - 1, k, 19, dim_E_12, dim_ICU_D_unconf_123)) + G_D_transmission * odin_sum4<real_t>(G_D, i - 1, i, j - 1, j, 0, k_G_D, k - 1, k, 19, dim_E_12, dim_G_D_123));
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_E_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_E_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_G_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_G_D_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(G_D[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_G_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_D_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_H_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_D_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_H_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_R_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_H_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_R_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_H_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_D_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_D_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_D_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_W_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_D_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_W_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_R_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_W_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_R_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_W_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_pre_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_pre_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_pre_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_pre_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_A_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_A_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_1; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_C_1_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_C_1[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_C_1_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_C_2_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_C_2[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_C_2_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_P_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_P_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_progress_tmp[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_RS[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_T_PCR_pos_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pos_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_PCR_pos[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_T_PCR_pos_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_T_PCR_pre_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_PCR_pre[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_T_PCR_pre_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_sero_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_T_sero_pos_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pos_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_sero_pos[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_T_sero_pos_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_D_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_D_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_D_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_R_conf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_R_unconf_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_R_progress);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= 2; ++j) {
      state_next[103 + i - 1 + 19 * (j - 1)] = vaccine_n_candidates[19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= 2; ++j) {
      vaccine_probability_doses[i - 1 + 19 * (j - 1)] = ((static_cast<int>(step) >= dim_vaccine_dose_step_3 || vaccine_n_candidates[19 * (j - 1) + i - 1] == 0 ? 0 : vaccine_dose_step[dim_vaccine_dose_step_12 * (step + 1 - 1) + dim_vaccine_dose_step_1 * (j - 1) + i - 1] / (real_t) vaccine_n_candidates[19 * (j - 1) + i - 1]));
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_D_conf_progress[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_D_conf_progress[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_D_unconf_progress[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_D_unconf_progress[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_R_conf_progress[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_R_conf_progress[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_R_unconf_progress[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_H_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_R_unconf_progress[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_conf_progress[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_conf_progress[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_unconf_progress[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_pre_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_unconf_progress[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_I_C_2[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = n_I_C_1_progress[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k_C_1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_C_2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_C_2[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = n_I_C_2_progress[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_C_2[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = aux_I_C_2[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_progress[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_D_conf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k_ICU_W_D - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_D_conf_progress[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_D_conf_progress[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_D_unconf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k_ICU_W_D - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_D_unconf_progress[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_D_unconf_progress[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_R_conf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k_ICU_W_R - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_R_conf_progress[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_R_conf_progress[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_R_unconf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k_ICU_W_R - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_R_unconf_progress[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_R_unconf_progress[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    delta_D_hosp[i - 1] = odin_sum4<real_t>(n_H_D_unconf_progress, i - 1, i, 0, dim_strain_transmission, k_H_D - 1, k_H_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(n_H_D_conf_progress, i - 1, i, 0, dim_strain_transmission, k_H_D - 1, k_H_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(n_ICU_D_unconf_progress, i - 1, i, 0, dim_strain_transmission, k_ICU_D - 1, k_ICU_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(n_ICU_D_conf_progress, i - 1, i, 0, dim_strain_transmission, k_ICU_D - 1, k_ICU_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(n_W_D_unconf_progress, i - 1, i, 0, dim_strain_transmission, k_W_D - 1, k_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_D_unconf_123) + odin_sum4<real_t>(n_W_D_conf_progress, i - 1, i, 0, dim_strain_transmission, k_W_D - 1, k_W_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_W_D_unconf_123);
  }
  for (int i = 1; i <= 19; ++i) {
    delta_D_non_hosp[i - 1] = odin_sum4<real_t>(n_G_D_progress, i - 1, i, 0, dim_strain_transmission, k_G_D - 1, k_G_D, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_G_D_123);
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_conf_to_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_ICU_pre_conf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1]), p_ICU_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_unconf_to_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_ICU_pre_unconf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1]), p_ICU_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_R[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_progress[dim_I_C_2_123 * (k - 1) + dim_E_12 * (k_C_2 - 1) + 19 * (j - 1) + i - 1]), 1 - p_H_by_age[i - 1] * rel_p_hosp_if_sympt[19 * (k - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_progress_capped[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = odin_min(n_R_progress_tmp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], odin_min(T_sero_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], T_PCR_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]));
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= 2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_T_sero_pre_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(T_sero_pre[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_T_sero_pre_progress[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_PCR_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pos_123 * (l - 1)] = T_PCR_pos[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_T_PCR_pos_progress[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_T_PCR_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pos_123 * (l - 1)] = new_T_PCR_pos[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_T_PCR_pre_progress[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k_PCR_pre - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_PCR_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_PCR_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pos_123 * (l - 1)] = new_T_PCR_pos[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_T_PCR_pos_progress[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= 19; ++j) {
      for (int k = 1; k <= dim_strain_transmission; ++k) {
        s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = m[19 * (j - 1) + i - 1] * odin_sum3<real_t>(I_with_diff_trans, j - 1, j, k - 1, k, 0, dim_rel_susceptibility_2, 19, dim_E_12);
      }
    }
  }
  for (int i = 1; i <= n_age_groups; ++i) {
    for (int j = 1; j <= n_groups; ++j) {
      for (int k = 1; k <= dim_strain_transmission; ++k) {
        s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = beta * s_ij[361 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = (n_age_groups + 1); i <= n_groups; ++i) {
    for (int j = 1; j <= n_age_groups; ++j) {
      for (int k = 1; k <= dim_strain_transmission; ++k) {
        s_ij[i - 1 + 19 * (j - 1) + 361 * (k - 1)] = beta * s_ij[361 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      vaccine_probability[i - 1 + 19 * (j - 1)] = 1 - std::exp(- vaccine_progression_rate_base[19 * (j - 1) + i - 1] * dt);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    int j = index_dose[0];
    vaccine_probability[i - 1 + 19 * (j - 1)] = vaccine_probability_doses[19 * 0 + i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    int j = index_dose[1];
    vaccine_probability[i - 1 + 19 * (j - 1)] = vaccine_probability_doses[19 * 1 + i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_conf_to_ICU_D_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_D_conf_progress[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_D_conf_progress[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_unconf_to_ICU_D_unconf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_D_unconf_progress[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_D_unconf_progress[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  real_t delta_D_carehomes_tot = delta_D_non_hosp[18];
  real_t delta_D_comm_tot = odin_sum1<real_t>(delta_D_non_hosp, 0, 18);
  real_t delta_D_hosp_tot = odin_sum1<real_t>(delta_D_hosp, 0, 19);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      lambda[i - 1 + 19 * (j - 1)] = odin_sum3<real_t>(s_ij, i - 1, i, 0, 19, j - 1, j, 19, 361);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_D_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_H_R_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_conf_to_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_ICU_pre_conf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_conf_to_ICU_D_conf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_unconf_to_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_ICU_pre_unconf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_unconf_to_ICU_D_unconf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_W_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_pre_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_G_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_progress[dim_I_C_2_123 * (k - 1) + dim_E_12 * (k_C_2 - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_G_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = (model_pcr_and_serology == 1 ? n_R_progress_capped[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] : n_R_progress_tmp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_T_sero_pre_to_T_sero_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(odin_sum4<real_t>(n_T_sero_pre_progress, i - 1, i, j - 1, j, 0, 2, k - 1, k, 19, dim_E_12, dim_T_sero_pre_123)), p_sero_pos[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_D_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_W_R_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_I_C_2[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = I_C_2[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_I_C_2[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          p_E_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = vaccine_probability[19 * (l - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          p_I_A_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = vaccine_probability[19 * (l - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          p_I_P_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = vaccine_probability[19 * (l - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        p_R_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = vaccine_probability[19 * (k - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      p_S_next_vacc_class[i - 1 + 19 * (j - 1)] = vaccine_probability[19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    state_next[27 + i - 1] = D_hosp[i - 1] + delta_D_hosp[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    state_next[46 + i - 1] = D_non_hosp[i - 1] + delta_D_non_hosp[i - 1];
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_T_PCR_pos + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pos_123 * (l - 1)] = new_T_PCR_pos[dim_T_PCR_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  state_next[26] = odin_sum4<real_t>(new_T_PCR_pos, 1, 18, 0, dim_strain_transmission, 0, k_PCR_pos, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_T_PCR_pos_123);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_tmp_vaccine_probability + i - 1 + 19 * (j - 1)] = vaccine_probability[19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_G_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = n_I_C_2_to_G_D[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_G_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_G_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = n_G_D_progress[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_G_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_G_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = aux_G_D[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_G_D_progress[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_conf_to_ICU_W_D_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_D_conf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_D_conf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_unconf_to_ICU_W_D_unconf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_D_unconf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_D_unconf_progress[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_EE_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(n_E_progress[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_E_next_vacc_class[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_E_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_E_progress[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_E_next_vacc_class[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_D_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_conf_to_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_ICU_pre_conf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_conf_to_ICU_D_conf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_conf_to_ICU_W_D_conf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_ICU_pre_unconf_to_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_ICU_pre_unconf_progress[dim_ICU_pre_unconf_123 * (k - 1) + dim_E_12 * (k_ICU_pre - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_unconf_to_ICU_D_unconf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_unconf_to_ICU_W_D_unconf[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_II_A_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_A_progress[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_A_next_vacc_class[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_II_P_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_P_progress[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_P_next_vacc_class[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_A_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_A_progress[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_A_next_vacc_class[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_hosp[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_I_C_2_progress[dim_I_C_2_123 * (k - 1) + dim_E_12 * (k_C_2 - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_G_D[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_I_P_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_P_progress[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_I_P_next_vacc_class[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_RS_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_R_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_next_vacc_class_tmp[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_R_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_sero_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_sero_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pos_123 * (l - 1)] = T_sero_pos[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_T_sero_pos_progress[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_T_sero_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pos_123 * (l - 1)] = new_T_sero_pos[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_T_sero_pre_to_T_sero_pos[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_sero_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_sero_pos[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pos_123 * (l - 1)] = new_T_sero_pos[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_T_sero_pos_progress[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_D_unconf_to_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = aux_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_D_unconf_to_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_W_R_unconf_to_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = aux_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_W_R_unconf_to_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      p_SE[i - 1 + 19 * (j - 1)] = 1 - std::exp(- odin_sum2<real_t>(lambda, i - 1, i, 0, dim_strain_transmission, 19) * rel_susceptibility[dim_rel_susceptibility_1 * (j - 1) + i - 1] * dt);
    }
  }
  state_next[16] = (fmodr<real_t>(step, steps_per_day) == 0 ? delta_D_carehomes_tot : D_carehomes_inc + delta_D_carehomes_tot);
  state_next[15] = D_carehomes_tot + delta_D_carehomes_tot;
  state_next[14] = (fmodr<real_t>(step, steps_per_day) == 0 ? delta_D_comm_tot : D_comm_inc + delta_D_comm_tot);
  state_next[13] = D_comm_tot + delta_D_comm_tot;
  state_next[17] = (fmodr<real_t>(step, steps_per_day) == 0 ? delta_D_hosp_tot : D_hosp_inc + delta_D_hosp_tot);
  state_next[12] = D_hosp_tot + delta_D_hosp_tot;
  state_next[18] = D_tot + delta_D_hosp_tot + delta_D_comm_tot + delta_D_carehomes_tot;
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_I_C_2 + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_2_123 * (l - 1)] = new_I_C_2[dim_I_C_2_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      state_next[offset_variable_prob_strain + i - 1 + 19 * (j - 1)] = lambda[19 * (j - 1) + i - 1] / (real_t) odin_sum2<real_t>(lambda, i - 1, i, 0, dim_strain_transmission, 19);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_conf_to_ICU_W_R_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_R_conf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_R_conf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_unconf_to_ICU_W_R_unconf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_R_unconf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_R_unconf_progress[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_EE[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = n_E_progress[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_EE_next_vacc_class[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_EI_A_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_EE_next_vacc_class[dim_E_123 * (k - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1]), 1 - p_C[i - 1] * rel_p_sympt[19 * (k - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_D_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_II_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = n_I_A_progress[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_A_next_vacc_class[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_II_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = n_I_P_progress[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_P_next_vacc_class[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_ICU_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_to_hosp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_ICU_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_RS[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_RS_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_next_vacc_class_capped[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = odin_min(n_R_next_vacc_class_tmp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], odin_min(T_sero_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1], T_PCR_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]));
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      n_S_progress_tot[i - 1 + 19 * (j - 1)] = dust::distr::rbinom(rng_state, std::round(S[19 * (j - 1) + i - 1]), p_SE[19 * (j - 1) + i - 1]);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_G_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_G_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = G_D[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_G_D[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_D_unconf_to_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = aux_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_D_unconf_to_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_sero_pos; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_T_sero_pos + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pos_123 * (l - 1)] = new_T_sero_pos[dim_T_sero_pos_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_W_D_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = new_W_D_conf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_W_D_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_D_unconf_123 * (l - 1)] = new_W_D_unconf[dim_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_W_R_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = new_W_R_conf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_W_R_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_W_R_unconf_123 * (l - 1)] = new_W_R_unconf[dim_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    state_next[65 + i - 1] = cum_admit_by_age[i - 1] + odin_sum3<real_t>(n_I_C_2_to_hosp, i - 1, i, 0, dim_strain_transmission, 0, dim_rel_susceptibility_2, 19, dim_E_12);
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_E_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_E_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(n_E_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_E, j - 1, j, 19, dim_E_12, dim_E_123) + odin_sum4<real_t>(n_EE_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_E, j - 1, j, 19, dim_E_12, dim_E_123);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_I_A_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_I_A_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(n_I_A_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_A, j - 1, j, 19, dim_E_12, dim_I_A_123) + odin_sum4<real_t>(n_II_A_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_A, j - 1, j, 19, dim_E_12, dim_I_A_123);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_I_P_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_I_P_vaccinated[19 * (j - 1) + i - 1] + odin_sum4<real_t>(n_I_P_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_P, j - 1, j, 19, dim_E_12, dim_I_P_123) + odin_sum4<real_t>(n_II_P_next_vacc_class, i - 1, i, 0, dim_strain_transmission, 0, k_P, j - 1, j, 19, dim_E_12, dim_I_P_123);
    }
  }
  state_next[19] = odin_sum4<real_t>(new_T_sero_pos, 3, 13, 0, dim_strain_transmission, 0, k_sero_pos, 0, dim_rel_susceptibility_2, 19, dim_E_12, dim_T_sero_pos_123);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      int l = 1;
      aux_I_C_1[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = n_II_P[dim_I_P_123 * 0 + dim_E_12 * (k_P - 1) + 19 * (j - 1) + i - 1] + n_II_P_next_vacc_class[dim_I_P_123 * (n_vacc_classes - 1) + dim_E_12 * (k_P - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        aux_I_C_1[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = n_II_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k_P - 1) + 19 * (j - 1) + i - 1] + n_II_P_next_vacc_class[dim_I_P_123 * (l - 1 - 1) + dim_E_12 * (k_P - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_C_1; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_C_1[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = n_I_C_1_progress[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_1; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_C_1[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = aux_I_C_1[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_C_1_progress[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_EI_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_EE[dim_E_123 * (k - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1]), 1 - p_C[i - 1] * rel_p_sympt[19 * (k - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_EI_P_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_EE_next_vacc_class[dim_E_123 * (k - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] - n_EI_A_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          n_ICU_W_R_unconf_to_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(aux_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_test);
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_to_ICU_pre[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_star_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_R_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = (model_pcr_and_serology == 1 ? n_R_next_vacc_class_capped[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] : n_R_next_vacc_class_tmp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    int j = 1;
    for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
      n_S_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_S_progress_tot[19 * (k - 1) + i - 1]), lambda[19 * 0 + i - 1] / (real_t) odin_sum2<real_t>(lambda, i - 1, i, 0, dim_strain_transmission, 19));
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 2; j <= n_strains; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_S_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_S_progress_tot[19 * (k - 1) + i - 1] - odin_sum3<real_t>(n_S_progress, i - 1, i, 0, j - 1, k - 1, k, 19, dim_E_12)), lambda[19 * (j - 1) + i - 1] / (real_t) odin_sum2<real_t>(lambda, i - 1, i, j - 1, n_strains, 19));
      }
    }
  }
  {
     int i = 4;
     for (int j = 2; j <= n_strains; ++j) {
       int k = 1;
       n_S_progress[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = odin_min(n_S_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + strain_seed, n_S_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + S[19 * (k - 1) + i - 1] - odin_sum3<real_t>(n_S_progress, i - 1, i, 0, dim_strain_transmission, k - 1, k, 19, dim_E_12));
     }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      int l = 1;
      n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(n_EE[dim_E_123 * 0 + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (n_vacc_classes - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1]), p_sero_pre_1);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = dust::distr::rbinom(rng_state, std::round(n_EE[dim_E_123 * (l - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (l - 1 - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1]), p_sero_pre_1);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 2;
      int l = 1;
      n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = n_EE[dim_E_123 * 0 + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (n_vacc_classes - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] - n_com_to_T_sero_pre[dim_T_sero_pre_123 * 0 + dim_E_12 * 0 + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 2;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        n_com_to_T_sero_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = n_EE[dim_E_123 * (l - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (l - 1 - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] - n_com_to_T_sero_pre[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_hosp_non_ICU[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_I_C_2_to_hosp[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_ICU_pre[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_W_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_D_unconf_to_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = aux_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_D_unconf_to_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_G_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_G_D + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_G_D_123 * (l - 1)] = new_G_D[dim_G_D_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_D_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = new_ICU_D_conf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_D_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_D_unconf_123 * (l - 1)] = new_ICU_D_unconf[dim_ICU_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = n_EI_A[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = n_II_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_A_next_vacc_class[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_A_next_vacc_class[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        int l = 1;
        aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_I_A_next_vacc_class[dim_I_A_123 * (n_vacc_classes - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_I_A_next_vacc_class[dim_I_A_123 * (l - 1 - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      int l = 1;
      aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EI_A_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EI_A_next_vacc_class[dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_A; ++k) {
        int l = 1;
        aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_A_next_vacc_class[dim_I_A_123 * (n_vacc_classes - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_A; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_A_next_vacc_class[dim_I_A_123 * (l - 1 - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  real_t delta_new_conf = odin_sum1<real_t>(n_H_D_unconf_to_conf, 0, dim_H_D_unconf) + odin_sum1<real_t>(n_H_R_unconf_to_conf, 0, dim_H_R_unconf) + odin_sum1<real_t>(n_ICU_pre_unconf_to_conf, 0, dim_ICU_pre_unconf) + odin_sum1<real_t>(n_ICU_D_unconf_to_conf, 0, dim_ICU_D_unconf) + odin_sum1<real_t>(n_ICU_W_R_unconf_to_conf, 0, dim_ICU_W_R_unconf) + odin_sum1<real_t>(n_ICU_W_D_unconf_to_conf, 0, dim_ICU_W_D_unconf) + odin_sum1<real_t>(n_W_R_unconf_to_conf, 0, dim_W_R_unconf) + odin_sum1<real_t>(n_W_D_unconf_to_conf, 0, dim_W_D_unconf);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_EI_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_EE[dim_E_123 * (k - 1) + dim_E_12 * (k_E - 1) + 19 * (j - 1) + i - 1] - n_EI_A[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_H_D[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_hosp_non_ICU[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_H_D_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_SE_next_vacc_class[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_S_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_S_next_vacc_class[19 * (k - 1) + i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      n_S_next_vacc_class[i - 1 + 19 * (j - 1)] = dust::distr::rbinom(rng_state, std::round(S[19 * (j - 1) + i - 1] - odin_sum3<real_t>(n_S_progress, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12)), p_S_next_vacc_class[19 * (j - 1) + i - 1]);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_W_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_W_R_unconf_to_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = aux_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_W_R_unconf_to_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_ICU_pre_unconf_to_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_ICU_pre_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = new_ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_ICU_pre_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_ICU_pre_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = aux_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_ICU_pre_unconf_to_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_ICU_pre_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = new_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_ICU_pre[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_ICU_pre_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_1; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_I_C_1[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = I_C_1[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_I_C_1[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        new_R[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_A[dim_I_A_123 * (k - 1) + dim_E_12 * (k_A - 1) + 19 * (j - 1) + i - 1] + n_I_C_2_to_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_R_conf_progress[dim_H_R_unconf_123 * (k - 1) + dim_E_12 * (k_H_R - 1) + 19 * (j - 1) + i - 1] + n_H_R_unconf_progress[dim_H_R_unconf_123 * (k - 1) + dim_E_12 * (k_H_R - 1) + 19 * (j - 1) + i - 1] + n_W_R_conf_progress[dim_W_R_unconf_123 * (k - 1) + dim_E_12 * (k_W_R - 1) + 19 * (j - 1) + i - 1] + n_W_R_unconf_progress[dim_W_R_unconf_123 * (k - 1) + dim_E_12 * (k_W_R - 1) + 19 * (j - 1) + i - 1] - n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_R_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      new_R[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_R[dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_II_A_next_vacc_class[dim_I_A_123 * (n_vacc_classes - 1) + dim_E_12 * (k_A - 1) + 19 * (j - 1) + i - 1] + n_R_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= n_vacc_classes; ++k) {
        new_R[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_A_next_vacc_class[dim_I_A_123 * (k - 1 - 1) + dim_E_12 * (k_A - 1) + 19 * (j - 1) + i - 1] + n_R_next_vacc_class[dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        new_T_PCR_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = T_PCR_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_T_PCR_pos_progress[dim_T_PCR_pos_123 * (k - 1) + dim_E_12 * (k_PCR_pos - 1) + 19 * (j - 1) + i - 1] - model_pcr_and_serology * n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      new_T_PCR_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_PCR_neg[dim_E_12 * 0 + 19 * (j - 1) + i - 1] + model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= n_vacc_classes; ++k) {
        new_T_PCR_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_PCR_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_PCR_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pre_123 * (l - 1)] = T_PCR_pre[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_T_PCR_pre_progress[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_T_PCR_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pre_123 * (l - 1)] = new_T_PCR_pre[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_S_progress[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_PCR_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_PCR_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pre_123 * (l - 1)] = new_T_PCR_pre[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_T_PCR_pre_progress[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        new_T_sero_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = T_sero_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + odin_sum4<real_t>(n_T_sero_pre_progress, i - 1, i, j - 1, j, 0, 2, k - 1, k, 19, dim_E_12, dim_T_sero_pre_123) - n_T_sero_pre_to_T_sero_pos[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_T_sero_pos_progress[dim_T_sero_pos_123 * (k - 1) + dim_E_12 * (k_sero_pos - 1) + 19 * (j - 1) + i - 1] - model_pcr_and_serology * n_R_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      new_T_sero_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_sero_neg[dim_E_12 * 0 + 19 * (j - 1) + i - 1] + model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= n_vacc_classes; ++k) {
        new_T_sero_neg[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_sero_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + model_pcr_and_serology * n_R_next_vacc_class[dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= 2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_T_sero_pre[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = T_sero_pre[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_com_to_T_sero_pre[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_T_sero_pre_progress[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_W_D_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = new_ICU_W_D_conf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_W_D_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_D_unconf_123 * (l - 1)] = new_ICU_W_D_unconf[dim_ICU_W_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  state_next[3] = cum_infections + odin_sum1<real_t>(n_S_progress, 0, dim_R);
  for (int i = 1; i <= dim_strain_transmission; ++i) {
    state_next[141 + i - 1] = cum_infections_per_strain[i - 1] + odin_sum3<real_t>(n_S_progress, 0, 19, i - 1, i, 0, dim_rel_susceptibility_2, 19, dim_E_12);
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_R_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_R_vaccinated[19 * (j - 1) + i - 1] + odin_sum3<real_t>(n_R_next_vacc_class, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12) + odin_sum3<real_t>(n_RS_next_vacc_class, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12);
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = n_EI_P[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = n_II_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_II_P_next_vacc_class[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_P_next_vacc_class[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        int l = 1;
        aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_I_P_next_vacc_class[dim_I_P_123 * (n_vacc_classes - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_I_P_next_vacc_class[dim_I_P_123 * (l - 1 - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      int l = 1;
      aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EI_P_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EI_P_next_vacc_class[dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_P; ++k) {
        int l = 1;
        aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_P_next_vacc_class[dim_I_P_123 * (n_vacc_classes - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_P; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_II_P_next_vacc_class[dim_I_P_123 * (l - 1 - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_to_H_D[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_star_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_H_R[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_hosp_non_ICU[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_H_D[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= n_vacc_classes; ++k) {
        n_SE[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = n_S_progress[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_SE_next_vacc_class[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  real_t new_ICU_tot = odin_sum1<real_t>(new_ICU_W_R_conf, 0, dim_ICU_W_R_unconf) + odin_sum1<real_t>(new_ICU_W_D_conf, 0, dim_ICU_W_D_unconf) + odin_sum1<real_t>(new_ICU_D_conf, 0, dim_ICU_D_unconf);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_I_A[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      new_S[i - 1 + 19 * (j - 1)] = S[19 * (j - 1) + i - 1] + odin_sum3<real_t>(n_RS, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12) - odin_sum3<real_t>(n_S_progress, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12) - n_S_next_vacc_class[19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    int j = 1;
    new_S[i - 1 + 19 * (j - 1)] = new_S[19 * 0 + i - 1] + n_S_next_vacc_class[19 * (n_vacc_classes - 1) + i - 1] + odin_sum3<real_t>(n_RS_next_vacc_class, i - 1, i, 0, dim_strain_transmission, n_vacc_classes - 1, n_vacc_classes, 19, dim_E_12);
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 2; j <= n_vacc_classes; ++j) {
      new_S[i - 1 + 19 * (j - 1)] = new_S[19 * (j - 1) + i - 1] + n_S_next_vacc_class[19 * (j - 1 - 1) + i - 1] + odin_sum3<real_t>(n_RS_next_vacc_class, i - 1, i, 0, dim_strain_transmission, j - 1 - 1, j - 1, 19, dim_E_12);
    }
  }
  real_t new_sympt_cases = odin_sum1<real_t>(n_EI_P, 0, dim_R) + odin_sum1<real_t>(n_EI_P_next_vacc_class, 0, dim_R);
  real_t new_sympt_cases_non_variant_over25 = odin_sum3<real_t>(n_EI_P, 5, n_groups, 0, 1, 0, dim_rel_susceptibility_2, 19, dim_E_12) + odin_sum3<real_t>(n_EI_P_next_vacc_class, 5, n_groups, 0, 1, 0, dim_rel_susceptibility_2, 19, dim_E_12);
  real_t new_sympt_cases_over25 = odin_sum3<real_t>(n_EI_P, 5, n_groups, 0, dim_strain_transmission, 0, dim_rel_susceptibility_2, 19, dim_E_12) + odin_sum3<real_t>(n_EI_P_next_vacc_class, 5, n_groups, 0, dim_strain_transmission, 0, dim_rel_susceptibility_2, 19, dim_E_12);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_W_R_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = new_ICU_W_R_conf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_W_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_W_R_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_W_R_unconf_123 * (l - 1)] = new_ICU_W_R_unconf[dim_ICU_W_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_pre_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = new_ICU_pre_conf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_ICU_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_ICU_pre_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_ICU_pre_unconf_123 * (l - 1)] = new_ICU_pre_unconf[dim_ICU_pre_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_C_1; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_I_C_1 + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_C_1_123 * (l - 1)] = new_I_C_1[dim_I_C_1_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        state_next[offset_variable_R + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        state_next[offset_variable_T_PCR_neg + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_PCR_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_PCR_pre; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_T_PCR_pre + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_PCR_pre_123 * (l - 1)] = new_T_PCR_pre[dim_T_PCR_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        state_next[offset_variable_T_sero_neg + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = new_T_sero_neg[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= 2; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_T_sero_pre + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_T_sero_pre_123 * (l - 1)] = new_T_sero_pre[dim_T_sero_pre_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_cum_n_S_vaccinated + i - 1 + 19 * (j - 1)] = cum_n_S_vaccinated[19 * (j - 1) + i - 1] + n_S_next_vacc_class[19 * (j - 1) + i - 1] + odin_sum3<real_t>(n_SE_next_vacc_class, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12);
    }
  }
  state_next[5] = cum_new_conf + delta_new_conf;
  state_next[2] = (fmodr<real_t>(step, steps_per_day) == 0 ? delta_new_conf : new_conf_inc + delta_new_conf);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = n_SE[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = n_EE[dim_E_123 * (l - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_EE[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_EE_next_vacc_class[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_E_next_vacc_class[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        int l = 1;
        aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * 0 + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_E_next_vacc_class[dim_E_123 * (n_vacc_classes - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_E_next_vacc_class[dim_E_123 * (l - 1 - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      int l = 1;
      aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * 0 + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_SE_next_vacc_class[dim_E_12 * (n_vacc_classes - 1) + 19 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 2; l <= n_vacc_classes; ++l) {
        aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_SE_next_vacc_class[dim_E_12 * (l - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_E; ++k) {
        int l = 1;
        aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (n_vacc_classes - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 2; k <= k_E; ++k) {
        for (int l = 2; l <= n_vacc_classes; ++l) {
          aux_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_EE_next_vacc_class[dim_E_123 * (l - 1 - 1) + dim_E_12 * (k - 1 - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        n_I_C_2_to_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = dust::distr::rbinom(rng_state, std::round(n_I_C_2_to_H_R[dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1]), p_star_by_age[i - 1]);
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_D_unconf_to_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_H_D_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = new_H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_H_D_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_H_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = aux_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_D_unconf_to_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_H_D_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = new_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_H_D[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_H_D_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_I_P[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  state_next[9] = new_ICU_tot;
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_A; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_I_A + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_A_123 * (l - 1)] = new_I_A[dim_I_A_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_S + i - 1 + 19 * (j - 1)] = new_S[19 * (j - 1) + i - 1];
    }
  }
  state_next[20] = cum_sympt_cases + new_sympt_cases;
  state_next[22] = cum_sympt_cases_non_variant_over25 + new_sympt_cases_non_variant_over25;
  state_next[21] = cum_sympt_cases_over25 + new_sympt_cases_over25;
  state_next[23] = ((fmodr<real_t>(step, steps_per_day) == 0 ? new_sympt_cases : sympt_cases_inc + new_sympt_cases));
  state_next[25] = ((fmodr<real_t>(step, steps_per_day) == 0 ? new_sympt_cases_non_variant_over25 : sympt_cases_non_variant_over25_inc + new_sympt_cases_non_variant_over25));
  state_next[24] = ((fmodr<real_t>(step, steps_per_day) == 0 ? new_sympt_cases_over25 : sympt_cases_over25_inc + new_sympt_cases_over25));
  real_t delta_admit_conf = odin_sum1<real_t>(n_I_C_2_to_H_D_conf, 0, dim_R) + odin_sum1<real_t>(n_I_C_2_to_H_R_conf, 0, dim_R) + odin_sum1<real_t>(n_I_C_2_to_ICU_pre_conf, 0, dim_R);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_E[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + aux_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] + n_H_R_unconf_to_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_H_R_conf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = new_H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_H_R_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          new_H_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = aux_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1] - n_H_R_unconf_to_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      int k = 1;
      for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
        new_H_R_unconf[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = new_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * 0 + 19 * (j - 1) + i - 1] + n_I_C_2_to_H_R[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1] - n_I_C_2_to_H_R_conf[dim_E_12 * (l - 1) + 19 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_H_D_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = new_H_D_conf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_D; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_H_D_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_D_unconf_123 * (l - 1)] = new_H_D_unconf[dim_H_D_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_P; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_I_P + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_I_P_123 * (l - 1)] = new_I_P[dim_I_P_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= dim_rel_susceptibility_2; ++k) {
        I_weighted_strain[i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1)] = strain_transmission[j - 1] * (I_A_transmission * odin_sum4<real_t>(new_I_A, i - 1, i, j - 1, j, 0, k_A, k - 1, k, 19, dim_E_12, dim_I_A_123) + I_P_transmission * odin_sum4<real_t>(new_I_P, i - 1, i, j - 1, j, 0, k_P, k - 1, k, 19, dim_E_12, dim_I_P_123) + I_C_1_transmission * odin_sum4<real_t>(new_I_C_1, i - 1, i, j - 1, j, 0, k_C_1, k - 1, k, 19, dim_E_12, dim_I_C_1_123) + I_C_2_transmission * odin_sum4<real_t>(new_I_C_2, i - 1, i, j - 1, j, 0, k_C_2, k - 1, k, 19, dim_E_12, dim_I_C_2_123) + hosp_transmission * (odin_sum4<real_t>(new_ICU_pre_unconf, i - 1, i, j - 1, j, 0, k_ICU_pre, k - 1, k, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(new_ICU_pre_conf, i - 1, i, j - 1, j, 0, k_ICU_pre, k - 1, k, 19, dim_E_12, dim_ICU_pre_unconf_123) + odin_sum4<real_t>(new_H_R_unconf, i - 1, i, j - 1, j, 0, k_H_R, k - 1, k, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(new_H_R_conf, i - 1, i, j - 1, j, 0, k_H_R, k - 1, k, 19, dim_E_12, dim_H_R_unconf_123) + odin_sum4<real_t>(new_H_D_unconf, i - 1, i, j - 1, j, 0, k_H_D, k - 1, k, 19, dim_E_12, dim_H_D_unconf_123) + odin_sum4<real_t>(new_H_D_conf, i - 1, i, j - 1, j, 0, k_H_D, k - 1, k, 19, dim_E_12, dim_H_D_unconf_123)) + ICU_transmission * (odin_sum4<real_t>(new_ICU_W_R_unconf, i - 1, i, j - 1, j, 0, k_ICU_W_R, k - 1, k, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(new_ICU_W_R_conf, i - 1, i, j - 1, j, 0, k_ICU_W_R, k - 1, k, 19, dim_E_12, dim_ICU_W_R_unconf_123) + odin_sum4<real_t>(new_ICU_W_D_unconf, i - 1, i, j - 1, j, 0, k_ICU_W_D, k - 1, k, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(new_ICU_W_D_conf, i - 1, i, j - 1, j, 0, k_ICU_W_D, k - 1, k, 19, dim_E_12, dim_ICU_W_D_unconf_123) + odin_sum4<real_t>(new_ICU_D_unconf, i - 1, i, j - 1, j, 0, k_ICU_D, k - 1, k, 19, dim_E_12, dim_ICU_D_unconf_123) + odin_sum4<real_t>(new_ICU_D_conf, i - 1, i, j - 1, j, 0, k_ICU_D, k - 1, k, 19, dim_E_12, dim_ICU_D_unconf_123)) + G_D_transmission * odin_sum4<real_t>(new_G_D, i - 1, i, j - 1, j, 0, k_G_D, k - 1, k, 19, dim_E_12, dim_G_D_123));
      }
    }
  }
  real_t new_general_tot = odin_sum1<real_t>(new_ICU_pre_conf, 0, dim_ICU_pre_unconf) + odin_sum1<real_t>(new_H_R_conf, 0, dim_H_R_unconf) + odin_sum1<real_t>(new_H_D_conf, 0, dim_H_D_unconf) + odin_sum1<real_t>(new_W_R_conf, 0, dim_W_R_unconf) + odin_sum1<real_t>(new_W_D_conf, 0, dim_W_D_unconf);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_E; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_E + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_E_123 * (l - 1)] = new_E[dim_E_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_H_R_conf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = new_H_R_conf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_strain_transmission; ++j) {
      for (int k = 1; k <= k_H_R; ++k) {
        for (int l = 1; l <= dim_rel_susceptibility_2; ++l) {
          state_next[offset_variable_H_R_unconf + i - 1 + 19 * (j - 1) + dim_E_12 * (k - 1) + dim_H_R_unconf_123 * (l - 1)] = new_H_R_unconf[dim_H_R_unconf_123 * (l - 1) + dim_E_12 * (k - 1) + 19 * (j - 1) + i - 1];
        }
      }
    }
  }
  state_next[1] = (fmodr<real_t>(step, steps_per_day) == 0 ? delta_admit_conf : admit_conf_inc + delta_admit_conf);
  state_next[4] = cum_admit_conf + delta_admit_conf;
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= dim_rel_susceptibility_2; ++j) {
      state_next[offset_variable_I_weighted + i - 1 + 19 * (j - 1)] = odin_sum3<real_t>(I_weighted_strain, i - 1, i, 0, dim_strain_transmission, j - 1, j, 19, dim_E_12);
    }
  }
  state_next[10] = new_general_tot;
  state_next[11] = new_ICU_tot + new_general_tot;
}
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1) {
  real_t tot = 0.0;
  for (int j = from_j; j < to_j; ++j) {
    int jj = j * dim_x_1;
    for (int i = from_i; i < to_i; ++i) {
      tot += x[i + jj];
    }
  }
  return tot;
}
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum3(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int dim_x_1, int dim_x_12) {
  real_t tot = 0.0;
  for (int k = from_k; k < to_k; ++k) {
    int kk = k * dim_x_12;
    for (int j = from_j; j < to_j; ++j) {
      int jj = j * dim_x_1 + kk;
      for (int i = from_i; i < to_i; ++i) {
        tot += x[i + jj];
      }
    }
  }
  return tot;
}
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum4(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int from_l, int to_l, int dim_x_1, int dim_x_12, int dim_x_123) {
  real_t tot = 0.0;
  for (int l = from_l; l < to_l; ++l) {
    int ll = l * dim_x_123;
    for (int k = from_k; k < to_k; ++k) {
      int kk = k * dim_x_12 + ll;
      for (int j = from_j; j < to_j; ++j) {
        int jj = j * dim_x_1 + kk;
        for (int i = from_i; i < to_i; ++i) {
          tot += x[i + jj];
        }
      }
    }
  }
  return tot;
}
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_t, typename container>
HOSTDEVICE real_t odin_sum1(const container x, size_t from, size_t to) {
  real_t tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
template<>
dust::pars_t<carehomes> dust_pars<carehomes>(cpp11::list user) {
  typedef typename carehomes::real_t real_t;
  auto shared = std::make_shared<carehomes::shared_t>();
  carehomes::internal_t internal;
  internal.delta_D_hosp = std::vector<real_t>(19);
  internal.delta_D_non_hosp = std::vector<real_t>(19);
  shared->gamma_sero_pre = std::vector<real_t>(2);
  shared->initial_D_hosp = std::vector<real_t>(19);
  shared->initial_D_non_hosp = std::vector<real_t>(19);
  shared->initial_N_tot = std::vector<real_t>(19);
  shared->initial_cum_admit_by_age = std::vector<real_t>(19);
  shared->initial_tmp_vaccine_n_candidates = std::vector<real_t>(38);
  internal.p_G_D_by_age = std::vector<real_t>(19);
  internal.p_H_D_by_age = std::vector<real_t>(19);
  internal.p_H_by_age = std::vector<real_t>(19);
  internal.p_ICU_D_by_age = std::vector<real_t>(19);
  internal.p_ICU_by_age = std::vector<real_t>(19);
  shared->p_RS = std::vector<real_t>(19);
  internal.p_W_D_by_age = std::vector<real_t>(19);
  internal.p_star_by_age = std::vector<real_t>(19);
  internal.vaccine_n_candidates = std::vector<real_t>(38);
  internal.vaccine_probability_doses = std::vector<real_t>(38);
  shared->initial_D_carehomes_inc = 0;
  shared->initial_D_carehomes_tot = 0;
  shared->initial_D_comm_inc = 0;
  shared->initial_D_comm_tot = 0;
  for (int i = 1; i <= 19; ++i) {
    shared->initial_D_hosp[i - 1] = 0;
  }
  shared->initial_D_hosp_inc = 0;
  shared->initial_D_hosp_tot = 0;
  for (int i = 1; i <= 19; ++i) {
    shared->initial_D_non_hosp[i - 1] = 0;
  }
  shared->initial_D_tot = 0;
  shared->initial_ICU_tot = 0;
  for (int i = 1; i <= 19; ++i) {
    shared->initial_N_tot[i - 1] = 0;
  }
  shared->initial_N_tot2 = 0;
  shared->initial_N_tot3 = 0;
  shared->initial_admit_conf_inc = 0;
  for (int i = 1; i <= 19; ++i) {
    shared->initial_cum_admit_by_age[i - 1] = 0;
  }
  shared->initial_cum_admit_conf = 0;
  shared->initial_cum_infections = 0;
  shared->initial_cum_new_conf = 0;
  shared->initial_cum_sympt_cases = 0;
  shared->initial_cum_sympt_cases_non_variant_over25 = 0;
  shared->initial_cum_sympt_cases_over25 = 0;
  shared->initial_general_tot = 0;
  shared->initial_hosp_tot = 0;
  shared->initial_new_conf_inc = 0;
  shared->initial_react_pos = 0;
  shared->initial_sero_pos = 0;
  shared->initial_sympt_cases_inc = 0;
  shared->initial_sympt_cases_non_variant_over25_inc = 0;
  shared->initial_sympt_cases_over25_inc = 0;
  shared->initial_time = 0;
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= 2; ++j) {
      shared->initial_tmp_vaccine_n_candidates[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  shared->n_age_groups = 17;
  shared->n_doses = 2;
  shared->n_groups = 19;
  shared->G_D_transmission = NA_REAL;
  shared->ICU_transmission = NA_REAL;
  shared->I_A_transmission = NA_REAL;
  shared->I_C_1_transmission = NA_REAL;
  shared->I_C_2_transmission = NA_REAL;
  shared->I_P_transmission = NA_REAL;
  shared->N_tot_15_64 = NA_REAL;
  shared->N_tot_all = NA_REAL;
  shared->N_tot_over25 = NA_REAL;
  shared->N_tot_react = NA_REAL;
  shared->exp_noise = NA_REAL;
  shared->hosp_transmission = NA_REAL;
  shared->k_A = NA_INTEGER;
  shared->k_C_1 = NA_INTEGER;
  shared->k_C_2 = NA_INTEGER;
  shared->k_E = NA_INTEGER;
  shared->k_G_D = NA_INTEGER;
  shared->k_H_D = NA_INTEGER;
  shared->k_H_R = NA_INTEGER;
  shared->k_ICU_D = NA_INTEGER;
  shared->k_ICU_W_D = NA_INTEGER;
  shared->k_ICU_W_R = NA_INTEGER;
  shared->k_ICU_pre = NA_INTEGER;
  shared->k_P = NA_INTEGER;
  shared->k_PCR_pos = NA_INTEGER;
  shared->k_PCR_pre = NA_INTEGER;
  shared->k_W_D = NA_INTEGER;
  shared->k_W_R = NA_INTEGER;
  shared->k_sero_pos = NA_INTEGER;
  shared->kappa_ICU = NA_REAL;
  shared->kappa_admitted = NA_REAL;
  shared->kappa_all_admission = NA_REAL;
  shared->kappa_death = NA_REAL;
  shared->kappa_death_carehomes = NA_REAL;
  shared->kappa_death_comm = NA_REAL;
  shared->kappa_death_hosp = NA_REAL;
  shared->kappa_death_non_hosp = NA_REAL;
  shared->kappa_diagnoses = NA_REAL;
  shared->kappa_general = NA_REAL;
  shared->kappa_hosp = NA_REAL;
  shared->kappa_pillar2_cases = NA_REAL;
  shared->p_NC = NA_REAL;
  shared->phi_ICU = NA_REAL;
  shared->phi_admitted = NA_REAL;
  shared->phi_all_admission = NA_REAL;
  shared->phi_death_carehomes = NA_REAL;
  shared->phi_death_comm = NA_REAL;
  shared->phi_death_hosp = NA_REAL;
  shared->phi_diagnoses = NA_REAL;
  shared->phi_general = NA_REAL;
  shared->phi_hosp = NA_REAL;
  shared->phi_pillar2_cases = NA_REAL;
  shared->pillar2_sensitivity = NA_REAL;
  shared->pillar2_specificity = NA_REAL;
  shared->react_sensitivity = NA_REAL;
  shared->react_specificity = NA_REAL;
  shared->rho_pillar2_tests = NA_REAL;
  shared->sero_sensitivity = NA_REAL;
  shared->sero_specificity = NA_REAL;
  shared->steps_per_day = NA_INTEGER;
  shared->gamma_A = 0.10000000000000001;
  shared->gamma_C_1 = 0.10000000000000001;
  shared->gamma_C_2 = 0.10000000000000001;
  shared->gamma_E = 0.10000000000000001;
  shared->gamma_G_D = 0.10000000000000001;
  shared->gamma_H_D = 0.10000000000000001;
  shared->gamma_H_R = 0.10000000000000001;
  shared->gamma_ICU_D = 0.10000000000000001;
  shared->gamma_ICU_W_D = 0.10000000000000001;
  shared->gamma_ICU_W_R = 0.10000000000000001;
  shared->gamma_ICU_pre = 0.10000000000000001;
  shared->gamma_P = 0.10000000000000001;
  shared->gamma_PCR_pos = 0.10000000000000001;
  shared->gamma_PCR_pre = 0.10000000000000001;
  shared->gamma_U = 0.10000000000000001;
  shared->gamma_W_D = 0.10000000000000001;
  shared->gamma_W_R = 0.10000000000000001;
  shared->gamma_sero_pos = 0.10000000000000001;
  shared->gamma_sero_pre_1 = 0.10000000000000001;
  shared->gamma_sero_pre_2 = 0.10000000000000001;
  shared->model_pcr_and_serology_user = 1;
  shared->p_sero_pre_1 = 0.5;
  shared->G_D_transmission = user_get_scalar<real_t>(user, "G_D_transmission", shared->G_D_transmission, NA_REAL, NA_REAL);
  shared->ICU_transmission = user_get_scalar<real_t>(user, "ICU_transmission", shared->ICU_transmission, NA_REAL, NA_REAL);
  shared->I_A_transmission = user_get_scalar<real_t>(user, "I_A_transmission", shared->I_A_transmission, NA_REAL, NA_REAL);
  shared->I_C_1_transmission = user_get_scalar<real_t>(user, "I_C_1_transmission", shared->I_C_1_transmission, NA_REAL, NA_REAL);
  shared->I_C_2_transmission = user_get_scalar<real_t>(user, "I_C_2_transmission", shared->I_C_2_transmission, NA_REAL, NA_REAL);
  shared->I_P_transmission = user_get_scalar<real_t>(user, "I_P_transmission", shared->I_P_transmission, NA_REAL, NA_REAL);
  shared->N_tot_15_64 = user_get_scalar<real_t>(user, "N_tot_15_64", shared->N_tot_15_64, NA_REAL, NA_REAL);
  shared->N_tot_all = user_get_scalar<real_t>(user, "N_tot_all", shared->N_tot_all, NA_REAL, NA_REAL);
  shared->N_tot_over25 = user_get_scalar<real_t>(user, "N_tot_over25", shared->N_tot_over25, NA_REAL, NA_REAL);
  shared->N_tot_react = user_get_scalar<real_t>(user, "N_tot_react", shared->N_tot_react, NA_REAL, NA_REAL);
  std::array <int, 1> dim_beta_step;
  shared->beta_step = user_get_array_variable<real_t, 1>(user, "beta_step", shared->beta_step, dim_beta_step, NA_REAL, NA_REAL);
  shared->dim_beta_step = shared->beta_step.size();
  shared->exp_noise = user_get_scalar<real_t>(user, "exp_noise", shared->exp_noise, NA_REAL, NA_REAL);
  shared->gamma_A = user_get_scalar<real_t>(user, "gamma_A", shared->gamma_A, NA_REAL, NA_REAL);
  shared->gamma_C_1 = user_get_scalar<real_t>(user, "gamma_C_1", shared->gamma_C_1, NA_REAL, NA_REAL);
  shared->gamma_C_2 = user_get_scalar<real_t>(user, "gamma_C_2", shared->gamma_C_2, NA_REAL, NA_REAL);
  shared->gamma_E = user_get_scalar<real_t>(user, "gamma_E", shared->gamma_E, NA_REAL, NA_REAL);
  shared->gamma_G_D = user_get_scalar<real_t>(user, "gamma_G_D", shared->gamma_G_D, NA_REAL, NA_REAL);
  shared->gamma_H_D = user_get_scalar<real_t>(user, "gamma_H_D", shared->gamma_H_D, NA_REAL, NA_REAL);
  shared->gamma_H_R = user_get_scalar<real_t>(user, "gamma_H_R", shared->gamma_H_R, NA_REAL, NA_REAL);
  shared->gamma_ICU_D = user_get_scalar<real_t>(user, "gamma_ICU_D", shared->gamma_ICU_D, NA_REAL, NA_REAL);
  shared->gamma_ICU_W_D = user_get_scalar<real_t>(user, "gamma_ICU_W_D", shared->gamma_ICU_W_D, NA_REAL, NA_REAL);
  shared->gamma_ICU_W_R = user_get_scalar<real_t>(user, "gamma_ICU_W_R", shared->gamma_ICU_W_R, NA_REAL, NA_REAL);
  shared->gamma_ICU_pre = user_get_scalar<real_t>(user, "gamma_ICU_pre", shared->gamma_ICU_pre, NA_REAL, NA_REAL);
  shared->gamma_P = user_get_scalar<real_t>(user, "gamma_P", shared->gamma_P, NA_REAL, NA_REAL);
  shared->gamma_PCR_pos = user_get_scalar<real_t>(user, "gamma_PCR_pos", shared->gamma_PCR_pos, NA_REAL, NA_REAL);
  shared->gamma_PCR_pre = user_get_scalar<real_t>(user, "gamma_PCR_pre", shared->gamma_PCR_pre, NA_REAL, NA_REAL);
  shared->gamma_U = user_get_scalar<real_t>(user, "gamma_U", shared->gamma_U, NA_REAL, NA_REAL);
  shared->gamma_W_D = user_get_scalar<real_t>(user, "gamma_W_D", shared->gamma_W_D, NA_REAL, NA_REAL);
  shared->gamma_W_R = user_get_scalar<real_t>(user, "gamma_W_R", shared->gamma_W_R, NA_REAL, NA_REAL);
  shared->gamma_sero_pos = user_get_scalar<real_t>(user, "gamma_sero_pos", shared->gamma_sero_pos, NA_REAL, NA_REAL);
  shared->gamma_sero_pre_1 = user_get_scalar<real_t>(user, "gamma_sero_pre_1", shared->gamma_sero_pre_1, NA_REAL, NA_REAL);
  shared->gamma_sero_pre_2 = user_get_scalar<real_t>(user, "gamma_sero_pre_2", shared->gamma_sero_pre_2, NA_REAL, NA_REAL);
  shared->hosp_transmission = user_get_scalar<real_t>(user, "hosp_transmission", shared->hosp_transmission, NA_REAL, NA_REAL);
  shared->index_dose = user_get_array_fixed<int, 1>(user, "index_dose", shared->index_dose, {2}, NA_REAL, NA_REAL);
  shared->k_A = user_get_scalar<int>(user, "k_A", shared->k_A, NA_REAL, NA_REAL);
  shared->k_C_1 = user_get_scalar<int>(user, "k_C_1", shared->k_C_1, NA_REAL, NA_REAL);
  shared->k_C_2 = user_get_scalar<int>(user, "k_C_2", shared->k_C_2, NA_REAL, NA_REAL);
  shared->k_E = user_get_scalar<int>(user, "k_E", shared->k_E, NA_REAL, NA_REAL);
  shared->k_G_D = user_get_scalar<int>(user, "k_G_D", shared->k_G_D, NA_REAL, NA_REAL);
  shared->k_H_D = user_get_scalar<int>(user, "k_H_D", shared->k_H_D, NA_REAL, NA_REAL);
  shared->k_H_R = user_get_scalar<int>(user, "k_H_R", shared->k_H_R, NA_REAL, NA_REAL);
  shared->k_ICU_D = user_get_scalar<int>(user, "k_ICU_D", shared->k_ICU_D, NA_REAL, NA_REAL);
  shared->k_ICU_W_D = user_get_scalar<int>(user, "k_ICU_W_D", shared->k_ICU_W_D, NA_REAL, NA_REAL);
  shared->k_ICU_W_R = user_get_scalar<int>(user, "k_ICU_W_R", shared->k_ICU_W_R, NA_REAL, NA_REAL);
  shared->k_ICU_pre = user_get_scalar<int>(user, "k_ICU_pre", shared->k_ICU_pre, NA_REAL, NA_REAL);
  shared->k_P = user_get_scalar<int>(user, "k_P", shared->k_P, NA_REAL, NA_REAL);
  shared->k_PCR_pos = user_get_scalar<int>(user, "k_PCR_pos", shared->k_PCR_pos, NA_REAL, NA_REAL);
  shared->k_PCR_pre = user_get_scalar<int>(user, "k_PCR_pre", shared->k_PCR_pre, NA_REAL, NA_REAL);
  shared->k_W_D = user_get_scalar<int>(user, "k_W_D", shared->k_W_D, NA_REAL, NA_REAL);
  shared->k_W_R = user_get_scalar<int>(user, "k_W_R", shared->k_W_R, NA_REAL, NA_REAL);
  shared->k_sero_pos = user_get_scalar<int>(user, "k_sero_pos", shared->k_sero_pos, NA_REAL, NA_REAL);
  shared->kappa_ICU = user_get_scalar<real_t>(user, "kappa_ICU", shared->kappa_ICU, NA_REAL, NA_REAL);
  shared->kappa_admitted = user_get_scalar<real_t>(user, "kappa_admitted", shared->kappa_admitted, NA_REAL, NA_REAL);
  shared->kappa_all_admission = user_get_scalar<real_t>(user, "kappa_all_admission", shared->kappa_all_admission, NA_REAL, NA_REAL);
  shared->kappa_death = user_get_scalar<real_t>(user, "kappa_death", shared->kappa_death, NA_REAL, NA_REAL);
  shared->kappa_death_carehomes = user_get_scalar<real_t>(user, "kappa_death_carehomes", shared->kappa_death_carehomes, NA_REAL, NA_REAL);
  shared->kappa_death_comm = user_get_scalar<real_t>(user, "kappa_death_comm", shared->kappa_death_comm, NA_REAL, NA_REAL);
  shared->kappa_death_hosp = user_get_scalar<real_t>(user, "kappa_death_hosp", shared->kappa_death_hosp, NA_REAL, NA_REAL);
  shared->kappa_death_non_hosp = user_get_scalar<real_t>(user, "kappa_death_non_hosp", shared->kappa_death_non_hosp, NA_REAL, NA_REAL);
  shared->kappa_diagnoses = user_get_scalar<real_t>(user, "kappa_diagnoses", shared->kappa_diagnoses, NA_REAL, NA_REAL);
  shared->kappa_general = user_get_scalar<real_t>(user, "kappa_general", shared->kappa_general, NA_REAL, NA_REAL);
  shared->kappa_hosp = user_get_scalar<real_t>(user, "kappa_hosp", shared->kappa_hosp, NA_REAL, NA_REAL);
  shared->kappa_pillar2_cases = user_get_scalar<real_t>(user, "kappa_pillar2_cases", shared->kappa_pillar2_cases, NA_REAL, NA_REAL);
  shared->m = user_get_array_fixed<real_t, 2>(user, "m", shared->m, {19, 19}, NA_REAL, NA_REAL);
  shared->model_pcr_and_serology_user = user_get_scalar<real_t>(user, "model_pcr_and_serology_user", shared->model_pcr_and_serology_user, NA_REAL, NA_REAL);
  shared->p_C = user_get_array_fixed<real_t, 1>(user, "p_C", shared->p_C, {19}, NA_REAL, NA_REAL);
  std::array <int, 1> dim_p_G_D_step;
  shared->p_G_D_step = user_get_array_variable<real_t, 1>(user, "p_G_D_step", shared->p_G_D_step, dim_p_G_D_step, NA_REAL, NA_REAL);
  shared->dim_p_G_D_step = shared->p_G_D_step.size();
  std::array <int, 1> dim_p_H_D_step;
  shared->p_H_D_step = user_get_array_variable<real_t, 1>(user, "p_H_D_step", shared->p_H_D_step, dim_p_H_D_step, NA_REAL, NA_REAL);
  shared->dim_p_H_D_step = shared->p_H_D_step.size();
  std::array <int, 1> dim_p_H_step;
  shared->p_H_step = user_get_array_variable<real_t, 1>(user, "p_H_step", shared->p_H_step, dim_p_H_step, NA_REAL, NA_REAL);
  shared->dim_p_H_step = shared->p_H_step.size();
  std::array <int, 1> dim_p_ICU_D_step;
  shared->p_ICU_D_step = user_get_array_variable<real_t, 1>(user, "p_ICU_D_step", shared->p_ICU_D_step, dim_p_ICU_D_step, NA_REAL, NA_REAL);
  shared->dim_p_ICU_D_step = shared->p_ICU_D_step.size();
  std::array <int, 1> dim_p_ICU_step;
  shared->p_ICU_step = user_get_array_variable<real_t, 1>(user, "p_ICU_step", shared->p_ICU_step, dim_p_ICU_step, NA_REAL, NA_REAL);
  shared->dim_p_ICU_step = shared->p_ICU_step.size();
  shared->p_NC = user_get_scalar<real_t>(user, "p_NC", shared->p_NC, NA_REAL, NA_REAL);
  std::array <int, 1> dim_p_W_D_step;
  shared->p_W_D_step = user_get_array_variable<real_t, 1>(user, "p_W_D_step", shared->p_W_D_step, dim_p_W_D_step, NA_REAL, NA_REAL);
  shared->dim_p_W_D_step = shared->p_W_D_step.size();
  shared->p_sero_pos = user_get_array_fixed<real_t, 1>(user, "p_sero_pos", shared->p_sero_pos, {19}, NA_REAL, NA_REAL);
  shared->p_sero_pre_1 = user_get_scalar<real_t>(user, "p_sero_pre_1", shared->p_sero_pre_1, NA_REAL, NA_REAL);
  std::array <int, 1> dim_p_star_step;
  shared->p_star_step = user_get_array_variable<real_t, 1>(user, "p_star_step", shared->p_star_step, dim_p_star_step, NA_REAL, NA_REAL);
  shared->dim_p_star_step = shared->p_star_step.size();
  shared->phi_ICU = user_get_scalar<real_t>(user, "phi_ICU", shared->phi_ICU, NA_REAL, NA_REAL);
  shared->phi_admitted = user_get_scalar<real_t>(user, "phi_admitted", shared->phi_admitted, NA_REAL, NA_REAL);
  shared->phi_all_admission = user_get_scalar<real_t>(user, "phi_all_admission", shared->phi_all_admission, NA_REAL, NA_REAL);
  shared->phi_death_carehomes = user_get_scalar<real_t>(user, "phi_death_carehomes", shared->phi_death_carehomes, NA_REAL, NA_REAL);
  shared->phi_death_comm = user_get_scalar<real_t>(user, "phi_death_comm", shared->phi_death_comm, NA_REAL, NA_REAL);
  shared->phi_death_hosp = user_get_scalar<real_t>(user, "phi_death_hosp", shared->phi_death_hosp, NA_REAL, NA_REAL);
  shared->phi_diagnoses = user_get_scalar<real_t>(user, "phi_diagnoses", shared->phi_diagnoses, NA_REAL, NA_REAL);
  shared->phi_general = user_get_scalar<real_t>(user, "phi_general", shared->phi_general, NA_REAL, NA_REAL);
  shared->phi_hosp = user_get_scalar<real_t>(user, "phi_hosp", shared->phi_hosp, NA_REAL, NA_REAL);
  shared->phi_pillar2_cases = user_get_scalar<real_t>(user, "phi_pillar2_cases", shared->phi_pillar2_cases, NA_REAL, NA_REAL);
  shared->pillar2_sensitivity = user_get_scalar<real_t>(user, "pillar2_sensitivity", shared->pillar2_sensitivity, NA_REAL, NA_REAL);
  shared->pillar2_specificity = user_get_scalar<real_t>(user, "pillar2_specificity", shared->pillar2_specificity, NA_REAL, NA_REAL);
  shared->psi_G_D = user_get_array_fixed<real_t, 1>(user, "psi_G_D", shared->psi_G_D, {19}, NA_REAL, NA_REAL);
  shared->psi_H = user_get_array_fixed<real_t, 1>(user, "psi_H", shared->psi_H, {19}, NA_REAL, NA_REAL);
  shared->psi_H_D = user_get_array_fixed<real_t, 1>(user, "psi_H_D", shared->psi_H_D, {19}, NA_REAL, NA_REAL);
  shared->psi_ICU = user_get_array_fixed<real_t, 1>(user, "psi_ICU", shared->psi_ICU, {19}, NA_REAL, NA_REAL);
  shared->psi_ICU_D = user_get_array_fixed<real_t, 1>(user, "psi_ICU_D", shared->psi_ICU_D, {19}, NA_REAL, NA_REAL);
  shared->psi_W_D = user_get_array_fixed<real_t, 1>(user, "psi_W_D", shared->psi_W_D, {19}, NA_REAL, NA_REAL);
  shared->psi_star = user_get_array_fixed<real_t, 1>(user, "psi_star", shared->psi_star, {19}, NA_REAL, NA_REAL);
  shared->react_sensitivity = user_get_scalar<real_t>(user, "react_sensitivity", shared->react_sensitivity, NA_REAL, NA_REAL);
  shared->react_specificity = user_get_scalar<real_t>(user, "react_specificity", shared->react_specificity, NA_REAL, NA_REAL);
  std::array <int, 2> dim_rel_susceptibility;
  shared->rel_susceptibility = user_get_array_variable<real_t, 2>(user, "rel_susceptibility", shared->rel_susceptibility, dim_rel_susceptibility, NA_REAL, NA_REAL);
  shared->dim_rel_susceptibility = shared->rel_susceptibility.size();
  shared->dim_rel_susceptibility_1 = dim_rel_susceptibility[0];
  shared->dim_rel_susceptibility_2 = dim_rel_susceptibility[1];
  shared->rho_pillar2_tests = user_get_scalar<real_t>(user, "rho_pillar2_tests", shared->rho_pillar2_tests, NA_REAL, NA_REAL);
  shared->sero_sensitivity = user_get_scalar<real_t>(user, "sero_sensitivity", shared->sero_sensitivity, NA_REAL, NA_REAL);
  shared->sero_specificity = user_get_scalar<real_t>(user, "sero_specificity", shared->sero_specificity, NA_REAL, NA_REAL);
  shared->steps_per_day = user_get_scalar<int>(user, "steps_per_day", shared->steps_per_day, NA_REAL, NA_REAL);
  std::array <int, 1> dim_strain_seed_step;
  shared->strain_seed_step = user_get_array_variable<real_t, 1>(user, "strain_seed_step", shared->strain_seed_step, dim_strain_seed_step, NA_REAL, NA_REAL);
  shared->dim_strain_seed_step = shared->strain_seed_step.size();
  std::array <int, 1> dim_strain_transmission;
  shared->strain_transmission = user_get_array_variable<real_t, 1>(user, "strain_transmission", shared->strain_transmission, dim_strain_transmission, NA_REAL, NA_REAL);
  shared->dim_strain_transmission = shared->strain_transmission.size();
  std::array <int, 3> dim_vaccine_dose_step;
  shared->vaccine_dose_step = user_get_array_variable<real_t, 3>(user, "vaccine_dose_step", shared->vaccine_dose_step, dim_vaccine_dose_step, NA_REAL, NA_REAL);
  shared->dim_vaccine_dose_step = shared->vaccine_dose_step.size();
  shared->dim_vaccine_dose_step_1 = dim_vaccine_dose_step[0];
  shared->dim_vaccine_dose_step_2 = dim_vaccine_dose_step[1];
  shared->dim_vaccine_dose_step_3 = dim_vaccine_dose_step[2];
  shared->waning_rate = user_get_array_fixed<real_t, 1>(user, "waning_rate", shared->waning_rate, {19}, NA_REAL, NA_REAL);
  shared->dt = 1 / (real_t) shared->steps_per_day;
  {
     int i = 1;
     shared->gamma_sero_pre[i - 1] = shared->gamma_sero_pre_1;
  }
  {
     int i = 2;
     shared->gamma_sero_pre[i - 1] = shared->gamma_sero_pre_2;
  }
  shared->initial_beta_out = shared->beta_step[0];
  shared->model_pcr_and_serology = (shared->model_pcr_and_serology_user == 1 ? 1 : 0);
  shared->initial_cum_infections_per_strain = std::vector<real_t>(shared->dim_strain_transmission);
  shared->dim_E = 19 * shared->dim_strain_transmission * shared->k_E * shared->dim_rel_susceptibility_2;
  shared->dim_E_12 = 19 * shared->dim_strain_transmission;
  shared->dim_E_123 = 19 * shared->dim_strain_transmission * shared->k_E;
  shared->dim_G_D = 19 * shared->dim_strain_transmission * shared->k_G_D * shared->dim_rel_susceptibility_2;
  shared->dim_G_D_123 = 19 * shared->dim_strain_transmission * shared->k_G_D;
  shared->dim_H_D_unconf = 19 * shared->dim_strain_transmission * shared->k_H_D * shared->dim_rel_susceptibility_2;
  shared->dim_H_D_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_H_D;
  shared->dim_H_R_unconf = 19 * shared->dim_strain_transmission * shared->k_H_R * shared->dim_rel_susceptibility_2;
  shared->dim_H_R_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_H_R;
  shared->dim_ICU_D_unconf = 19 * shared->dim_strain_transmission * shared->k_ICU_D * shared->dim_rel_susceptibility_2;
  shared->dim_ICU_D_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_ICU_D;
  shared->dim_ICU_W_D_unconf = 19 * shared->dim_strain_transmission * shared->k_ICU_W_D * shared->dim_rel_susceptibility_2;
  shared->dim_ICU_W_D_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_ICU_W_D;
  shared->dim_ICU_W_R_unconf = 19 * shared->dim_strain_transmission * shared->k_ICU_W_R * shared->dim_rel_susceptibility_2;
  shared->dim_ICU_W_R_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_ICU_W_R;
  shared->dim_ICU_pre_unconf = 19 * shared->dim_strain_transmission * shared->k_ICU_pre * shared->dim_rel_susceptibility_2;
  shared->dim_ICU_pre_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_ICU_pre;
  shared->dim_I_A = 19 * shared->dim_strain_transmission * shared->k_A * shared->dim_rel_susceptibility_2;
  shared->dim_I_A_123 = 19 * shared->dim_strain_transmission * shared->k_A;
  shared->dim_I_C_1 = 19 * shared->dim_strain_transmission * shared->k_C_1 * shared->dim_rel_susceptibility_2;
  shared->dim_I_C_1_123 = 19 * shared->dim_strain_transmission * shared->k_C_1;
  shared->dim_I_C_2 = 19 * shared->dim_strain_transmission * shared->k_C_2 * shared->dim_rel_susceptibility_2;
  shared->dim_I_C_2_123 = 19 * shared->dim_strain_transmission * shared->k_C_2;
  shared->dim_I_P = 19 * shared->dim_strain_transmission * shared->k_P * shared->dim_rel_susceptibility_2;
  shared->dim_I_P_123 = 19 * shared->dim_strain_transmission * shared->k_P;
  shared->dim_R = 19 * shared->dim_strain_transmission * shared->dim_rel_susceptibility_2;
  shared->dim_T_PCR_pos = 19 * shared->dim_strain_transmission * shared->k_PCR_pos * shared->dim_rel_susceptibility_2;
  shared->dim_T_PCR_pos_123 = 19 * shared->dim_strain_transmission * shared->k_PCR_pos;
  shared->dim_T_PCR_pre = 19 * shared->dim_strain_transmission * shared->k_PCR_pre * shared->dim_rel_susceptibility_2;
  shared->dim_T_PCR_pre_123 = 19 * shared->dim_strain_transmission * shared->k_PCR_pre;
  shared->dim_T_sero_pos = 19 * shared->dim_strain_transmission * shared->k_sero_pos * shared->dim_rel_susceptibility_2;
  shared->dim_T_sero_pos_123 = 19 * shared->dim_strain_transmission * shared->k_sero_pos;
  shared->dim_T_sero_pre = 19 * shared->dim_strain_transmission * 2 * shared->dim_rel_susceptibility_2;
  shared->dim_T_sero_pre_123 = 19 * shared->dim_strain_transmission * 2;
  shared->dim_W_D_unconf = 19 * shared->dim_strain_transmission * shared->k_W_D * shared->dim_rel_susceptibility_2;
  shared->dim_W_D_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_W_D;
  shared->dim_W_R_unconf = 19 * shared->dim_strain_transmission * shared->k_W_R * shared->dim_rel_susceptibility_2;
  shared->dim_W_R_unconf_123 = 19 * shared->dim_strain_transmission * shared->k_W_R;
  shared->dim_cum_n_S_vaccinated = 19 * shared->dim_rel_susceptibility_2;
  shared->dim_s_ij = 19 * 19 * shared->dim_strain_transmission;
  shared->dim_vaccine_dose_step_12 = shared->dim_vaccine_dose_step_1 * shared->dim_vaccine_dose_step_2;
  for (int i = 1; i <= shared->dim_strain_transmission; ++i) {
    shared->initial_cum_infections_per_strain[i - 1] = 0;
  }
  shared->n_strains = shared->dim_strain_transmission;
  shared->offset_variable_cum_n_S_vaccinated = shared->dim_strain_transmission + 141;
  shared->p_E_progress = 1 - std::exp(- shared->gamma_E * shared->dt);
  shared->p_G_D_progress = 1 - std::exp(- shared->gamma_G_D * shared->dt);
  shared->p_H_D_progress = 1 - std::exp(- shared->gamma_H_D * shared->dt);
  shared->p_H_R_progress = 1 - std::exp(- shared->gamma_H_R * shared->dt);
  shared->p_ICU_D_progress = 1 - std::exp(- shared->gamma_ICU_D * shared->dt);
  shared->p_ICU_W_D_progress = 1 - std::exp(- shared->gamma_ICU_W_D * shared->dt);
  shared->p_ICU_W_R_progress = 1 - std::exp(- shared->gamma_ICU_W_R * shared->dt);
  shared->p_ICU_pre_progress = 1 - std::exp(- shared->gamma_ICU_pre * shared->dt);
  shared->p_I_A_progress = 1 - std::exp(- shared->gamma_A * shared->dt);
  shared->p_I_C_1_progress = 1 - std::exp(- shared->gamma_C_1 * shared->dt);
  shared->p_I_C_2_progress = 1 - std::exp(- shared->gamma_C_2 * shared->dt);
  shared->p_I_P_progress = 1 - std::exp(- shared->gamma_P * shared->dt);
  for (int i = 1; i <= 19; ++i) {
    shared->p_RS[i - 1] = 1 - std::exp(- shared->waning_rate[i - 1] * shared->dt);
  }
  shared->p_T_PCR_pos_progress = 1 - std::exp(- shared->gamma_PCR_pos * shared->dt);
  shared->p_T_PCR_pre_progress = 1 - std::exp(- shared->gamma_PCR_pre * shared->dt);
  shared->p_T_sero_pos_progress = 1 - std::exp(- shared->gamma_sero_pos * shared->dt);
  shared->p_W_D_progress = 1 - std::exp(- shared->gamma_W_D * shared->dt);
  shared->p_W_R_progress = 1 - std::exp(- shared->gamma_W_R * shared->dt);
  shared->p_test = 1 - std::exp(- shared->gamma_U * shared->dt);
  internal.I_weighted_strain = std::vector<real_t>(shared->dim_R);
  internal.I_with_diff_trans = std::vector<real_t>(shared->dim_R);
  internal.aux_E = std::vector<real_t>(shared->dim_E);
  internal.aux_G_D = std::vector<real_t>(shared->dim_G_D);
  internal.aux_H_D_conf = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.aux_H_D_unconf = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.aux_H_R_conf = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.aux_H_R_unconf = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.aux_ICU_D_conf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.aux_ICU_D_unconf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.aux_ICU_W_D_conf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.aux_ICU_W_D_unconf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.aux_ICU_W_R_conf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.aux_ICU_W_R_unconf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.aux_ICU_pre_conf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.aux_ICU_pre_unconf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.aux_I_A = std::vector<real_t>(shared->dim_I_A);
  internal.aux_I_C_1 = std::vector<real_t>(shared->dim_I_C_1);
  internal.aux_I_C_2 = std::vector<real_t>(shared->dim_I_C_2);
  internal.aux_I_P = std::vector<real_t>(shared->dim_I_P);
  internal.aux_W_D_conf = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.aux_W_D_unconf = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.aux_W_R_conf = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.aux_W_R_unconf = std::vector<real_t>(shared->dim_W_R_unconf);
  shared->initial_E = std::vector<real_t>(shared->dim_E);
  shared->initial_G_D = std::vector<real_t>(shared->dim_G_D);
  shared->initial_H_D_conf = std::vector<real_t>(shared->dim_H_D_unconf);
  shared->initial_H_D_unconf = std::vector<real_t>(shared->dim_H_D_unconf);
  shared->initial_H_R_conf = std::vector<real_t>(shared->dim_H_R_unconf);
  shared->initial_H_R_unconf = std::vector<real_t>(shared->dim_H_R_unconf);
  shared->initial_ICU_D_conf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  shared->initial_ICU_D_unconf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  shared->initial_ICU_W_D_conf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  shared->initial_ICU_W_D_unconf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  shared->initial_ICU_W_R_conf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  shared->initial_ICU_W_R_unconf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  shared->initial_ICU_pre_conf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  shared->initial_ICU_pre_unconf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  shared->initial_I_A = std::vector<real_t>(shared->dim_I_A);
  shared->initial_I_C_1 = std::vector<real_t>(shared->dim_I_C_1);
  shared->initial_I_C_2 = std::vector<real_t>(shared->dim_I_C_2);
  shared->initial_I_P = std::vector<real_t>(shared->dim_I_P);
  shared->initial_I_weighted = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_R = std::vector<real_t>(shared->dim_R);
  shared->initial_S = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_T_PCR_neg = std::vector<real_t>(shared->dim_R);
  shared->initial_T_PCR_pos = std::vector<real_t>(shared->dim_T_PCR_pos);
  shared->initial_T_PCR_pre = std::vector<real_t>(shared->dim_T_PCR_pre);
  shared->initial_T_sero_neg = std::vector<real_t>(shared->dim_R);
  shared->initial_T_sero_pos = std::vector<real_t>(shared->dim_T_sero_pos);
  shared->initial_T_sero_pre = std::vector<real_t>(shared->dim_T_sero_pre);
  shared->initial_W_D_conf = std::vector<real_t>(shared->dim_W_D_unconf);
  shared->initial_W_D_unconf = std::vector<real_t>(shared->dim_W_D_unconf);
  shared->initial_W_R_conf = std::vector<real_t>(shared->dim_W_R_unconf);
  shared->initial_W_R_unconf = std::vector<real_t>(shared->dim_W_R_unconf);
  shared->initial_cum_n_E_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_cum_n_I_A_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_cum_n_I_P_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_cum_n_R_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_cum_n_S_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_cum_n_vaccinated = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->initial_prob_strain = std::vector<real_t>(shared->dim_E_12);
  shared->initial_tmp_vaccine_probability = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  internal.lambda = std::vector<real_t>(shared->dim_E_12);
  internal.n_EE = std::vector<real_t>(shared->dim_E);
  internal.n_EE_next_vacc_class = std::vector<real_t>(shared->dim_E);
  internal.n_EI_A = std::vector<real_t>(shared->dim_R);
  internal.n_EI_A_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.n_EI_P = std::vector<real_t>(shared->dim_R);
  internal.n_EI_P_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.n_E_next_vacc_class = std::vector<real_t>(shared->dim_E);
  internal.n_E_progress = std::vector<real_t>(shared->dim_E);
  internal.n_G_D_progress = std::vector<real_t>(shared->dim_G_D);
  internal.n_H_D_conf_progress = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.n_H_D_unconf_progress = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.n_H_D_unconf_to_conf = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.n_H_R_conf_progress = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.n_H_R_unconf_progress = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.n_H_R_unconf_to_conf = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.n_ICU_D_conf_progress = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.n_ICU_D_unconf_progress = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.n_ICU_D_unconf_to_conf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.n_ICU_W_D_conf_progress = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.n_ICU_W_D_unconf_progress = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.n_ICU_W_D_unconf_to_conf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.n_ICU_W_R_conf_progress = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.n_ICU_W_R_unconf_progress = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.n_ICU_W_R_unconf_to_conf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.n_ICU_pre_conf_progress = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.n_ICU_pre_conf_to_ICU_D_conf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_conf_to_ICU_W_D_conf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_conf_to_ICU_W_R_conf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_unconf_progress = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.n_ICU_pre_unconf_to_ICU_D_unconf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_unconf_to_ICU_W_D_unconf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_unconf_to_ICU_W_R_unconf = std::vector<real_t>(shared->dim_R);
  internal.n_ICU_pre_unconf_to_conf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.n_II_A = std::vector<real_t>(shared->dim_I_A);
  internal.n_II_A_next_vacc_class = std::vector<real_t>(shared->dim_I_A);
  internal.n_II_P = std::vector<real_t>(shared->dim_I_P);
  internal.n_II_P_next_vacc_class = std::vector<real_t>(shared->dim_I_P);
  internal.n_I_A_next_vacc_class = std::vector<real_t>(shared->dim_I_A);
  internal.n_I_A_progress = std::vector<real_t>(shared->dim_I_A);
  internal.n_I_C_1_progress = std::vector<real_t>(shared->dim_I_C_1);
  internal.n_I_C_2_progress = std::vector<real_t>(shared->dim_I_C_2);
  internal.n_I_C_2_to_G_D = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_H_D = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_H_D_conf = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_H_R = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_H_R_conf = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_ICU_pre = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_ICU_pre_conf = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_R = std::vector<real_t>(shared->dim_R);
  internal.n_I_C_2_to_hosp = std::vector<real_t>(shared->dim_R);
  internal.n_I_P_next_vacc_class = std::vector<real_t>(shared->dim_I_P);
  internal.n_I_P_progress = std::vector<real_t>(shared->dim_I_P);
  internal.n_RS = std::vector<real_t>(shared->dim_R);
  internal.n_RS_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.n_R_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.n_R_next_vacc_class_capped = std::vector<real_t>(shared->dim_R);
  internal.n_R_next_vacc_class_tmp = std::vector<real_t>(shared->dim_R);
  internal.n_R_progress = std::vector<real_t>(shared->dim_R);
  internal.n_R_progress_capped = std::vector<real_t>(shared->dim_R);
  internal.n_R_progress_tmp = std::vector<real_t>(shared->dim_R);
  internal.n_SE = std::vector<real_t>(shared->dim_R);
  internal.n_SE_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.n_S_next_vacc_class = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  internal.n_S_progress = std::vector<real_t>(shared->dim_R);
  internal.n_S_progress_tot = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  internal.n_T_PCR_pos_progress = std::vector<real_t>(shared->dim_T_PCR_pos);
  internal.n_T_PCR_pre_progress = std::vector<real_t>(shared->dim_T_PCR_pre);
  internal.n_T_sero_pos_progress = std::vector<real_t>(shared->dim_T_sero_pos);
  internal.n_T_sero_pre_progress = std::vector<real_t>(shared->dim_T_sero_pre);
  internal.n_T_sero_pre_to_T_sero_pos = std::vector<real_t>(shared->dim_R);
  internal.n_W_D_conf_progress = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.n_W_D_unconf_progress = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.n_W_D_unconf_to_conf = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.n_W_R_conf_progress = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.n_W_R_unconf_progress = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.n_W_R_unconf_to_conf = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.n_com_to_T_sero_pre = std::vector<real_t>(shared->dim_T_sero_pre);
  internal.n_hosp_non_ICU = std::vector<real_t>(shared->dim_R);
  internal.new_E = std::vector<real_t>(shared->dim_E);
  internal.new_G_D = std::vector<real_t>(shared->dim_G_D);
  internal.new_H_D_conf = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.new_H_D_unconf = std::vector<real_t>(shared->dim_H_D_unconf);
  internal.new_H_R_conf = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.new_H_R_unconf = std::vector<real_t>(shared->dim_H_R_unconf);
  internal.new_ICU_D_conf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.new_ICU_D_unconf = std::vector<real_t>(shared->dim_ICU_D_unconf);
  internal.new_ICU_W_D_conf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.new_ICU_W_D_unconf = std::vector<real_t>(shared->dim_ICU_W_D_unconf);
  internal.new_ICU_W_R_conf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.new_ICU_W_R_unconf = std::vector<real_t>(shared->dim_ICU_W_R_unconf);
  internal.new_ICU_pre_conf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.new_ICU_pre_unconf = std::vector<real_t>(shared->dim_ICU_pre_unconf);
  internal.new_I_A = std::vector<real_t>(shared->dim_I_A);
  internal.new_I_C_1 = std::vector<real_t>(shared->dim_I_C_1);
  internal.new_I_C_2 = std::vector<real_t>(shared->dim_I_C_2);
  internal.new_I_P = std::vector<real_t>(shared->dim_I_P);
  internal.new_R = std::vector<real_t>(shared->dim_R);
  internal.new_S = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  internal.new_T_PCR_neg = std::vector<real_t>(shared->dim_R);
  internal.new_T_PCR_pos = std::vector<real_t>(shared->dim_T_PCR_pos);
  internal.new_T_PCR_pre = std::vector<real_t>(shared->dim_T_PCR_pre);
  internal.new_T_sero_neg = std::vector<real_t>(shared->dim_R);
  internal.new_T_sero_pos = std::vector<real_t>(shared->dim_T_sero_pos);
  internal.new_T_sero_pre = std::vector<real_t>(shared->dim_T_sero_pre);
  internal.new_W_D_conf = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.new_W_D_unconf = std::vector<real_t>(shared->dim_W_D_unconf);
  internal.new_W_R_conf = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.new_W_R_unconf = std::vector<real_t>(shared->dim_W_R_unconf);
  internal.p_E_next_vacc_class = std::vector<real_t>(shared->dim_E);
  internal.p_I_A_next_vacc_class = std::vector<real_t>(shared->dim_I_A);
  internal.p_I_P_next_vacc_class = std::vector<real_t>(shared->dim_I_P);
  internal.p_R_next_vacc_class = std::vector<real_t>(shared->dim_R);
  internal.p_SE = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  internal.p_S_next_vacc_class = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  shared->p_T_sero_pre_progress = std::vector<real_t>(shared->dim_T_sero_pre);
  internal.s_ij = std::vector<real_t>(shared->dim_s_ij);
  internal.vaccine_probability = std::vector<real_t>(shared->dim_cum_n_S_vaccinated);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_E; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_E[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_E_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_G_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_G_D[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_G_D_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_H_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_H_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_H_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_H_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_H_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_H_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_H_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_H_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_H_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_W_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_W_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_W_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_pre; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_pre_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_ICU_pre; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_ICU_pre_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_ICU_pre_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_A; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_I_A[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_A_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_C_1; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_I_C_1[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_1_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_C_2; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_I_C_2[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_C_2_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_P; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_I_P[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_I_P_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_I_weighted[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
        shared->initial_R[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = 0;
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_S[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
        shared->initial_T_PCR_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = 0;
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_PCR_pos; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_T_PCR_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pos_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_PCR_pre; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_T_PCR_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_PCR_pre_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->dim_rel_susceptibility_2; ++k) {
        shared->initial_T_sero_neg[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1)] = 0;
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_sero_pos; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_T_sero_pos[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pos_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= 2; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_T_sero_pre[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_W_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_W_D_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_W_D; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_W_D_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_D_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_W_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_W_R_conf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= shared->k_W_R; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->initial_W_R_unconf[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_W_R_unconf_123 * (l - 1)] = 0;
        }
      }
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_E_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_I_A_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_I_P_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_R_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_S_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_cum_n_vaccinated[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      shared->initial_prob_strain[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_rel_susceptibility_2; ++j) {
      shared->initial_tmp_vaccine_probability[i - 1 + 19 * (j - 1)] = 0;
    }
  }
  shared->n_vacc_classes = shared->dim_rel_susceptibility_2;
  shared->offset_variable_E = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + 141;
  shared->offset_variable_G_D = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + 141;
  shared->offset_variable_H_D_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + 141;
  shared->offset_variable_H_D_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + 141;
  shared->offset_variable_H_R_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + 141;
  shared->offset_variable_H_R_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + 141;
  shared->offset_variable_ICU_D_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + 141;
  shared->offset_variable_ICU_D_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + 141;
  shared->offset_variable_ICU_W_D_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + 141;
  shared->offset_variable_ICU_W_D_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + 141;
  shared->offset_variable_ICU_W_R_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + 141;
  shared->offset_variable_ICU_W_R_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + 141;
  shared->offset_variable_ICU_pre_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + 141;
  shared->offset_variable_ICU_pre_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + 141;
  shared->offset_variable_I_A = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + 141;
  shared->offset_variable_I_C_1 = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + 141;
  shared->offset_variable_I_C_2 = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + 141;
  shared->offset_variable_I_P = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + 141;
  shared->offset_variable_I_weighted = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + 141;
  shared->offset_variable_R = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + 141;
  shared->offset_variable_S = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_T_PCR_neg = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + 141;
  shared->offset_variable_T_PCR_pos = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_T_sero_pre + shared->dim_T_sero_pos + shared->dim_T_PCR_pre + 141;
  shared->offset_variable_T_PCR_pre = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_T_sero_pre + shared->dim_T_sero_pos + 141;
  shared->offset_variable_T_sero_neg = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_T_sero_pos = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + shared->dim_T_sero_pre + 141;
  shared->offset_variable_T_sero_pre = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + shared->dim_W_D_unconf + 141;
  shared->offset_variable_W_D_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + shared->dim_W_D_unconf + 141;
  shared->offset_variable_W_D_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + shared->dim_W_R_unconf + 141;
  shared->offset_variable_W_R_conf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + shared->dim_W_R_unconf + 141;
  shared->offset_variable_W_R_unconf = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_R + shared->dim_R + shared->dim_R + shared->dim_E + shared->dim_I_A + shared->dim_I_P + shared->dim_I_C_1 + shared->dim_I_C_2 + shared->dim_G_D + shared->dim_ICU_pre_unconf + shared->dim_ICU_pre_unconf + shared->dim_H_R_unconf + shared->dim_H_R_unconf + shared->dim_H_D_unconf + shared->dim_H_D_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_R_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_W_D_unconf + shared->dim_ICU_D_unconf + shared->dim_ICU_D_unconf + 141;
  shared->offset_variable_cum_n_E_vaccinated = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_cum_n_I_A_vaccinated = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_cum_n_I_P_vaccinated = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_cum_n_R_vaccinated = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_cum_n_vaccinated = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_prob_strain = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + 141;
  shared->offset_variable_tmp_vaccine_probability = shared->dim_strain_transmission + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_cum_n_S_vaccinated + shared->dim_E_12 + shared->dim_cum_n_S_vaccinated + 141;
  shared->rel_infectivity = user_get_array_fixed<real_t, 2>(user, "rel_infectivity", shared->rel_infectivity, {19, shared->dim_rel_susceptibility_2}, NA_REAL, NA_REAL);
  shared->rel_p_hosp_if_sympt = user_get_array_fixed<real_t, 2>(user, "rel_p_hosp_if_sympt", shared->rel_p_hosp_if_sympt, {19, shared->dim_rel_susceptibility_2}, NA_REAL, NA_REAL);
  shared->rel_p_sympt = user_get_array_fixed<real_t, 2>(user, "rel_p_sympt", shared->rel_p_sympt, {19, shared->dim_rel_susceptibility_2}, NA_REAL, NA_REAL);
  shared->vaccine_progression_rate_base = user_get_array_fixed<real_t, 2>(user, "vaccine_progression_rate_base", shared->vaccine_progression_rate_base, {19, shared->dim_rel_susceptibility_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= 19; ++i) {
    for (int j = 1; j <= shared->dim_strain_transmission; ++j) {
      for (int k = 1; k <= 2; ++k) {
        for (int l = 1; l <= shared->dim_rel_susceptibility_2; ++l) {
          shared->p_T_sero_pre_progress[i - 1 + 19 * (j - 1) + shared->dim_E_12 * (k - 1) + shared->dim_T_sero_pre_123 * (l - 1)] = 1 - std::exp(- shared->gamma_sero_pre[k - 1] * shared->dt);
        }
      }
    }
  }
  return dust::pars_t<carehomes>(shared, internal);
}
template <>
cpp11::sexp dust_info<carehomes>(const dust::pars_t<carehomes>& pars) {
  const carehomes::internal_t internal = pars.internal;
  const std::shared_ptr<const carehomes::shared_t> shared = pars.shared;
  cpp11::writable::strings nms({"time", "admit_conf_inc", "new_conf_inc", "cum_infections", "cum_admit_conf", "cum_new_conf", "beta_out", "N_tot2", "N_tot3", "ICU_tot", "general_tot", "hosp_tot", "D_hosp_tot", "D_comm_tot", "D_comm_inc", "D_carehomes_tot", "D_carehomes_inc", "D_hosp_inc", "D_tot", "sero_pos", "cum_sympt_cases", "cum_sympt_cases_over25", "cum_sympt_cases_non_variant_over25", "sympt_cases_inc", "sympt_cases_over25_inc", "sympt_cases_non_variant_over25_inc", "react_pos", "D_hosp", "D_non_hosp", "cum_admit_by_age", "N_tot", "tmp_vaccine_n_candidates", "cum_infections_per_strain", "cum_n_S_vaccinated", "cum_n_E_vaccinated", "cum_n_I_A_vaccinated", "cum_n_I_P_vaccinated", "cum_n_R_vaccinated", "cum_n_vaccinated", "S", "prob_strain", "I_weighted", "tmp_vaccine_probability", "T_sero_neg", "R", "T_PCR_neg", "E", "I_A", "I_P", "I_C_1", "I_C_2", "G_D", "ICU_pre_unconf", "ICU_pre_conf", "H_R_unconf", "H_R_conf", "H_D_unconf", "H_D_conf", "ICU_W_R_unconf", "ICU_W_R_conf", "ICU_W_D_unconf", "ICU_W_D_conf", "ICU_D_unconf", "ICU_D_conf", "W_R_unconf", "W_R_conf", "W_D_unconf", "W_D_conf", "T_sero_pre", "T_sero_pos", "T_PCR_pre", "T_PCR_pos"});
  cpp11::writable::list dim(72);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({1});
  dim[6] = cpp11::writable::integers({1});
  dim[7] = cpp11::writable::integers({1});
  dim[8] = cpp11::writable::integers({1});
  dim[9] = cpp11::writable::integers({1});
  dim[10] = cpp11::writable::integers({1});
  dim[11] = cpp11::writable::integers({1});
  dim[12] = cpp11::writable::integers({1});
  dim[13] = cpp11::writable::integers({1});
  dim[14] = cpp11::writable::integers({1});
  dim[15] = cpp11::writable::integers({1});
  dim[16] = cpp11::writable::integers({1});
  dim[17] = cpp11::writable::integers({1});
  dim[18] = cpp11::writable::integers({1});
  dim[19] = cpp11::writable::integers({1});
  dim[20] = cpp11::writable::integers({1});
  dim[21] = cpp11::writable::integers({1});
  dim[22] = cpp11::writable::integers({1});
  dim[23] = cpp11::writable::integers({1});
  dim[24] = cpp11::writable::integers({1});
  dim[25] = cpp11::writable::integers({1});
  dim[26] = cpp11::writable::integers({1});
  dim[27] = cpp11::writable::integers({19});
  dim[28] = cpp11::writable::integers({19});
  dim[29] = cpp11::writable::integers({19});
  dim[30] = cpp11::writable::integers({19});
  dim[31] = cpp11::writable::integers({19, 2});
  dim[32] = cpp11::writable::integers({shared->dim_strain_transmission});
  dim[33] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[34] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[35] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[36] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[37] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[38] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[39] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[40] = cpp11::writable::integers({19, shared->dim_strain_transmission});
  dim[41] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[42] = cpp11::writable::integers({19, shared->dim_rel_susceptibility_2});
  dim[43] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->dim_rel_susceptibility_2});
  dim[44] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->dim_rel_susceptibility_2});
  dim[45] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->dim_rel_susceptibility_2});
  dim[46] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_E, shared->dim_rel_susceptibility_2});
  dim[47] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_A, shared->dim_rel_susceptibility_2});
  dim[48] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_P, shared->dim_rel_susceptibility_2});
  dim[49] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_C_1, shared->dim_rel_susceptibility_2});
  dim[50] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_C_2, shared->dim_rel_susceptibility_2});
  dim[51] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_G_D, shared->dim_rel_susceptibility_2});
  dim[52] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_pre, shared->dim_rel_susceptibility_2});
  dim[53] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_pre, shared->dim_rel_susceptibility_2});
  dim[54] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_H_R, shared->dim_rel_susceptibility_2});
  dim[55] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_H_R, shared->dim_rel_susceptibility_2});
  dim[56] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_H_D, shared->dim_rel_susceptibility_2});
  dim[57] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_H_D, shared->dim_rel_susceptibility_2});
  dim[58] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_W_R, shared->dim_rel_susceptibility_2});
  dim[59] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_W_R, shared->dim_rel_susceptibility_2});
  dim[60] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_W_D, shared->dim_rel_susceptibility_2});
  dim[61] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_W_D, shared->dim_rel_susceptibility_2});
  dim[62] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_D, shared->dim_rel_susceptibility_2});
  dim[63] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_ICU_D, shared->dim_rel_susceptibility_2});
  dim[64] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_W_R, shared->dim_rel_susceptibility_2});
  dim[65] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_W_R, shared->dim_rel_susceptibility_2});
  dim[66] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_W_D, shared->dim_rel_susceptibility_2});
  dim[67] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_W_D, shared->dim_rel_susceptibility_2});
  dim[68] = cpp11::writable::integers({19, shared->dim_strain_transmission, 2, shared->dim_rel_susceptibility_2});
  dim[69] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_sero_pos, shared->dim_rel_susceptibility_2});
  dim[70] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_PCR_pre, shared->dim_rel_susceptibility_2});
  dim[71] = cpp11::writable::integers({19, shared->dim_strain_transmission, shared->k_PCR_pos, shared->dim_rel_susceptibility_2});
  dim.names() = nms;
  cpp11::writable::list index(72);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = cpp11::writable::integers({6});
  index[6] = cpp11::writable::integers({7});
  index[7] = cpp11::writable::integers({8});
  index[8] = cpp11::writable::integers({9});
  index[9] = cpp11::writable::integers({10});
  index[10] = cpp11::writable::integers({11});
  index[11] = cpp11::writable::integers({12});
  index[12] = cpp11::writable::integers({13});
  index[13] = cpp11::writable::integers({14});
  index[14] = cpp11::writable::integers({15});
  index[15] = cpp11::writable::integers({16});
  index[16] = cpp11::writable::integers({17});
  index[17] = cpp11::writable::integers({18});
  index[18] = cpp11::writable::integers({19});
  index[19] = cpp11::writable::integers({20});
  index[20] = cpp11::writable::integers({21});
  index[21] = cpp11::writable::integers({22});
  index[22] = cpp11::writable::integers({23});
  index[23] = cpp11::writable::integers({24});
  index[24] = cpp11::writable::integers({25});
  index[25] = cpp11::writable::integers({26});
  index[26] = cpp11::writable::integers({27});
  index[27] = integer_sequence(28, 19);
  index[28] = integer_sequence(47, 19);
  index[29] = integer_sequence(66, 19);
  index[30] = integer_sequence(85, 19);
  index[31] = integer_sequence(104, 38);
  index[32] = integer_sequence(142, shared->dim_strain_transmission);
  index[33] = integer_sequence(shared->offset_variable_cum_n_S_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[34] = integer_sequence(shared->offset_variable_cum_n_E_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[35] = integer_sequence(shared->offset_variable_cum_n_I_A_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[36] = integer_sequence(shared->offset_variable_cum_n_I_P_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[37] = integer_sequence(shared->offset_variable_cum_n_R_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[38] = integer_sequence(shared->offset_variable_cum_n_vaccinated + 1, shared->dim_cum_n_S_vaccinated);
  index[39] = integer_sequence(shared->offset_variable_S + 1, shared->dim_cum_n_S_vaccinated);
  index[40] = integer_sequence(shared->offset_variable_prob_strain + 1, shared->dim_E_12);
  index[41] = integer_sequence(shared->offset_variable_I_weighted + 1, shared->dim_cum_n_S_vaccinated);
  index[42] = integer_sequence(shared->offset_variable_tmp_vaccine_probability + 1, shared->dim_cum_n_S_vaccinated);
  index[43] = integer_sequence(shared->offset_variable_T_sero_neg + 1, shared->dim_R);
  index[44] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[45] = integer_sequence(shared->offset_variable_T_PCR_neg + 1, shared->dim_R);
  index[46] = integer_sequence(shared->offset_variable_E + 1, shared->dim_E);
  index[47] = integer_sequence(shared->offset_variable_I_A + 1, shared->dim_I_A);
  index[48] = integer_sequence(shared->offset_variable_I_P + 1, shared->dim_I_P);
  index[49] = integer_sequence(shared->offset_variable_I_C_1 + 1, shared->dim_I_C_1);
  index[50] = integer_sequence(shared->offset_variable_I_C_2 + 1, shared->dim_I_C_2);
  index[51] = integer_sequence(shared->offset_variable_G_D + 1, shared->dim_G_D);
  index[52] = integer_sequence(shared->offset_variable_ICU_pre_unconf + 1, shared->dim_ICU_pre_unconf);
  index[53] = integer_sequence(shared->offset_variable_ICU_pre_conf + 1, shared->dim_ICU_pre_unconf);
  index[54] = integer_sequence(shared->offset_variable_H_R_unconf + 1, shared->dim_H_R_unconf);
  index[55] = integer_sequence(shared->offset_variable_H_R_conf + 1, shared->dim_H_R_unconf);
  index[56] = integer_sequence(shared->offset_variable_H_D_unconf + 1, shared->dim_H_D_unconf);
  index[57] = integer_sequence(shared->offset_variable_H_D_conf + 1, shared->dim_H_D_unconf);
  index[58] = integer_sequence(shared->offset_variable_ICU_W_R_unconf + 1, shared->dim_ICU_W_R_unconf);
  index[59] = integer_sequence(shared->offset_variable_ICU_W_R_conf + 1, shared->dim_ICU_W_R_unconf);
  index[60] = integer_sequence(shared->offset_variable_ICU_W_D_unconf + 1, shared->dim_ICU_W_D_unconf);
  index[61] = integer_sequence(shared->offset_variable_ICU_W_D_conf + 1, shared->dim_ICU_W_D_unconf);
  index[62] = integer_sequence(shared->offset_variable_ICU_D_unconf + 1, shared->dim_ICU_D_unconf);
  index[63] = integer_sequence(shared->offset_variable_ICU_D_conf + 1, shared->dim_ICU_D_unconf);
  index[64] = integer_sequence(shared->offset_variable_W_R_unconf + 1, shared->dim_W_R_unconf);
  index[65] = integer_sequence(shared->offset_variable_W_R_conf + 1, shared->dim_W_R_unconf);
  index[66] = integer_sequence(shared->offset_variable_W_D_unconf + 1, shared->dim_W_D_unconf);
  index[67] = integer_sequence(shared->offset_variable_W_D_conf + 1, shared->dim_W_D_unconf);
  index[68] = integer_sequence(shared->offset_variable_T_sero_pre + 1, shared->dim_T_sero_pre);
  index[69] = integer_sequence(shared->offset_variable_T_sero_pos + 1, shared->dim_T_sero_pos);
  index[70] = integer_sequence(shared->offset_variable_T_PCR_pre + 1, shared->dim_T_PCR_pre);
  index[71] = integer_sequence(shared->offset_variable_T_PCR_pos + 1, shared->dim_T_PCR_pos);
  index.names() = nms;
  size_t len = shared->offset_variable_T_PCR_pos + shared->dim_T_PCR_pos;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
template <>
carehomes::data_t dust_data<carehomes>(cpp11::list data) {
  typedef carehomes::real_t real_t;
  return carehomes::data_t{
      cpp11::as_cpp<real_t>(data["icu"]),
      cpp11::as_cpp<real_t>(data["general"]),
      cpp11::as_cpp<real_t>(data["hosp"]),
      cpp11::as_cpp<real_t>(data["deaths_hosp"]),
      cpp11::as_cpp<real_t>(data["deaths_comm"]),
      cpp11::as_cpp<real_t>(data["deaths_carehomes"]),
      cpp11::as_cpp<real_t>(data["deaths"]),
      cpp11::as_cpp<real_t>(data["deaths_non_hosp"]),
      cpp11::as_cpp<real_t>(data["admitted"]),
      cpp11::as_cpp<real_t>(data["diagnoses"]),
      cpp11::as_cpp<real_t>(data["all_admission"]),
      cpp11::as_cpp<real_t>(data["npos_15_64"]),
      cpp11::as_cpp<real_t>(data["ntot_15_64"]),
      cpp11::as_cpp<real_t>(data["pillar2_pos"]),
      cpp11::as_cpp<real_t>(data["pillar2_tot"]),
      cpp11::as_cpp<real_t>(data["pillar2_cases"]),
      cpp11::as_cpp<real_t>(data["pillar2_over25_pos"]),
      cpp11::as_cpp<real_t>(data["pillar2_over25_tot"]),
      cpp11::as_cpp<real_t>(data["pillar2_over25_cases"]),
      cpp11::as_cpp<real_t>(data["react_pos"]),
      cpp11::as_cpp<real_t>(data["react_tot"]),
      cpp11::as_cpp<real_t>(data["strain_non_variant"]),
      cpp11::as_cpp<real_t>(data["strain_tot"])
    };
}

SEXP dust_carehomes_alloc(cpp11::list r_pars, bool pars_multi, size_t step,
                         cpp11::sexp r_n_particles, size_t n_threads,
                         cpp11::sexp r_seed, cpp11::sexp device_id) {
  return dust_alloc<carehomes>(r_pars, pars_multi, step, r_n_particles,
                               n_threads, r_seed, device_id);
}

SEXP dust_carehomes_run(SEXP ptr, size_t step_end, bool device) {
  return dust_run<carehomes>(ptr, step_end, device);
}

SEXP dust_carehomes_simulate(SEXP ptr, cpp11::sexp step_end) {
  return dust_simulate<carehomes>(ptr, step_end);
}

SEXP dust_carehomes_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust_set_index<carehomes>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step) {
  dust_set_state<carehomes>(ptr, r_state, r_step);
  return R_NilValue;
}

SEXP dust_carehomes_reset(SEXP ptr, cpp11::list r_pars, size_t step) {
  return dust_reset<carehomes>(ptr, r_pars, step);
}

SEXP dust_carehomes_state(SEXP ptr, SEXP r_index) {
  return dust_state<carehomes>(ptr, r_index);
}

size_t dust_carehomes_step(SEXP ptr) {
  return dust_step<carehomes>(ptr);
}

void dust_carehomes_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust_reorder<carehomes>(ptr, r_index);
}

SEXP dust_carehomes_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust_resample<carehomes>(ptr, r_weights);
}

SEXP dust_carehomes_set_pars(SEXP ptr, cpp11::list r_pars) {
  return dust_set_pars<carehomes>(ptr, r_pars);
}

SEXP dust_carehomes_rng_state(SEXP ptr, bool first_only) {
  return dust_rng_state<carehomes>(ptr, first_only);
}

SEXP dust_carehomes_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust_set_rng_state<carehomes>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_carehomes_set_data(SEXP ptr, cpp11::list data) {
  dust_set_data<carehomes>(ptr, data);
  return R_NilValue;
}

SEXP dust_carehomes_compare_data(SEXP ptr) {
  return dust_compare_data<carehomes>(ptr);
}

SEXP dust_carehomes_filter(SEXP ptr, bool save_history) {
  return dust_filter<carehomes>(ptr, save_history);
}

cpp11::sexp dust_carehomes_capabilities() {
  return dust_capabilities<carehomes>();
}

void dust_carehomes_set_n_threads(SEXP ptr, int n_threads) {
  return dust_set_n_threads<carehomes>(ptr, n_threads);
}

int dust_carehomes_n_state(SEXP ptr) {
  return dust_n_state<carehomes>(ptr);
}

cpp11::sexp dust_carehomes_device_info() {
  return dust_device_info<carehomes>();
}
