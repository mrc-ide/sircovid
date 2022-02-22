// Can't use std::min on GPU
template <typename T>
__host__ __device__
T min(const T& a, const T&b) {
  return a < b ? a : b;
}

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type ll_nbinom(real_type data, real_type model, real_type kappa,
                    real_type exp_noise, rng_state_type& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_type mu = model +
    dust::random::exponential<real_type>(rng_state, exp_noise);
  return dust::density::negative_binomial_mu(data, kappa, mu, true);
}

template <typename real_type>
__host__ __device__
real_type ll_binom(real_type data_x, real_type data_size,
                   real_type model_prob) {
  if (std::isnan(data_x) || std::isnan(data_size)) {
    return 0;
  }
  return dust::density::binomial(data_x, data_size, model_prob, true);
}

template <typename real_type>
__host__ __device__
real_type ll_betabinom(real_type data_x, real_type data_size,
                       real_type model_prob, real_type rho) {
  if (std::isnan(data_x) || std::isnan(data_size)) {
    return 0;
  }
  return dust::density::beta_binomial(data_x, data_size, model_prob, rho, true);
}

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type test_prob_pos(real_type pos, real_type neg, real_type sensitivity,
                        real_type specificity, real_type exp_noise,
                        rng_state_type& rng_state) {
  // We add some exponential noise to the number of positives and negatives
  // to help ensure prob_pos is not 0 or 1. If e.g. prob_pos were 0 and there
  // were individuals who tested positive, this would result in a weight of 0
  // for a particle. If all particles have weights of 0, the particle filter
  // breaks. The exponential noise produces small non-zero weights in these
  // circumstances to prevent the particle filter from breaking.
  pos += dust::random::exponential<real_type>(rng_state, exp_noise);
  neg += dust::random::exponential<real_type>(rng_state, exp_noise);
  return (sensitivity * pos + (1 - specificity) * neg) / (pos + neg);
}

// Lots of data items here. Even though some of these are counts (in
// fact, most of them are counts) we'll pull them in as floating point
// numbers (double or float) because that will help with the NA
// handling. Once this is working reliably we could cast these
// explicitly as ints, but there is not likely to be any major gain
// here.
//
// [[odin.dust::compare_data(icu = real_type)]]
// [[odin.dust::compare_data(general = real_type)]]
// [[odin.dust::compare_data(hosp = real_type)]]
// [[odin.dust::compare_data(deaths_hosp = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_0_49 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_50_54 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_55_59 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_60_64 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_65_69 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_70_74 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_75_79 = real_type)]]
// [[odin.dust::compare_data(deaths_hosp_80_plus = real_type)]]
// [[odin.dust::compare_data(deaths_comm = real_type)]]
// [[odin.dust::compare_data(deaths_carehomes = real_type)]]
// [[odin.dust::compare_data(deaths = real_type)]]
// [[odin.dust::compare_data(deaths_non_hosp = real_type)]]
// [[odin.dust::compare_data(admitted = real_type)]]
// [[odin.dust::compare_data(diagnoses = real_type)]]
// [[odin.dust::compare_data(all_admission = real_type)]]
// [[odin.dust::compare_data(all_admission_0_9 = real_type)]]
// [[odin.dust::compare_data(all_admission_10_19 = real_type)]]
// [[odin.dust::compare_data(all_admission_20_29 = real_type)]]
// [[odin.dust::compare_data(all_admission_30_39 = real_type)]]
// [[odin.dust::compare_data(all_admission_40_49 = real_type)]]
// [[odin.dust::compare_data(all_admission_50_59 = real_type)]]
// [[odin.dust::compare_data(all_admission_60_69 = real_type)]]
// [[odin.dust::compare_data(all_admission_70_79 = real_type)]]
// [[odin.dust::compare_data(all_admission_80_plus = real_type)]]
// [[odin.dust::compare_data(sero_pos_15_64_1 = real_type)]]
// [[odin.dust::compare_data(sero_tot_15_64_1 = real_type)]]
// [[odin.dust::compare_data(sero_pos_15_64_2 = real_type)]]
// [[odin.dust::compare_data(sero_tot_15_64_2 = real_type)]]
// [[odin.dust::compare_data(pillar2_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_over25_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_over25_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_over25_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_under15_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_under15_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_under15_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_15_24_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_15_24_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_15_24_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_25_49_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_25_49_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_25_49_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_50_64_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_50_64_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_50_64_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_65_79_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_65_79_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_65_79_cases = real_type)]]
// [[odin.dust::compare_data(pillar2_80_plus_pos = real_type)]]
// [[odin.dust::compare_data(pillar2_80_plus_tot = real_type)]]
// [[odin.dust::compare_data(pillar2_80_plus_cases = real_type)]]
// [[odin.dust::compare_data(react_pos = real_type)]]
// [[odin.dust::compare_data(react_tot = real_type)]]
// [[odin.dust::compare_data(strain_non_variant = real_type)]]
// [[odin.dust::compare_data(strain_tot = real_type)]]
// [[odin.dust::compare_data(strain_over25_non_variant = real_type)]]
// [[odin.dust::compare_data(strain_over25_tot = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
  compare(const typename T::real_type * state,
          const typename T::data_type& data,
          const typename T::internal_type internal,
          std::shared_ptr<const typename T::shared_type> shared,
          typename T::rng_state_type& rng_state) {
    typedef typename T::real_type real_type;

    // State variables; these largely correspond to the quantities in data
    const real_type model_icu = odin(ICU_tot);
    const real_type model_general = odin(general_tot);
    const real_type model_hosp = model_icu + model_general;
    const real_type model_deaths_carehomes = odin(D_carehomes_inc);
    const real_type model_deaths_comm = odin(D_comm_inc);
    const real_type model_deaths_hosp = odin(D_hosp_inc);
    const real_type model_deaths_hosp_0_49 = odin(D_hosp_0_49_inc);
    const real_type model_deaths_hosp_50_54 = odin(D_hosp_50_54_inc);
    const real_type model_deaths_hosp_55_59 = odin(D_hosp_55_59_inc);
    const real_type model_deaths_hosp_60_64 = odin(D_hosp_60_64_inc);
    const real_type model_deaths_hosp_65_69 = odin(D_hosp_65_69_inc);
    const real_type model_deaths_hosp_70_74 = odin(D_hosp_70_74_inc);
    const real_type model_deaths_hosp_75_79 = odin(D_hosp_75_79_inc);
    const real_type model_deaths_hosp_80_plus = odin(D_hosp_80_plus_inc);
    const real_type model_admitted = odin(admit_conf_inc);
    const real_type model_diagnoses = odin(new_conf_inc);
    const real_type model_all_admission_0_9 = odin(all_admission_0_9_conf_inc);
    const real_type model_all_admission_10_19 =
      odin(all_admission_10_19_conf_inc);
    const real_type model_all_admission_20_29 =
      odin(all_admission_20_29_conf_inc);
    const real_type model_all_admission_30_39 =
      odin(all_admission_30_39_conf_inc);
    const real_type model_all_admission_40_49 =
      odin(all_admission_40_49_conf_inc);
    const real_type model_all_admission_50_59 =
      odin(all_admission_50_59_conf_inc);
    const real_type model_all_admission_60_69 =
      odin(all_admission_60_69_conf_inc);
    const real_type model_all_admission_70_79 =
      odin(all_admission_70_79_conf_inc);
    const real_type model_all_admission_80_plus =
      odin(all_admission_80_plus_conf_inc);
    const real_type model_all_admission = model_admitted + model_diagnoses;
    const real_type model_sero_pos_1 = odin(sero_pos_1);
    const real_type model_sero_pos_2 = odin(sero_pos_2);
    const real_type model_sympt_cases = odin(sympt_cases_inc);
    const real_type model_sympt_cases_over25 = odin(sympt_cases_over25_inc);
    const real_type model_sympt_cases_under15 = odin(sympt_cases_under15_inc);
    const real_type model_sympt_cases_15_24 = odin(sympt_cases_15_24_inc);
    const real_type model_sympt_cases_25_49 = odin(sympt_cases_25_49_inc);
    const real_type model_sympt_cases_50_64 = odin(sympt_cases_50_64_inc);
    const real_type model_sympt_cases_65_79 = odin(sympt_cases_65_79_inc);
    const real_type model_sympt_cases_80_plus = odin(sympt_cases_80_plus_inc);
    const real_type model_sympt_cases_non_variant =
      odin(sympt_cases_non_variant_inc);
    const real_type model_sympt_cases_non_variant_over25 =
      odin(sympt_cases_non_variant_over25_inc);
    const real_type model_react_pos = odin(react_pos);

    const int time_p_3 = static_cast<int>(odin(time) + 3);

    const real_type p_NC_today_under15 =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_under15) : odin(p_NC_under15);
    const real_type p_NC_today_15_24 =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_15_24) : odin(p_NC_15_24);
    const real_type p_NC_today_25_49 =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_25_49) : odin(p_NC_25_49);
    const real_type p_NC_today_50_64 =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_50_64) : odin(p_NC_50_64);
    const real_type p_NC_today_65_79 =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_65_79) : odin(p_NC_65_79);
    const real_type p_NC_today_80_plus =
      (time_p_3 % 7 < 2) ?
      odin(p_NC_weekend_80_plus) : odin(p_NC_80_plus);

    const real_type pillar2_under15_negs =
      p_NC_today_under15 * (odin(N_tot_under15) - model_sympt_cases_under15);
    const real_type model_pillar2_under15_prob_pos =
      test_prob_pos(model_sympt_cases_under15,
                    pillar2_under15_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_15_24_negs =
      p_NC_today_15_24 * (odin(N_tot_15_24) - model_sympt_cases_15_24);
    const real_type model_pillar2_15_24_prob_pos =
      test_prob_pos(model_sympt_cases_15_24,
                    pillar2_15_24_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_25_49_negs =
      p_NC_today_25_49 * (odin(N_tot_25_49) - model_sympt_cases_25_49);
    const real_type model_pillar2_25_49_prob_pos =
      test_prob_pos(model_sympt_cases_25_49,
                    pillar2_25_49_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_50_64_negs =
      p_NC_today_50_64 * (odin(N_tot_50_64) - model_sympt_cases_50_64);
    const real_type model_pillar2_50_64_prob_pos =
      test_prob_pos(model_sympt_cases_50_64,
                    pillar2_50_64_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_65_79_negs =
      p_NC_today_65_79 * (odin(N_tot_65_79) - model_sympt_cases_65_79);
    const real_type model_pillar2_65_79_prob_pos =
      test_prob_pos(model_sympt_cases_65_79,
                    pillar2_65_79_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_80_plus_negs =
      p_NC_today_80_plus * (odin(N_tot_80_plus) - model_sympt_cases_80_plus);
    const real_type model_pillar2_80_plus_prob_pos =
      test_prob_pos(model_sympt_cases_80_plus,
                    pillar2_80_plus_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_over25_negs =
      pillar2_25_49_negs + pillar2_50_64_negs +
      pillar2_65_79_negs + pillar2_80_plus_negs;
    const real_type model_pillar2_over25_prob_pos =
      test_prob_pos(model_sympt_cases_over25,
                    pillar2_over25_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type pillar2_negs =
      pillar2_under15_negs + pillar2_15_24_negs + pillar2_over25_negs;
    const real_type model_pillar2_prob_pos =
      test_prob_pos(model_sympt_cases,
                    pillar2_negs,
                    odin(pillar2_sensitivity),
                    odin(pillar2_specificity),
                    odin(exp_noise),
                    rng_state);

    const real_type phi_pillar2_cases_today_under15 =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_under15) : odin(phi_pillar2_cases_under15);
    const real_type phi_pillar2_cases_today_15_24 =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_15_24) : odin(phi_pillar2_cases_15_24);
    const real_type phi_pillar2_cases_today_25_49 =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_25_49) : odin(phi_pillar2_cases_25_49);
    const real_type phi_pillar2_cases_today_50_64 =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_50_64) : odin(phi_pillar2_cases_50_64);
    const real_type phi_pillar2_cases_today_65_79 =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_65_79) : odin(phi_pillar2_cases_65_79);
    const real_type phi_pillar2_cases_today_80_plus =
      (time_p_3 % 7 < 2) ?
      odin(phi_pillar2_cases_weekend_80_plus) : odin(phi_pillar2_cases_80_plus);

    const real_type model_pillar2_under15_cases =
      phi_pillar2_cases_today_under15 * model_sympt_cases_under15;
    const real_type model_pillar2_15_24_cases =
      phi_pillar2_cases_today_15_24 * model_sympt_cases_15_24;
    const real_type model_pillar2_25_49_cases =
      phi_pillar2_cases_today_25_49 * model_sympt_cases_25_49;
    const real_type model_pillar2_50_64_cases =
      phi_pillar2_cases_today_50_64 * model_sympt_cases_50_64;
    const real_type model_pillar2_65_79_cases =
      phi_pillar2_cases_today_65_79 * model_sympt_cases_65_79;
    const real_type model_pillar2_80_plus_cases =
      phi_pillar2_cases_today_80_plus * model_sympt_cases_80_plus;
    const real_type model_pillar2_over25_cases =
      model_pillar2_25_49_cases + model_pillar2_50_64_cases +
      model_pillar2_65_79_cases + model_pillar2_80_plus_cases;
    const real_type model_pillar2_cases =
      model_pillar2_under15_cases + model_pillar2_15_24_cases +
      model_pillar2_over25_cases;

    const real_type model_react_pos_capped =
      min(model_react_pos, odin(N_tot_react));
    const real_type model_react_prob_pos =
      test_prob_pos(model_react_pos_capped,
                    odin(N_tot_react) - model_react_pos_capped,
                    odin(react_sensitivity),
                    odin(react_specificity),
                    odin(exp_noise),
                    rng_state);

    // serology assay 1
    const real_type model_sero_pos_1_capped =
      min(model_sero_pos_1, odin(N_tot_15_64));
    const real_type model_sero_prob_pos_1 =
      test_prob_pos(model_sero_pos_1_capped,
                    odin(N_tot_15_64) - model_sero_pos_1_capped,
                    odin(sero_sensitivity_1),
                    odin(sero_specificity_1),
                    odin(exp_noise),
                    rng_state);

    // serology assay 2
    const real_type model_sero_pos_2_capped =
      min(model_sero_pos_2, odin(N_tot_15_64));
    const real_type model_sero_prob_pos_2 =
      test_prob_pos(model_sero_pos_2_capped,
                    odin(N_tot_15_64) - model_sero_pos_2_capped,
                    odin(sero_sensitivity_2),
                    odin(sero_specificity_2),
                    odin(exp_noise),
                    rng_state);

    // Strain
    real_type strain_sensitivity = 1.0;
    real_type strain_specificity = 1.0;
    const real_type model_strain_prob_pos =
      test_prob_pos(model_sympt_cases_non_variant,
                    model_sympt_cases -
                      model_sympt_cases_non_variant,
                      strain_sensitivity,
                      strain_specificity,
                      odin(exp_noise),
                      rng_state);
    const real_type model_strain_over25_prob_pos =
      test_prob_pos(model_sympt_cases_non_variant_over25,
                    model_sympt_cases_over25 -
                      model_sympt_cases_non_variant_over25,
                      strain_sensitivity,
                      strain_specificity,
                      odin(exp_noise),
                      rng_state);

    // Note that in ll_nbinom, the purpose of exp_noise is to allow a
    // non-zero probability when the model value is 0 and the observed
    // value is non-zero (i.e. there is overreporting)
    const real_type ll_icu =
      ll_nbinom(data.icu, odin(phi_ICU) * model_icu,
                odin(kappa_ICU), odin(exp_noise), rng_state);
    const real_type ll_general =
      ll_nbinom(data.general, odin(phi_general) * model_general,
                odin(kappa_general), odin(exp_noise), rng_state);
    const real_type ll_hosp =
      ll_nbinom(data.hosp, odin(phi_hosp) * model_hosp,
                odin(kappa_hosp), odin(exp_noise), rng_state);

    // We will compute one of the following:
    // 1. ll_deaths_hosp, ll_deaths_carehomes and ll_deaths_comm
    // 2. ll_deaths_hosp and ll_deaths_non_hosp
    // 3. ll_deaths
    const real_type ll_deaths_hosp_0_49 =
      ll_nbinom(data.deaths_hosp_0_49,
                odin(phi_death_hosp) * model_deaths_hosp_0_49,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_50_54 =
      ll_nbinom(data.deaths_hosp_50_54,
                odin(phi_death_hosp) * model_deaths_hosp_50_54,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_55_59 =
      ll_nbinom(data.deaths_hosp_55_59,
                odin(phi_death_hosp) * model_deaths_hosp_55_59,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_60_64 =
      ll_nbinom(data.deaths_hosp_60_64,
                odin(phi_death_hosp) * model_deaths_hosp_60_64,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_65_69 =
      ll_nbinom(data.deaths_hosp_65_69,
                odin(phi_death_hosp) * model_deaths_hosp_65_69,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_70_74 =
      ll_nbinom(data.deaths_hosp_70_74,
                odin(phi_death_hosp) * model_deaths_hosp_70_74,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_75_79 =
      ll_nbinom(data.deaths_hosp_75_79,
                odin(phi_death_hosp) * model_deaths_hosp_75_79,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp_80_plus =
      ll_nbinom(data.deaths_hosp_80_plus,
                odin(phi_death_hosp) * model_deaths_hosp_80_plus,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths_hosp =
      ll_nbinom(data.deaths_hosp, odin(phi_death_hosp) * model_deaths_hosp,
                odin(kappa_death_hosp), odin(exp_noise), rng_state);

    const real_type ll_deaths_carehomes =
      ll_nbinom(data.deaths_carehomes,
                odin(phi_death_carehomes) * model_deaths_carehomes,
                odin(kappa_death_carehomes), odin(exp_noise), rng_state);
    const real_type ll_deaths_comm =
      ll_nbinom(data.deaths_comm, odin(phi_death_comm) * model_deaths_comm,
                odin(kappa_death_comm), odin(exp_noise), rng_state);
    const real_type ll_deaths_non_hosp =
      ll_nbinom(data.deaths_non_hosp,
                odin(phi_death_carehomes) * model_deaths_carehomes +
                  odin(phi_death_comm) * model_deaths_comm,
                  odin(kappa_death_non_hosp), odin(exp_noise), rng_state);
    const real_type ll_deaths =
      ll_nbinom(data.deaths,
                odin(phi_death_hosp) * model_deaths_hosp +
                  odin(phi_death_carehomes) * model_deaths_carehomes +
                  odin(phi_death_comm) * model_deaths_comm,
                  odin(kappa_death), odin(exp_noise), rng_state);

    const real_type ll_admitted =
      ll_nbinom(data.admitted, odin(phi_admitted) * model_admitted,
                odin(kappa_admitted), odin(exp_noise), rng_state);
    const real_type ll_diagnoses =
      ll_nbinom(data.diagnoses, odin(phi_diagnoses) * model_diagnoses,
                odin(kappa_diagnoses), odin(exp_noise), rng_state);
    const real_type ll_all_admission =
      ll_nbinom(data.all_admission,
                odin(phi_all_admission) * model_all_admission,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_0_9 =
      ll_nbinom(data.all_admission_0_9,
                odin(phi_all_admission) * model_all_admission_0_9,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_10_19 =
      ll_nbinom(data.all_admission_10_19,
                odin(phi_all_admission) * model_all_admission_10_19,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_20_29 =
      ll_nbinom(data.all_admission_20_29,
                odin(phi_all_admission) * model_all_admission_20_29,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_30_39 =
      ll_nbinom(data.all_admission_30_39,
                odin(phi_all_admission) * model_all_admission_30_39,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_40_49 =
      ll_nbinom(data.all_admission_40_49,
                odin(phi_all_admission) * model_all_admission_40_49,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_50_59 =
      ll_nbinom(data.all_admission_50_59,
                odin(phi_all_admission) * model_all_admission_50_59,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_60_69 =
      ll_nbinom(data.all_admission_60_69,
                odin(phi_all_admission) * model_all_admission_60_69,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_70_79 =
      ll_nbinom(data.all_admission_70_79,
                odin(phi_all_admission) * model_all_admission_70_79,
                odin(kappa_all_admission), odin(exp_noise), rng_state);
    const real_type ll_all_admission_80_plus =
      ll_nbinom(data.all_admission_80_plus,
                odin(phi_all_admission) * model_all_admission_80_plus,
                odin(kappa_all_admission), odin(exp_noise), rng_state);

    const real_type ll_serology_1 =
      ll_binom(data.sero_pos_15_64_1, data.sero_tot_15_64_1,
               model_sero_prob_pos_1);
    const real_type ll_serology_2 =
      ll_binom(data.sero_pos_15_64_2, data.sero_tot_15_64_2,
               model_sero_prob_pos_2);

    const real_type ll_pillar2_tests =
      ll_betabinom(data.pillar2_pos, data.pillar2_tot,
                   model_pillar2_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_over25_tests =
      ll_betabinom(data.pillar2_over25_pos, data.pillar2_over25_tot,
                   model_pillar2_over25_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_under15_tests =
      ll_betabinom(data.pillar2_under15_pos, data.pillar2_under15_tot,
                   model_pillar2_under15_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_15_24_tests =
      ll_betabinom(data.pillar2_15_24_pos, data.pillar2_15_24_tot,
                   model_pillar2_15_24_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_25_49_tests =
      ll_betabinom(data.pillar2_25_49_pos, data.pillar2_25_49_tot,
                   model_pillar2_25_49_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_50_64_tests =
      ll_betabinom(data.pillar2_50_64_pos, data.pillar2_50_64_tot,
                   model_pillar2_50_64_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_65_79_tests =
      ll_betabinom(data.pillar2_65_79_pos, data.pillar2_65_79_tot,
                   model_pillar2_65_79_prob_pos, odin(rho_pillar2_tests));
    const real_type ll_pillar2_80_plus_tests =
      ll_betabinom(data.pillar2_80_plus_pos, data.pillar2_80_plus_tot,
                   model_pillar2_80_plus_prob_pos, odin(rho_pillar2_tests));

    const real_type ll_pillar2_under15_cases =
      ll_nbinom(data.pillar2_under15_cases,
                model_pillar2_under15_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_15_24_cases =
      ll_nbinom(data.pillar2_15_24_cases,
                model_pillar2_15_24_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_25_49_cases =
      ll_nbinom(data.pillar2_25_49_cases,
                model_pillar2_25_49_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_50_64_cases =
      ll_nbinom(data.pillar2_50_64_cases,
                model_pillar2_50_64_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_65_79_cases =
      ll_nbinom(data.pillar2_65_79_cases,
                model_pillar2_65_79_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_80_plus_cases =
      ll_nbinom(data.pillar2_80_plus_cases,
                model_pillar2_80_plus_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_over25_cases =
      ll_nbinom(data.pillar2_over25_cases,
                model_pillar2_over25_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);
    const real_type ll_pillar2_cases =
      ll_nbinom(data.pillar2_cases,
                model_pillar2_cases,
                odin(kappa_pillar2_cases), odin(exp_noise), rng_state);

    const real_type ll_react =
      ll_binom(data.react_pos, data.react_tot,
               model_react_prob_pos);
    const real_type ll_strain =
      ll_binom(data.strain_non_variant, data.strain_tot,
               model_strain_prob_pos);
    const real_type ll_strain_over25 =
      ll_binom(data.strain_over25_non_variant, data.strain_over25_tot,
               model_strain_over25_prob_pos);

    return ll_icu + ll_general + ll_hosp + ll_deaths_hosp_0_49 +
      ll_deaths_hosp_50_54 + ll_deaths_hosp_55_59 + ll_deaths_hosp_60_64 +
      ll_deaths_hosp_65_69 + ll_deaths_hosp_70_74 + ll_deaths_hosp_75_79 +
      ll_deaths_hosp_80_plus + ll_deaths_hosp + ll_deaths_carehomes +
      ll_deaths_comm + ll_deaths_non_hosp + ll_deaths + ll_admitted +
      ll_diagnoses + ll_all_admission +
      ll_all_admission_0_9 + ll_all_admission_10_19 + ll_all_admission_20_29 +
      ll_all_admission_30_39 + ll_all_admission_40_49 + ll_all_admission_50_59 +
      ll_all_admission_60_69 + ll_all_admission_70_79 +
      ll_all_admission_80_plus + ll_serology_1 + ll_serology_2 +
      ll_pillar2_tests + ll_pillar2_over25_tests + ll_pillar2_under15_tests +
      ll_pillar2_15_24_tests + ll_pillar2_25_49_tests + ll_pillar2_50_64_tests +
      ll_pillar2_65_79_tests + ll_pillar2_80_plus_tests +
      ll_pillar2_cases + ll_pillar2_over25_cases + ll_pillar2_under15_cases +
      ll_pillar2_15_24_cases + ll_pillar2_25_49_cases + ll_pillar2_50_64_cases +
      ll_pillar2_65_79_cases + ll_pillar2_80_plus_cases +
      ll_react + ll_strain + ll_strain_over25;
  }
