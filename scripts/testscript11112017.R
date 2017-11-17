
odin::odin_package(".") # looks for any models inside inst/odin
devtools::load_all()

PrEP_fitting = data.frame(time = c(2016, 2017),
                          group = c("S1a", "S1c"),
                          lower = c(50, 40),
                          point = c(55, 400),
                          upper = c(61, 66)


)


# test fits ---------------------------------------------------------------


frac_N_discard_points_test = data.frame(variable = c("Pro FSW"),
                                        min = c(0),
                                        max = c(1))
ART_data_points_test = data.frame(time = c(2014),
                                  Lower = c(0),
                                  Upper = c(1),
                                  variable = c("Pro FSW"))
prev_points_test = data.frame(time = c(2015),
                              variable = c(rep("Pro FSW", 1)),
                              value = c(0),
                              lower = c(0),

                              upper = c(1))

Ntot_data_points_test = data.frame(time = c(1992, 2002, 2013, 2020, 2030),
                                   point = c(10, 10, 10, 10, 10),
                                   lower = c(10, 10, 10, 10, 10),
                                   upper = c(10000000000, 10000000000, 10000000000, 10000000000, 10000000000),
                                   colour = c("data", "data", "data", "predicted", "predicted"))
# # -----------------------------------------------------------------------


number_of_prep_samples = 5

prep_ranges <- rbind(
  eP1a = c(0.47, 0.98),
  eP1b = c(0.2, 0.4),
  psia = c(1, 4),
  psib = c(1, 4),
  prep_offering_rate = c(1, 10),
  prep_dropout = c(1, 4),
  PrEPOnOff = c(1,1)

)



best_pars_combined = result[[1]]



res_best_runs_after_prep_fit = lapply(best_pars_combined, function(x) {



  # x[rownames(ranges)]

  x = lapply(x[rownames(ranges)], function(y) {
    if(length(y) == 9)
      y = y[1] else y
  })

  x = x[-which(names(x) %in% rownames(prep_ranges))]



  combined_ranges = cbind(unlist(x[rownames(ranges)]), unlist(x[rownames(ranges)]))

  combined_ranges = rbind(combined_ranges, prep_ranges)


  # i now have the ranges, have to put them into the function below
  #test fitting to nothing first to see if works
  #then add prep fitting


  # res_after_prep = run_model_with_fit(number_simulations = number_of_prep_samples, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq, best_set = best_set, time = time,
  #                                            ranges = combined_ranges, outputs = outputs,
  #                                            prev_points = prev_points_test, frac_N_discard_points = frac_N_discard_points_test,
  #                                            Ntot_data_points = Ntot_data_points_test, ART_data_points = ART_data_points_test, PrEP_fitting = PrEP_fitting)
  res_after_prep = run_model_with_fit_multiple(batch_size = number_of_prep_samples, number_simulations = number_of_prep_samples, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq, best_set = best_set, time = time,
                                               ranges = combined_ranges, outputs = outputs,
                                               prev_points = prev_points_test, frac_N_discard_points = frac_N_discard_points_test,
                                               Ntot_data_points = Ntot_data_points_test, ART_data_points = ART_data_points_test, PrEP_fitting = PrEP_fitting)

  # return(list(res_after_prep[[1]], res_after_prep[[3]], res_after_prep[[6]]))
  return(list(res_after_prep[[1]], res_after_prep[[2]], res_after_prep[[5]]))

  # return(combined_ranges)

}
)

# 1. best fit to all data
lapply(res_best_runs_after_prep_fit, function(x) x[[1]])

best_pars_after_prep = lapply(res_best_runs_after_prep_fit, function(x) {
  prepfitsout = x[[3]]
  print(x[[3]])
  bla = which(prepfitsout == min(prepfitsout))

  return(x[[2]][bla])

}
)

# 2. just to show best prep parameters
best_pars_after_prep_just_rep = lapply(best_pars_after_prep, function(x) x[[1]][rownames(prep_ranges)])

# 3. check that it fits to all the data



# 4. sample DALY weights and costs


best_pars_after_prep = lapply(best_pars_after_prep, function(x) x = x[[1]])


number_of_cost_DALY_samples = length(best_pars_after_prep)

cost_DALY_ranges <- rbind(
  cost_ART_initiation_study_FSW = c(173, 303),
  # cost_ART_initiation_gov_FSW = c(40, 50),

  cost_1_year_ART_study_FSW = c(100, 110),
  # cost_1_year_ART_gov_FSW = c(100, 110),


  cost_1_year_ART_rest = c(0.00, 66.38),


  # cost_PrEP_initiation = c(25, 29),
  # cost_1_year_PrEP_perfect = c(60, 70),
  # cost_1_year_PrEP_intermediate = c(50, 60),
  # cost_1_year_PrEP_non = c(40, 50),


  W0 = c(1, 1),
  W1 = c(1 - 0.111, 1 - 0.052),
  W3 = c(1 - 0.377, 1 - 0.184),
  W2 = c(1 - 0.743, 1 - 0.406)



)


CEA_outputs = unique(c("W0", "W1", "W2", "W3", "Number_DALY_W1","Number_DALY_W2", "Number_DALY_W3", "FSW_On_PrEP_all_cats", "PrEPinitiations", "PrEPinitiations1a",
                       "PrEPinitiations1b", "PrEPinitiations1c", "pc_susceptible_FSW_On_PrEP", "pc_all_FSW_On_PrEP", "Number_Susceptibles",
                       "HIV_positive_On_ART", "HIV_positive_Diagnosed_Off_ART", "Primary_Off_ART",
                       "CD4_above_500_Off_ART", "CD4_350_500_Off_ART", "CD4_200_350_Off_ART", "CD4_below_200_Off_ART", "cumuDeaths_On_ART", "HIV_positive", "ec", "cumuARTinitiations","cumuARTREinitiations", "rate_leave_pro_FSW","tau_intervention",
                "testing_prob", "tau", "N", "S0", "S1a", "S1b", "S1c", "S1d", "I01", "I11", "I02", "I03", "I04",
                "I05", "I22", "I23", "I24", "I25", "I32", "I33", "I34", "I35",  "I42", "I43", "I44", "I45", "prev",
                "frac_N", "Ntot", "epsilon", "rate_leave_client", "alphaItot", "prev_FSW", "prev_LowFSW", "prev_client",
                "prev_men", "prev_women", "c_comm_balanced", "c_noncomm_balanced", "who_believe_comm", "ART_coverage_FSW",
                "ART_coverage_men", "ART_coverage_women", "ART_coverage_all", "rho", "n_comm", "n_noncomm", "fc_comm",
                "fc_noncomm", "N", "cumuHIVDeaths", "lambda_sum_0", "lambda_sum_1a", "lambda_sum_1b", "lambda_sum_1c",
                "lambda_sum_1d", "S0", "S1a", "S1b", "S1c", "S1d", "OnPrEP1a", "OnPrEP1b",
                "OnPrEP1c", "ART_eligible_CD4_above_500", "ART_eligible_CD4_350_500","ART_eligible_CD4_200_349","ART_eligible_CD4_below_200",
                "cumuAllDeaths", "cumuHIVDeaths", "cumuARTinitiations", "cumuARTREinitiations",
                "OnPrEP", "ART_sex_ratio", "pc_S1b", "pc_S1a", "pc_S1c", "cumuInf",
                "intervention_ART_increase", "testing_prob", "rho_intervention",
                "ART_eligible_CD4_above_500", "ART_eligible_CD4_350_500", "ART_eligible_CD4_200_349",
                "ART_eligible_CD4_below_200", "new_people_in_group_FSW_only", "rate_move_out", "rate_move_in",
                "FSW_out", "FSW_in", "zeta", "tau", "prep_offering_rate", "intervention_testing_increase", "sigma",
                "PrEPOnOff", "prev", "frac_N", "Ntot", "epsilon", "rate_leave_client", "alphaItot", "prev_FSW",
                "prev_LowFSW", "prev_client", "prev_men", "prev_women", "c_comm_balanced", "c_noncomm_balanced",
                "who_believe_comm", "ART_coverage_FSW", "ART_coverage_men", "ART_coverage_women", "ART_coverage_all",
                "rho", "n_comm", "n_noncomm", "fc_comm", "fc_noncomm", "N", "cumuHIVDeaths", "lambda_0", "lambda_1a",
                "lambda_1b", "lambda_1c", "lambda_1d"))
#
#
# CEA_outputs = c("prev", "HIV_positive", "lambda_0", "lambda_1a",
#                 "lambda_1b", "lambda_1c", "lambda_1d",
#                 )

epi_start = 1986
# epi_end = 2030
epi_end = 2017

time <- seq(epi_start, epi_end, length.out = (epi_end - epi_start + 0.5)*2)

mid_years_vec = time %% 1 != 0





res_best_runs_after_prep_fit_sampling_costs_n_DALYs = lapply(best_pars_after_prep, function(x) {

  unique_pars = unique(c(rownames(ranges), rownames(prep_ranges)))

  # # x[rownames(ranges)]
  #
  x = lapply(x[unique_pars], function(y) {
    if(length(y) == 9)
      y = y[1] else y
  })

  combined_ranges = cbind(unlist(x[unique_pars]), unlist(x[unique_pars]))
  combined_ranges = rbind(combined_ranges, cost_DALY_ranges)


  result <- cotonou::run_model(number_simulations = number_of_cost_DALY_samples, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq,
                               years_seq = years_seq, best_set = best_set, time = time, ranges = combined_ranges, outputs = CEA_outputs)


  return(result[[2]]) # returning outcomes

}
)



total_DALYs = lapply(res_best_runs_after_prep_fit_sampling_costs_n_DALYs, function(x) {
  lapply(x, function(y) {

    # ____  _____ __  __  _____     _______   _____ ______        __   ___  _   _ _____ ____ ___ ____  _____
    # |  _ \| ____|  \/  |/ _ \ \   / / ____| |  ___/ ___\ \      / /  / _ \| | | |_   _/ ___|_ _|  _ \| ____|
    # | |_) |  _| | |\/| | | | \ \ / /|  _|   | |_  \___ \\ \ /\ / /  | | | | | | | | | \___ \| || | | |  _|
    # |  _ <| |___| |  | | |_| |\ V / | |___  |  _|  ___) |\ V  V /   | |_| | |_| | | |  ___) | || |_| | |___
    # |_| \_\_____|_|  |_|\___/  \_/  |_____| |_|   |____/  \_/\_/     \___/ \___/  |_| |____/___|____/|_____|
    #


    # y$W2[1]

    y$Number_Susceptibles[mid_years_vec,-9] * y$W0[1] +
      y$Number_DALY_W1[mid_years_vec,-9] * y$W1[1] +
      y$Number_DALY_W2[mid_years_vec,-9] * y$W2[1] +
      y$Number_DALY_W3[mid_years_vec,-9] * y$W3[1]

  })
})
total_DALYs
## should check results here















# checking

#
# which_original_fit = 1
# res_test = lapply(res_best_runs_after_prep_fit[[which_original_fit]][[2]], cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)
#
#
#
# S0_indiv = data.frame(time, t(do.call(rbind, lapply(res_test, function(x) x$S0[,1]))))
# S1a_indiv = data.frame(time, t(do.call(rbind, lapply(res_test, function(x) x$S1a[,1]))))
# S1b_indiv = data.frame(time, t(do.call(rbind, lapply(res_test, function(x) x$S1b[,1]))))
# S1c_indiv = data.frame(time, t(do.call(rbind, lapply(res_test, function(x) x$S1c[,1]))))
# S1d_indiv = data.frame(time, t(do.call(rbind, lapply(res_test, function(x) x$S1d[,1]))))
#
# all_prep_cats = rbind(S1a_indiv, S1b_indiv, S1c_indiv)
#
# all_prep_cats$group = rep(c("S1a", "S1b", "S1c"), each = length(time))
#
# all_prep_cats_melted = reshape2::melt(all_prep_cats, id.vars = c("time", "group"))
#
# colnames(all_prep_cats_melted) = c("time", "group", "variable", "point")
#
#
# ggplot() + geom_line(data = all_prep_cats_melted, aes(x = time, y = point, factor = variable, colour = group)) +
#   theme_bw()+ geom_point(data = PrEP_fitting, aes(x = time, y = point, colour = group)) + facet_wrap(~variable)+
#   geom_blank(data = prep_axes, aes( y = value))
#
#
#
#
#
#
# # for each original fit, taking the pars which best fit to prep
# unlist(lapply(lapply(res_best_runs_after_prep_fit, function(x) x[[3]]), function(x) {
#   # return(x[which(x == min(x))])
#   return(which(x == min(x)))
# }))
# #
# #
# #
# #
# #
# # lapply(res_best_runs_after_prep_fit, function(x) x[[2]])[[2]][[4]][rownames(prep_ranges)]
# #
#



