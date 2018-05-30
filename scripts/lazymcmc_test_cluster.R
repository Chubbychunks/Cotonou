# good_previous_fit = get(load("best_pars_combined_0904_top.Rdata"))
#
# good_previous_fit_variedpars = good_previous_fit[[1]][rownames(ranges)]
#
# lazymcmc_for_cluster(good_previous_fit_variedpars = good_previous_fit_variedpars,
#                      ranges = ranges,
#                      best_set = best_set, time = time,
#                      par_seq = par_seq, condom_seq = condom_seq,
#                      groups_seq = groups_seq, years_seq = years_seq, outputs = outputs,
#                      prev_points = prev_points_FSW_Cotonou_centrale_lower_bound,
#                      frac_N_discard_points = frac_N_discard_points_no_FSW_LB,
#                      Ntot_data_points = Ntot_data_points,
#                      ART_data_points = ART_data_points_lazymcmc,
#                      PrEP_fitting = PrEP_fitting)



lazymcmc_for_cluster <- function(good_previous_fit_variedpars, ranges, best_set, time,
                                 par_seq, condom_seq, groups_seq, years_seq, outputs,
                                 prev_points,
                                 frac_N_discard_points,
                                 Ntot_data_points,
                                 ART_data_points,
                                 PrEP_fitting) {

  # CREATING PARTAB WHICH IS THE INPUT DATAFRAME
  parTab_create = data.frame(
    names = rownames(ranges),
    values = ranges[,1],
    fixed = 0,
    steps = 0.1,
    lower_bound = ranges[,1],
    upper_bound = ranges[,2]

  )


  j=0
  for(i in 1:length(rownames(ranges)))
  {
    if(parTab_create[i,"names"] %in% names(good_previous_fit_variedpars))
    {
      parTab_create[i, "values"] = as.numeric(good_previous_fit_variedpars[as.character(parTab_create[i,"names"])][[1]][1]);
      j=j+1
    }


  }

  parTab_create[parTab_create$lower_bound ==parTab_create$upper_bound, "fixed"] = 1
  rownames(parTab_create) = NULL

  ###

  create_lik <- function(parTab, ranges, best_set, time,
                         par_seq, condom_seq, groups_seq, years_seq, outputs,
                         prev_points,
                         frac_N_discard_points,
                         Ntot_data_points,
                         ART_data_points,
                         PrEP_fitting,
                         PRIOR_FUNC,...){


    par_names <- rownames(ranges)

    ## using the `...` bit.
    likelihood_func <- function(pars){

      ranges_new = cbind(pars, pars)
      rownames(ranges_new) = par_names

      parameters = cotonou::lhs_parameters(1, set_pars = best_set, Ncat = 9, time = time,
                                           ranges = ranges_new, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)

      res = cotonou::return_outputs(parameters[[1]], cotonou::main_model, time = time, outputs = outputs)

      lik = cotonou::likelihood_lazymcmc(res, time = time,
                                         prev_points = prev_points_FSW_Cotonou_centrale_lower_bound,
                                         frac_N_discard_points = frac_N_discard_points_no_FSW_LB,
                                         Ntot_data_points = Ntot_data_points,
                                         ART_data_points = ART_data_points_lazymcmc,
                                         PrEP_fitting = PrEP_fitting)

      return(lik)
    }
    return(likelihood_func)
  }

  f <- create_lik(parTab,ranges, best_set, time,
                  par_seq, condom_seq, groups_seq, years_seq, outputs,
                  prev_points,
                  frac_N_discard_points,
                  Ntot_data_points,
                  ART_data_points,
                  PrEP_fitting)
  f(parTab$values)





}
