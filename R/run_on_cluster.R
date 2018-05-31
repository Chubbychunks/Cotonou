#' @export
#' @useDynLib cotonou
just_parameters <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs) {


  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)
  # end of parameters --------------------------------------------------------------

  return(parameters)

}



#' @export
#' @useDynLib cotonou
quantile_95 <- function(x) return(quantile(x, probs = c(0.025, 0.5, 0.975)))


#' @export
#' @useDynLib cotonou
return_outputs <- function(p, gen, time, outputs, solving_method = "lsoda") {
  mod <- gen(user = p)
  all_results <- mod$transform_variables(mod$run(time,  rtol = 10^-4,  atol = 10^-4, method = solving_method))
  # counter <<- counter + 1
  # if (counter %% 10 == 0)
  #   print(counter)


  return(all_results[outputs])
}

#' @export
#' @useDynLib cotonou
return_all_outputs <- function(p, gen, time) {
  mod <- gen(user = p)
  return(mod$transform_variables(mod$run(time)))
}

#' @export
#' @useDynLib cotonou
likelihood_rough <- function(x, time, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting) {
  the_prev = data.frame(time, x$prev_FSW, x$prev_LowFSW, x$prev_client, x$prev_women, x$prev_men)
  names(the_prev) = c("time", "Pro FSW", "Low-level FSW", "Clients", "Women", "Men")

  the_frac_N = data.frame(time, x$frac_N[,c(1, 5, 7, 8)], x$frac_N[,1] + x$frac_N[,2], x$frac_N[,2]/ x$frac_N[,1])
  names(the_frac_N) = c("time", "Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio")


  message = "nothing"

  likelihood_count <- 0
  ##
  frac_count <- 0

  if(all(!is.na(unlist(the_frac_N)))) {
    for (i in 1:length(frac_N_discard_points[,1]))
    {
      if(all(the_frac_N[, as.character(frac_N_discard_points[i, "variable"])] < frac_N_discard_points[i, "max"]) &&
         all(the_frac_N[,  as.character(frac_N_discard_points[i, "variable"])] > frac_N_discard_points[i, "min"])) {
        frac_count <- frac_count + 1
      }
    }
  }

  # if the run doesn't fit within the Ntot CIs, then set the frac_count to 0 so the run doesn't pass
  if(!(all(!is.na(x$Ntot)) &&
       x$Ntot[which(time == 1992)] < Ntot_data_points[1, "upper"] && x$Ntot[which(time == 1992)] > Ntot_data_points[1, "lower"] &&
       x$Ntot[which(time == 2002)] < Ntot_data_points[2, "upper"] && x$Ntot[which(time == 2002)] > Ntot_data_points[2, "lower"] &&
       x$Ntot[which(time == 2013)] < Ntot_data_points[3, "upper"] && x$Ntot[which(time == 2013)] > Ntot_data_points[3, "lower"])) {
    frac_count <- 0
    message = "Ntot"
  }

  if(frac_count != length(frac_N_discard_points[,1]))
    message = "frac count"

  prev_fits = c()

  if(frac_count == length(frac_N_discard_points[,1])) {
    # prevalence
    for(i in 1:length(prev_points[,1]))
    {

      point = subset(the_prev, time == prev_points[i, "time"], select = as.character(prev_points[i, "variable"]))
      # point = the_prev[the_prev$time == prev_points[i, "time"], as.character(prev_points[i, "variable"])]
      if(!is.na(point)) {{if((point < prev_points[i, "upper"]) && (point > prev_points[i, "lower"]))
      {
        likelihood_count <- likelihood_count + 1
        prev_fits <- c(prev_fits, i)
      }}}
    }
  }



  if("All" %in% levels(ART_data_points$variable))
  {

    # ART_ratio = x$Women_on_ART/x$Men_on_ART
    # if(all(ART_ratio[!is.na(ART_ratio)] > 1 & ART_ratio[!is.na(ART_ratio)] < 2))
    #   likelihood_count <- likelihood_count + 1

    ART_data_points_allgroups = ART_data_points[ART_data_points$variable == "All",]
    # fitting to ART cov
    if(likelihood_count > 0)
    {
      if(all(!is.na(x$ART_coverage_all))){

        for(i in 1:length(ART_data_points_allgroups[,1]))
        {
          the_time = ART_data_points_allgroups[i, "time"]
          if(x$ART_coverage_all[which(time == the_time)] > ART_data_points_allgroups[i, "Lower"] && x$ART_coverage_all[which(time == the_time)] < ART_data_points_allgroups[i, "Upper"]) {
            likelihood_count <- likelihood_count + 1
          }
        }
      }
    }
  }

  if("Pro FSW" %in% levels(ART_data_points$variable))
  {



    ART_data_points_FSW = ART_data_points[ART_data_points$variable == "Pro FSW",]
    # fitting to ART cov
    if(likelihood_count > 0)
    {
      if(all(!is.na(x$ART_coverage_FSW))){
        for(i in 1:length(ART_data_points_FSW[,1]))
        {
          the_time = ART_data_points_FSW[i, "time"]
          if(x$ART_coverage_FSW[which(time == the_time)] > ART_data_points_FSW[i, "Lower"] && x$ART_coverage_FSW[which(time == the_time)] < ART_data_points_FSW[i, "Upper"]) {
            likelihood_count <- likelihood_count + 1
          }
        }
      }
    }
  }



  if("Numbers FSW" %in% levels(ART_data_points$variable))
  {

    ART_data_points_FSW = ART_data_points[ART_data_points$variable == "Numbers FSW",]
    # fitting to ART cov
    if(likelihood_count > 0)
    {
      if(all(!is.na(x$ART_coverage_FSW))){
        for(i in 1:length(ART_data_points_FSW[,1]))
        {
          the_time = ART_data_points_FSW[i, "time"]


          if(x$HIV_positive_On_ART[which(time == the_time),1]  > ART_data_points_FSW[i, "Lower"] && x$HIV_positive_On_ART[which(time == the_time),1] < ART_data_points_FSW[i, "Upper"]) {

            likelihood_count <- likelihood_count + 1


          }
        }
      }
    }





  }






  prep_tasp_fit = Inf
  if(!is.null(PrEP_fitting) && sum(time %% 1 != 0) > 1)
  {





    # PrEP

    if(length(time == 589)) {
      total_on_prep = data.frame(time, x["FSW_On_PrEP_all_cats"])

      PY_PrEP =
        total_on_prep[which(time == 2015 + 1/12),2]/12 +
        total_on_prep[which(time == 2015 + 2/12),2]/12 +
        total_on_prep[which(time == 2015 + 3/12),2]/12 +
        total_on_prep[which(time == 2015 + 4/12),2]/12 +
        total_on_prep[which(time == 2015 + 5/12),2]/12 +
        total_on_prep[which(time == 2015 + 6/12),2]/12 +
        total_on_prep[which(time == 2015 + 7/12),2]/12 +
        total_on_prep[which(time == 2015 + 8/12),2]/12 +
        total_on_prep[which(time == 2015 + 9/12),2]/12 +
        total_on_prep[which(time == 2015 + 10/12),2]/12 +
        total_on_prep[which(time == 2015 + 11/12),2]/12 +
        total_on_prep[which(time == 2016),2]/12 +

        total_on_prep[which(time == 2016 + 1/12),2]/12 +
        total_on_prep[which(time == 2016 + 2/12),2]/12 +
        total_on_prep[which(time == 2016 + 3/12),2]/12 +
        total_on_prep[which(time == 2016 + 4/12),2]/12 +
        total_on_prep[which(time == 2016 + 5/12),2]/12 +
        total_on_prep[which(time == 2016 + 6/12),2]/12 +
        total_on_prep[which(time == 2016 + 7/12),2]/12 +
        total_on_prep[which(time == 2016 + 8/12),2]/12 +
        total_on_prep[which(time == 2016 + 9/12),2]/12 +
        total_on_prep[which(time == 2016 + 10/12),2]/12 +
        total_on_prep[which(time == 2016 + 11/12),2]/12 +
        total_on_prep[which(time == 2017),2]/12




      PrEPinitiations = x["PrEPinitiations"][[1]][,1]

      # data.frame(time, PrEPinitiations)

      Total_PrEPinitiations = PrEPinitiations[which(time == 2017)] -
        PrEPinitiations[which(time == 2015)]


      number_on_prep_end_of_study = total_on_prep[which(time == 2017),2]


      # TasP

      TasP_initiations = x["TasPinitiations"][[1]][,1][which(time == 2017)] - x["TasPinitiations"][[1]][,1][which(time == 2015)]

      Number_on_ART_end_of_study = x["HIV_positive_On_ART"][[1]][,1][which(time == 2017)]

      Number_on_ART_end_of_year1 = x["HIV_positive_On_ART"][[1]][,1][which(time == 2016)]



      prep_tasp_fit = (PY_PrEP-250)^2 + (Total_PrEPinitiations - 256)^2 + (number_on_prep_end_of_study - 121)^2 +
        (TasP_initiations - 107)^2 + (Number_on_ART_end_of_study - 137)^2 + (Number_on_ART_end_of_year1 - 122)^2


    }







    # PY_PrEP = total_on_prep[which(time == 2015.5),2] +
    #   total_on_prep[which(time == 2016.5),2]
    #
    # prep_fit = (PY_PrEP-250)^2


    #
    #     prep_fit = 0;
    #
    #     for(i in 1:length(PrEP_fitting[,1]))
    #     {
    #       time = PrEP_fitting[i, "time"]
    #       if(PrEP_fitting[i, "group"] == "S1a")
    #         prep_fit = prep_fit + (S1a[S1a$time == time, 1] - PrEP_fitting[i, "point"])^2
    #       if(PrEP_fitting[i, "group"] == "S1b")
    #         prep_fit = prep_fit + (S1b[S1b$time == time, 1] - PrEP_fitting[i, "point"])^2
    #       if(PrEP_fitting[i, "group"] == "S1c")
    #         prep_fit = prep_fit + (S1c[S1c$time == time, 1] - PrEP_fitting[i, "point"])^2
    # }


  }




  return (list(likelihood_count, prev_fits, message, prep_tasp_fit))
  # return (list(likelihood_count, frac_count))

}

#' @export
#' @useDynLib cotonou
run_model_for_tests <- function(number_simulations, time, parameters) {

  return(lapply(parameters, return_all_outputs, main_model, time = time))

}

#' @export
#' @useDynLib cotonou
run_model <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs) {


  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)
  # end of parameters --------------------------------------------------------------


  res = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)


  return(list(parameters, res))

}

#' @export
#' @useDynLib cotonou
run_model_with_fit <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {


  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


  # this is the slowest part - can I turn this into
  res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)



  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)

  sorted_likelihood_list = sort(unlist(lapply(likelihood_list, function(x) x[[1]])))

  best_runs = which(unlist(lapply(likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))

  out <- res[best_runs]

  message_list <- unlist(lapply(likelihood_list, function(x) x[[3]]))

  prep_fit <- unlist(lapply(likelihood_list, function(x) x[[4]]))

  return(list(parameters[best_runs], likelihood_list, out, best_runs, message_list, prep_fit))


}

#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {


  # LHS to create parameter sets
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


  # this is the slowest part - simulating model
  res = parallel::parLapply(NULL, parameters, return_outputs, main_model, time = time, outputs = outputs)



  # model fitting
  likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)

  sorted_likelihood_list = sort(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])))

  best_runs = which(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))

  out <- res[best_runs]

  return(list(parameters[best_runs], likelihood_list, out, best_runs))


}







#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster_multiple <- function(batch_size, number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {




  best_fit_pars = list()
  max_fit = 1

  ### minus 1
  best_fit_pars_minus_1 = list()
  max_fit_minus_1 = 0

  # results_list = list()
  for(i in 1:(number_simulations/batch_size))
  {
    # LHS to create parameter sets
    parameters <- cotonou::lhs_parameters_parallel(batch_size, set_pars = best_set, Ncat = 9, time = time,
                                                   ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


    # pars = parameters[(batch_size * (i - 1) + 1):(batch_size * i)]

    # this is the slowest part - simulating model
    res = parallel::parLapply(NULL, parameters, cotonou::return_outputs, main_model, time = time, outputs = outputs)


    if(all(lapply(lapply(res, function(x) x$prev[,1]), length) == length(time)))
    {
      # model fitting
      likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)

      sorted_likelihood_list = sort(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])))

      best_runs = which(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))



      if(max(sorted_likelihood_list) > max_fit)
      {

        ### minus 1
        best_fit_pars_minus_1 = best_fit_pars
        max_fit_minus_1 = max_fit


        max_fit = max(sorted_likelihood_list)
        best_fit_pars = parameters[best_runs]

      } else if(max(sorted_likelihood_list) == max_fit)
      {
        best_fit_pars[(length(best_fit_pars) + 1 ):(length(best_fit_pars) + length(best_runs))] <- parameters[best_runs]
      }

      ### minus 1
      if(max(sorted_likelihood_list) == max_fit_minus_1)
      {
        best_fit_pars_minus_1[(length(best_fit_pars_minus_1) + 1 ):(length(best_fit_pars_minus_1) + length(best_runs))] <- parameters[best_runs]

      }
    }

    print(max_fit)
    print(c(100*i/(number_simulations/batch_size), "%"))
    # print(max(sorted_likelihood_list))

    gc()
  }

  return(list(max_fit, best_fit_pars, max_fit_minus_1, best_fit_pars_minus_1, number_simulations))

  # return(list(parameters[best_runs], likelihood_list, out, best_runs))
}

#' @export
#' @useDynLib cotonou
run_model_with_fit_multiple <- function(batch_size, number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting, solving_method = "lsoda") {




  best_fit_pars = list()
  max_fit = 1


  prep_out = c()

  ### minus 1
  best_fit_pars_minus_1 = list()
  max_fit_minus_1 = 0



  # results_list = list()
  for(i in 1:(number_simulations/batch_size))
  {

    worked = F
    while(worked == F)
    {
      # LHS to create parameter sets
        parameters <- cotonou::lhs_parameters(batch_size, set_pars = best_set, Ncat = 9, time = time,
                                              ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


      # pars = parameters[(batch_size * (i - 1) + 1):(batch_size * i)]

      # this is the slowest part - simulating model
      res = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs, solving_method = solving_method)

      if(all(lapply(lapply(res, function(x) x$prev[,1]), length) == length(time)))
        worked = T

    }





    if(all(lapply(lapply(res, function(x) x$prev[,1]), length) == length(time)))

    {# model fitting
      likelihood_list = lapply(res, cotonou::likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)

      sorted_likelihood_list = sort(unlist(lapply(likelihood_list, function(x) x[[1]])))

      best_runs = which(unlist(lapply(likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))


      prep_fit <- unlist(lapply(likelihood_list, function(x) x[[4]]))



      if(max(sorted_likelihood_list) > max_fit)
      {
        ### minus 1
        best_fit_pars_minus_1 = best_fit_pars
        max_fit_minus_1 = max_fit


        max_fit = max(sorted_likelihood_list)
        best_fit_pars = parameters[best_runs]

        prep_out = prep_fit[best_runs]

      } else if(max(sorted_likelihood_list) == max_fit)
      {
        best_fit_pars[(length(best_fit_pars) + 1 ):(length(best_fit_pars) + length(best_runs))] <- parameters[best_runs]
        prep_out[(length(prep_out) + 1 ):(length(prep_out) + length(best_runs))] <- prep_fit[best_runs]

      }

      ### minus 1
      if(max(sorted_likelihood_list) == max_fit_minus_1)
      {
        best_fit_pars_minus_1[(length(best_fit_pars_minus_1) + 1 ):(length(best_fit_pars_minus_1) + length(best_runs))] <- parameters[best_runs]

      }

    }
    print(max_fit)
    print(c(100*i/(number_simulations/batch_size), "%"))

    best_fit_pars_test <<- best_fit_pars


    gc()
  }

  return(list(max_fit, best_fit_pars, max_fit_minus_1, best_fit_pars_minus_1, prep_out))

  # return(list(parameters[best_runs], likelihood_list, out, best_runs))
}




#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster_pars_done <- function(parameters, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {


  # this is the slowest part - simulating model
  # res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)
  res = parallel::parLapply(NULL, parameters, return_outputs, main_model, time = time, outputs = outputs)

  # model fitting
  # likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points)
  likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)

  # sorted_likelihood_list = sort(unlist(lapply(likelihood_list, function(x) x[[1]])))
  sorted_likelihood_list = sort(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])))

  # best_runs = which(unlist(lapply(likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))
  best_runs = which(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))

  # out <- res[best_runs]

  # return(list(parameters[best_runs], likelihood_list, out, best_runs))
  return(list(max(sorted_likelihood_list), parameters[best_runs]))


}


#' @export
#' @useDynLib cotonou
run_model_with_fit_for_correlations_with_sets_already <- function(parameters, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {


  res = lapply(parameters, function(x) {return_outputs(x, gen = main_model, time = time, outputs = outputs)})



  # likelihood_list = unlist(lapply(res, likelihood_rough, time = time, prev_points = prev_points))
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}


#' @export
#' @useDynLib cotonou
run_model_with_fit_for_correlations <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {

  counter = 0

  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)
  # end of parameters --------------------------------------------------------------


  res = lapply(parameters, function(x) {

    counter <<- counter + 1
    if (counter %% 100 == 0)
      cat(paste("simulation:", counter))

    return_outputs(x, gen = main_model, time = time, outputs = outputs)})



  # likelihood_list = unlist(lapply(res, likelihood_rough, time = time, prev_points = prev_points))
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}

#' @export
#' @useDynLib cotonou
run_model_with_fit_for_correlations_cluster <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting = PrEP_fitting) {


  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)
  # end of parameters --------------------------------------------------------------


  res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)



  # likelihood_list = unlist(lapply(res, likelihood_rough, time = time, prev_points = prev_points))
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points, PrEP_fitting = PrEP_fitting)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}


#' @export
#' @useDynLib cotonou
likelihood_lazymcmc <- function(x, time, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points, PrEP_fitting) {




  the_N = data.frame(time, x$N[,1], rowSums(x$N[,c(5, 6, 8)]), rowSums(x$N[,c(1, 2, 3, 4, 7)]))

  the_HIV_pos = data.frame(time, x$HIV_positive[,1], rowSums(x$HIV_positive[,c(5, 6, 8)]), rowSums(x$HIV_positive[,c(1, 2, 3, 4, 7)]))

  the_On_ART = data.frame(time, x$HIV_positive_On_ART[,1], rowSums(x$HIV_positive_On_ART[,c(5, 6, 8)]), rowSums(x$HIV_positive_On_ART[,c(1, 2, 3, 4, 7)]))



  the_2012_inc_FSW = x$lambda_sum_0[which(time == 2013),1]


  names(the_N) = c("time", "Pro FSW", "Men", "Women")
  names(the_HIV_pos) = c("time", "Pro FSW", "Men", "Women")

  names(the_On_ART) = c("time", "Pro FSW", "Men", "Women")

  the_prev = data.frame(time, x$prev_FSW, x$prev_LowFSW, x$prev_client, x$prev_women, x$prev_men)
  names(the_prev) = c("time", "Pro FSW", "Low-level FSW", "Clients", "Women", "Men")

  the_frac_N = data.frame(time, x$frac_N[,c(1, 5, 7, 8)], x$frac_N[,1] + x$frac_N[,2], x$frac_N[,2]/ x$frac_N[,1])
  names(the_frac_N) = c("time", "Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio")

  ##
  frac_count <- 0

  if(all(!is.na(unlist(the_frac_N)))) {
    for (i in 1:length(frac_N_discard_points[,1]))
    {
      if(all(the_frac_N[, as.character(frac_N_discard_points[i, "variable"])] < frac_N_discard_points[i, "max"]) &&
         all(the_frac_N[,  as.character(frac_N_discard_points[i, "variable"])] > frac_N_discard_points[i, "min"])) {
        frac_count <- frac_count + 1
      }
    }
  }

  # if the run doesn't fit within the Ntot CIs, then set the frac_count to 0 so the run doesn't pass
  if(!(all(!is.na(x$Ntot)) &&
       x$Ntot[which(time == 1992)] < Ntot_data_points[1, "upper"] && x$Ntot[which(time == 1992)] > Ntot_data_points[1, "lower"] &&
       x$Ntot[which(time == 2002)] < Ntot_data_points[2, "upper"] && x$Ntot[which(time == 2002)] > Ntot_data_points[2, "lower"] &&
       x$Ntot[which(time == 2013)] < Ntot_data_points[3, "upper"] && x$Ntot[which(time == 2013)] > Ntot_data_points[3, "lower"])) {
    frac_count <- 0
    message = "Ntot"
  }

  # if(frac_count != length(frac_N_discard_points[,1]))



  lik = -Inf
  # print(paste0("frac_count ",frac_count, "lik ", lik))
# print(the_N)

  if(frac_count == length(frac_N_discard_points[,1])) {

    lik <- 0


    # prevalence
    for(i in 1:length(prev_points[,1]))
    {
      HIV_pos = subset(the_HIV_pos, time == prev_points[i, "time"], select = as.character(prev_points[i, "variable"]))
      N = subset(the_N, time == prev_points[i, "time"], select = as.character(prev_points[i, "variable"]))

      # point = subset(the_prev, time == prev_points[i, "time"], select = as.character(prev_points[i, "variable"]))
      # point = the_prev[the_prev$time == prev_points[i, "time"], as.character(prev_points[i, "variable"])]
      if(!is.na(HIV_pos)) {


        lik = lik + dbinom(x = as.numeric(prev_points[i, "x"]), size = as.numeric(prev_points[i, "N"]), prob = as.numeric(HIV_pos)/as.numeric(N), log = T)
        # lik = lik + dbinom(x = prev_points[i, "x"], size = round(as.numeric(N)), prob = as.numeric(HIV_pos/N), log = T)
        # lik = lik + dbinom(x = (prev_points[i, "value"]*as.numeric(N)/100), size = as.numeric(N), prob = prev_points[i, "value"]/100, log = T)

      }
    }


    # ART_data_points_FSW = ART_data_points[ART_data_points$variable == "Pro FSW",]
    # fitting to ART cov
    if(all(!is.na(x$ART_coverage_FSW))){
      for(i in 1:length(ART_data_points[,1]))
      {


        the_time = ART_data_points[i, "time"]
        N = subset(the_N, time == the_time, select = as.character(ART_data_points[i, "variable"]))
        On_ART = subset(the_On_ART, time == the_time, select = as.character(ART_data_points[i, "variable"]))

        # print(lik)
        lik = lik + dbinom(x = ART_data_points[i, "x"], size = round(as.numeric(N)), prob = (as.numeric(On_ART)/as.numeric(N)), log = T)


      }
    }

    lik = lik + dbinom(x = as.numeric(6), size = as.numeric(425), prob = the_2012_inc_FSW, log = T)








  }

  # print(paste(frac_count, "frac"))

  return(lik)


}

#' @export
#' @useDynLib cotonou
lazymcmc_for_cluster <- function(good_previous_fit_variedpars, ranges, best_set, time,
                                 par_seq, condom_seq, groups_seq, years_seq, outputs,
                                 prev_points,
                                 frac_N_discard_points,
                                 Ntot_data_points,
                                 ART_data_points,
                                 PrEP_fitting,
                                 iterations,
                                 popt,
                                 opt_freq,
                                 thin,
                                 burnin,
                                 adaptive_period,
                                 save_block) {

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

  # create_lik <- function(parTab, ranges, best_set, time,
  #                        par_seq, condom_seq, groups_seq, years_seq, outputs,
  #                        prev_points,
  #                        frac_N_discard_points,
  #                        Ntot_data_points,
  #                        ART_data_points,
  #                        PrEP_fitting,
  #                        PRIOR_FUNC,...)
  create_lik <- function(parTab, data, PRIOR_FUNC,...)


  {


    par_names <- rownames(ranges)

    # print(ranges)

    ## using the `...` bit.
    likelihood_func <- function(pars){

      # print(ranges)

      ranges_new = cbind(pars, pars)
      rownames(ranges_new) = par_names

      parameters = cotonou::lhs_parameters(1, set_pars = best_set, Ncat = 9, time = time,
                                           ranges = ranges_new, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)

      res = cotonou::return_outputs(parameters[[1]], cotonou::main_model, time = time, outputs = outputs)

      lik = cotonou::likelihood_lazymcmc(res, time = time,
                                         prev_points = prev_points,
                                         frac_N_discard_points = frac_N_discard_points,
                                         Ntot_data_points = Ntot_data_points,
                                         ART_data_points = ART_data_points,
                                         PrEP_fitting = PrEP_fitting)

      return(lik)
    }
    return(likelihood_func)
  }

  f <- create_lik(parTab = parTab, data = NULL, PRIOR_FUNC = NULL, ranges, best_set, time,
                  par_seq, condom_seq, groups_seq, years_seq, outputs,
                  prev_points,
                  frac_N_discard_points,
                  Ntot_data_points,
                  ART_data_points,
                  PrEP_fitting)

  if(is.numeric(f(parTab$values))) {
    mcmcPars <- c("iterations"=iterations,popt=popt,opt_freq=opt_freq,thin=thin,burnin=burnin,adaptive_period=adaptive_period,save_block=save_block)
    res <- lazymcmc::run_MCMC(parTab,
                              mcmcPars = mcmcPars,
                              filename = "test",
                              CREATE_POSTERIOR_FUNC = create_lik,
                              mvrPars = NULL, PRIOR_FUNC = NULL, OPT_TUNING = 0.1,
                              ranges = ranges,
                              best_set = best_set, time = time,
                              par_seq = par_seq, condom_seq = condom_seq,
                              groups_seq = groups_seq, years_seq = years_seq, outputs = outputs,
                              prev_points = prev_points,
                              frac_N_discard_points = frac_N_discard_points,
                              Ntot_data_points = Ntot_data_points,
                              ART_data_points = ART_data_points,
                              PrEP_fitting = PrEP_fitting)

  } else {
    res <- "failed"
  }

  return(res)



}
