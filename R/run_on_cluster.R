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
return_outputs <- function(p, gen, time, outputs) {
  mod <- gen(user = p)
  all_results <- mod$transform_variables(mod$run(time))
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
likelihood_rough <- function(x, time, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {
  the_prev = data.frame(time, x$prev_FSW, x$prev_LowFSW, x$prev_client, x$prev_women, x$prev_men)
  names(the_prev) = c("time", "Pro FSW", "Low-level FSW", "Clients", "Women", "Men")

  the_frac_N = data.frame(time, x$frac_N[,c(1, 5, 7, 8)], x$frac_N[,1] + x$frac_N[,2], x$frac_N[,2]/ x$frac_N[,1])
  names(the_frac_N) = c("time", "Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio")


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

    ART_ratio = x$Women_on_ART/x$Men_on_ART
    if(all(ART_ratio[!is.na(ART_ratio)] > 1 & ART_ratio[!is.na(ART_ratio)] < 2))
      likelihood_count <- likelihood_count + 1

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



  return (list(likelihood_count, prev_fits, message))
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
run_model_with_fit <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {


  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


  # this is the slowest part - can I turn this into
  res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)



  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)

  sorted_likelihood_list = sort(unlist(lapply(likelihood_list, function(x) x[[1]])))

  best_runs = which(unlist(lapply(likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))

  out <- res[best_runs]

  message_list <- unlist(lapply(likelihood_list, function(x) x[[3]]))

  return(list(parameters[best_runs], likelihood_list, out, best_runs, message_list))


}

#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {


  # LHS to create parameter sets
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


  # this is the slowest part - simulating model
  res = parallel::parLapply(NULL, parameters, return_outputs, main_model, time = time, outputs = outputs)



  # model fitting
  likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)

  sorted_likelihood_list = sort(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])))

  best_runs = which(unlist(parallel::parLapply(NULL, likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))

  out <- res[best_runs]

  return(list(parameters[best_runs], likelihood_list, out, best_runs))


}







#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster_multiple <- function(batch_size, number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {




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

    # model fitting
    likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)

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
run_model_with_fit_multiple <- function(batch_size, number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {




  best_fit_pars = list()
  max_fit = 1



  ### minus 1
  best_fit_pars_minus_1 = list()
  max_fit_minus_1 = 0



  # results_list = list()
  for(i in 1:(number_simulations/batch_size))
  {
    # LHS to create parameter sets
    parameters <- cotonou::lhs_parameters(batch_size, set_pars = best_set, Ncat = 9, time = time,
                                          ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)


    # pars = parameters[(batch_size * (i - 1) + 1):(batch_size * i)]

    # this is the slowest part - simulating model
    res = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)

    # model fitting
    likelihood_list = lapply(res, cotonou::likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)

    sorted_likelihood_list = sort(unlist(lapply(likelihood_list, function(x) x[[1]])))

    best_runs = which(unlist(lapply(likelihood_list, function(x) x[[1]])) == max(sorted_likelihood_list))


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


    print(max_fit)
    print(c(100*i/(number_simulations/batch_size), "%"))

    best_fit_pars_test <<- best_fit_pars


    gc()
  }

  return(list(max_fit, best_fit_pars, max_fit_minus_1, best_fit_pars_minus_1))

  # return(list(parameters[best_runs], likelihood_list, out, best_runs))
}




#' @export
#' @useDynLib cotonou
run_model_with_fit_cluster_pars_done <- function(parameters, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {


  # this is the slowest part - simulating model
  # res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)
  res = parallel::parLapply(NULL, parameters, return_outputs, main_model, time = time, outputs = outputs)

  # model fitting
  # likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points)
  likelihood_list = parallel::parLapply(NULL, res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)

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
run_model_with_fit_for_correlations_with_sets_already <- function(parameters, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {


  res = lapply(parameters, function(x) {return_outputs(x, gen = main_model, time = time, outputs = outputs)})



  # likelihood_list = unlist(lapply(res, likelihood_rough, time = time, prev_points = prev_points))
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}


#' @export
#' @useDynLib cotonou
run_model_with_fit_for_correlations <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {

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
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}

#' @export
#' @useDynLib cotonou
run_model_with_fit_for_correlations_cluster <- function(number_simulations, par_seq, condom_seq, groups_seq, years_seq, best_set, time, ranges, outputs, prev_points, frac_N_discard_points, Ntot_data_points, ART_data_points) {


  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = ranges, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)
  # end of parameters --------------------------------------------------------------


  res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)



  # likelihood_list = unlist(lapply(res, likelihood_rough, time = time, prev_points = prev_points))
  likelihood_list = lapply(res, likelihood_rough, time = time, prev_points = prev_points, frac_N_discard_points = frac_N_discard_points, Ntot_data_points = Ntot_data_points, ART_data_points = ART_data_points)


  # return(list(time, prev_points, res))

  return(list(parameters, res, likelihood_list))


}
