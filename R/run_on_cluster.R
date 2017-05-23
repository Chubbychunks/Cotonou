#' @export
#' @useDynLib cotonou
return_outputs <- function(p, gen, time, outputs) {
  mod <- gen(user = p)
  all_results <- mod$transform_variables(mod$run(time))
  #   all_results[c("prev", "c_comm_balanced", "c_noncomm_balanced", "c_comm", "c_noncomm", "epsilon")]
  return(all_results[outputs])
}

#' @export
#' @useDynLib cotonou
likelihood_rough <- function(x) {
  the_prev = data.frame(time, x$prev_FSW, x$prev_LowFSW, x$prev_client, x$prev_women, x$prev_men)
  names(the_prev) = c("time", "Pro FSW", "Low-level FSW", "Clients", "Women", "Men")

  likelihood_count <- 0

  for(i in 1:length(prev_points[,1]))
  {
    # likelihood_count <- likelihood_count +

    point = subset(the_prev, time == prev_points[i, "time"], select = as.character(prev_points[i, "variable"]))
    if(!is.na(point)) {if((point < prev_points[i, "upper"]) && (point > prev_points[i, "lower"]))
    {
      # print(prev_points[i, c("time", "variable")]);
      likelihood_count <- likelihood_count + 1
    }}
  }



  return (likelihood_count)

}

#' @export
#' @useDynLib cotonou
run_on_cluster <- function(number_simulations) {


  # parameters --------------------------------------------------------------
  parameters <- cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                        ranges = rbind(

                                          init_clientN_from_PCR = c(1,1),
                                          # NO HIV, CONSTANT POP GROWTH RATE
                                          epsilon_1985 = c(0.08, 0.08),
                                          epsilon_1992 = c(0.08, 0.08),
                                          epsilon_2002 = c(0.08, 0.08),
                                          epsilon_2013 = c(0.08, 0.08),
                                          epsilon_2016 = c(0.08, 0.08),

                                          # epsilon_1985 = c(0.059, 0.059),
                                          # epsilon_1992 = c(0.059, 0.059),
                                          # epsilon_2002 = c(0.059, 0.059),
                                          # epsilon_2013 = c(0.059, 0.059),
                                          # epsilon_2016 = c(0.059, 0.059),

                                          # muF = c(0.05, 0.05),
                                          # muM = c(0.06, 0.06),

                                          muF = c(0.0295, 0.0295),
                                          muM = c(0.0315, 0.0315),

                                          # betaMtoF_noncomm = c(0.00144, 0.00626),

                                          betaMtoF_noncomm = c(0, 0),
                                          frac_women_ProFSW = c(0.004, 0.004),
                                          # frac_women_ProFSW = c(0.0024, 0.0067),
                                          # frac_women_LowFSW = c(0.0024, 0.0067),
                                          # frac_women_exFSW = c(0.0024, 0.0067),
                                          frac_men_client = c(0.6, 0.6),
                                          # frac_women_virgin = 0.1,
                                          # frac_men_virgin = 0.1



                                          RR_beta_GUD = c(1.43, 19.58),
                                          RR_beta_FtM = c(0.5, 2),

                                          c_comm_1993_ProFSW = c(1000, 1800),
                                          c_comm_2005_ProFSW = c(250, 600),
                                          c_comm_1998_Client = c(7, 12),
                                          c_comm_2015_Client = c(6, 12),

                                          c_noncomm_1998_Client = c(1, 3),
                                          c_noncomm_2015_Client = c(2, 6),

                                          who_believe_comm = c(0, 1),

                                          rate_leave_pro_FSW = c(0.4347826, 0.4347826),
                                          rate_leave_low_FSW = c(0.4347826, 0.4347826),
                                          rate_leave_client = c(0.5, 0.5),
                                          rate_enter_sexual_pop = c(0.3571429, 0.3571429),




                                          fc_y_comm_1993_ProFSW_Client = c(0.535, 0.687),
                                          fc_y_comm_2002_ProFSW_Client = c(0.872, 0.933),
                                          fc_y_comm_1998_ProFSW_Client = c(0.872, 0.933), # fake

                                          fc_y_noncomm_1985_ProFSW_Client = c(0.27, 0.43),
                                          fc_y_noncomm_2016_ProFSW_Client = c(0.27, 0.43),

                                          fc_y_noncomm_1998_GPM_GPF = c(0.0326087, 0.241404781),
                                          fc_y_noncomm_2016_GPM_GPF = c(0.0326087, 0.251404781)


                                        ))
  # end of parameters --------------------------------------------------------------


  # res = lapply(parameters, f, main_model, time = seq(1986, 2030, 1))
  res = lapply(parameters, return_outputs, main_model, time = time, outputs = outputs)



  # prev_points -------------------------------------------------------------
  prev_points = data.frame(time = c(1986, 1987, 1988, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015,
                                    1998, 2002, 2005, 2008, 2012, 2015,
                                    1998, 2008, 2011,
                                    1998, 2008, 2011,
                                    2012, 2015),
                           variable = c(rep("Pro FSW", 11),
                                        rep("Clients", 6),
                                        rep("Women", 3),
                                        rep("Men", 3),
                                        rep("Low-level FSW", 2)),
                           value = c(3.3, 8.2, 19.2, 53.3, 48.7, 40.6, 38.9, 34.8, 29.3, 27.4, 18.7,
                                     100*0.084, 9, 6.9, 5.8, 100*0.028, 100*0.016,
                                     100*0.035, 100*0.04, 2.2,
                                     100*0.033, 100*0.02, 1.6,
                                     100*0.167, 100*0.065),
                           lower = c(3.3, 8.2, 19.2, 48.02, 43.02, 36.58, 31.97, 30.42, 24.93, 23.01, 15.71,
                                     100*0.05898524, 100*0.068218538, 100*0.04293149, 100*0.034772131, 100*0.012660836, 100*0.006039259,
                                     100*0.024181624, 100*0.030073668, 100*0.012980254,
                                     100*0.022857312, 100*0.012427931, 100*0.007517563,
                                     100*0.091838441, 100*0.026704897),
                           upper = c(3.3, 8.2, 19.2, 58.48, 54.42, 44.67, 46.27, 39.38, 33.88, 32.23, 22.01,
                                     100*0.11561791, 100*0.115608811, 100*0.105215792, 100*0.090216628, 100*0.051602442, 100*0.035338436,
                                     100*0.047726245, 100*0.052817187, 100*0.035296286,
                                     100*0.047183668, 100*0.029774338, 100*0.028546718,
                                     100*0.268127672, 100*0.130153465))
  prev_points_all = prev_points
  prev_points = prev_points[-c(1,2,3),]


  # prev_points -------------------------------------------------------------


  # best runs etc -----------------------------------------------------------

  likelihood_list = unlist(lapply(res, likelihood_rough))
  sorted_likelihood_list = sort(likelihood_list)

  # table(sorted_likelihood_list)

  best_runs = which(unlist(lapply(res, likelihood_rough)) == max(sorted_likelihood_list))

  out <- res[best_runs]

  return(out)
}
