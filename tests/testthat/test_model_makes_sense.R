context("model makes sense?")


#  ___ ____ ____  _   _ _____ ____
# |_ _/ ___/ ___|| | | | ____/ ___|
#  | |\___ \___ \| | | |  _| \___ \
#  | | ___) |__) | |_| | |___ ___) |
# |___|____/____/ \___/|_____|____/
#

# sometimes it turns out that betaMtF_noncomm is 0 in parameters, so the test doesn't work
# movement onto prep makes the replaceDeaths thing not work...


# not sure if want lhs_parameters or generate_parameters

test_that("example", {
  expect_equal(1 + 1, 2, tolerance = 1e-6)
  expect_true(1 + 1 == 2)
  expect_error(sqrt(lm), "non-numeric")
  # expect_identical(3, sqrt(3)^2)
})

#result = run_model(parameters, main_model, time, output_vars = c("Ntot", "prev_client"))

#GENERAL TESTS

#Ncat works?
# test_that("Ncat", {
#   for (Ncat in c(2, 10))
#   {
#     parameters <- lhs_parameters(1, Ncat = Ncat, best_set = best_set_default)
#     result = run_model(parameters, main_model, time)
#     expect_equal(ncol(result$S0), Ncat)
#   }
# })

# ALL COMPARTMENTS ARE POSITIVE
test_that("all compartments positive", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("S[0-9]|I[0-9]", names(result)))]
  expect_true(all(unlist(xx) > -10^-3))
})

# CUMULATIVE INFECTIONS ALWAYS POSITIVE

test_that("cumulative infections", {
  parameters <- lhs_parameters(1, S1b_init = rep_len(101, 9), S1c_init = rep_len(100, 9), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  expect_true(all(diff(xx[[1]][,1]) >= 0))
  expect_true(all(diff(xx[[1]][,2]) >= 0))

})

# NO SEEDING OF EPIDEMIC
###################################################################################################################################
###################################################################################################################################

# no infected, no incidence?
test_that("no incidence, no cumulative infections", {

  pars = list(prev_init_FSW = 0, prev_init_rest = 0)

  parameters <- lhs_parameters(1, forced_pars = modifyList(pars, list(time = time_default)), S1a_init = rep(100,9), S1b_init = rep(100,9), S1c_init = rep(100,9), S1b_init = rep_len(101, 9), S1c_init = rep_len(100, 9), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default[-which(rownames(ranges_default) %in% names(pars)),])
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9]", names(result)))]
  expect_true(all(unlist(xx) == 0))
  expect_equal(ncol(result$S0), 9)

  xx <- result[c(grep("cumuInf", names(result)))]
  expect_true(all(unlist(xx) == 0))

})


# GROWTH RATE AND DEMOGRAPHY
###################################################################################################################################
###################################################################################################################################


# add test that there are the correct amounts of omega, mu etc


test_that("omega adds to 1", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)

  expect_equal(sum(parameters[[1]]$omega), 1)
})

test_that("omega keeps consistent population?", {

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default,
                               forced_pars = list(omega = c(0.01, 0.02, 0.3, 0.1, 0.12, 0.25, 0.1, 0.1, 0), beta_comm = c(0,0,0,0,0,0,0,0,0), beta_noncomm = c(0,0,0,0,0,0,0,0,0),
                                                  S0_init = c(100*0.01, 100*0.02, 100*0.3, 100*0.1, 100*0.12, 100*0.25, 100*0.1, 100*0.1, 100*0),I01_init = c(100*0.01, 100*0.02, 100*0.3, 100*0.1, 100*0.12, 100*0.25, 100*0.1, 100*0.1, 100*0),
                                                  time = time_default, replaceDeaths = 1, movement = 0))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[grep("frac_N", names(result))] # grepping all the Ss and Is


  expect_true(all(abs(diff(xx$frac_N))<10^-12))
  expect_equal(as.numeric(xx$frac_N[1,]), as.numeric(xx$frac_N[2,]))
})


test_that("omega keeps consistent population even with HIV? recruitment to PrEP is 0", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default, forced_pars = list(omega = c(0.01, 0.02, 0.3, 0.1, 0.12, 0.25, 0.1, 0.1, 0),
                                                                                                                                                                                                                                                            S0_init = c(100*0.01, 100*0.02, 100*0.3, 100*0.1, 100*0.12, 100*0.25, 100*0.1, 100*0.1, 100*0),
                                                                                                                                                                                                                                                            I01_init = c(100*0.01, 100*0.02, 100*0.3, 100*0.1, 100*0.12, 100*0.25, 100*0.1, 100*0.1, 100*0),
                                                                                                                                                                                                                                                            time = time_default, replaceDeaths = 1, movement = 0, eP1a = rep(0, 9), eP1b = rep(0, 9), eP1c = rep(0, 9)), set_null = list("zetaa_y", "zetab_y", "zetac_y"))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[grep("frac_N", names(result))] # grepping all the Ss and Is


  expect_true(all(abs(diff(xx$frac_N))<10^-12))
  expect_equal(as.numeric(xx$frac_N[1,]), as.numeric(xx$frac_N[2,]))
})

test_that("growth rate zero if replacing deaths and if no zeta - otherwise it slightly doesn't work :(", {
  parameters <- lhs_parameters(1, replaceDeaths = 1, movement = 0, epsilon_y = c(0,0,0,0,0), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default, set_null = list("zetaa_y", "zetab_y", "zetac_y"))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]



  xx <- result[grep("^[SI]", names(result))] # grepping all the Ss and Is
  N <- rowSums(do.call(cbind, xx))

  # are all increments in N equal to 0?
  expect_true(all(abs(diff(N)) < 10^-2))
})

test_that("growth rate increases", {
  parameters <- lhs_parameters(1, epsilon_y = c(0.1,0.1,0.1,0.1,0.1), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[grep("^[SI]", names(result))] # grepping all the Ss and Is
  N <- rowSums(do.call(cbind, xx))

  # test 2: are all increments in N positive AND are the increments getting bigger?
  # expect_true(all(diff(N) > 0) && all(diff(diff(N)) > 0))
  expect_true(all(diff(N) > 0))
})

test_that("growth rate decreases", {

  parameters <- lhs_parameters(1, epsilon_y = c(-0.01,-0.01,-0.01,-0.01,-0.01), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[grep("^[SI]", names(result))] # grepping all the Ss and Is
  N <- rowSums(do.call(cbind, xx))

  # test 2: are all increments in N negative AND are the increments getting bigger?
  # expect_true(all(diff(N) < 0) && all(diff(diff(N)) > 0))

  # just all increments in N negative
  expect_true(all(diff(N) < 0))

})





# LAMBDAS
# if prep is useless, then cumulative infections should be equal no matter what prep adherence is
test_that("useless prep", {

  parameters <- lhs_parameters(1, I11_init = rep(1000, 9), I01_init = rep(1000, 9), zetaa_y = matrix(rep(0, 45), ncol = 9), zetab_y = matrix(rep(0, 45), ncol = 9), zetac_y = matrix(rep(0, 45), ncol = 9),
                               zetaa_t = c(1985, 2013, 2015, 2016, 2020), zetab_t = c(1985, 2013, 2015, 2016, 2020), zetac_t = c(1985, 2013, 2015, 2016, 2020),
                               eP0 = rep(0, 9), eP1a = rep(0, 9), eP1b = rep(0, 9), eP1c = rep(0, 9), eP1d = rep(0, 9), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[grep("cumuInf", names(result))] # grepping all the Ss
  N1 <- rowSums(do.call(cbind, xx))

  parameters2 <- modifyList(parameters, list(zetaa_y = matrix(rep(c(0, 0, 0, 0, 0.1, 0.1, 0, 0, 0), 5), byrow = F, ncol = 9), zetab_y = matrix(rep(c(0, 0, 0, 0, 0.1, 0.1, 0, 0, 0), 5), byrow = F, ncol = 9), zetac_y = matrix(rep(c(0, 0, 0, 0, 0.1, 0.1, 0, 0, 0), 5), byrow = F, ncol = 9)))

  result2 = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx2 <- result2[grep("cumuInf", names(result2))] # grepping all the Ss
  N2 <- rowSums(do.call(cbind, xx2))

  # NOTE FOR THIS TEST THAT IT ONLY WORKS IF PREP UPTAKE DOESNT HAPPEN EARLY!!!!! IT AFFECTS HOW THE POPULATION GROWS...

  #   which(unlist(parameters) - unlist(parameters2) != 0)

  expect_true(all(abs(N1 - N2) < 10^-2))
})


test_that("useful prep", {

  parameters <- lhs_parameters(1, PrEPOnOff = 1, I11_init = rep(1000, 9), I01_init = rep(1000, 9), zetaa_y = matrix(rep(0.1, 45), ncol = 9), zetab_y = matrix(rep(0.1, 45), ncol = 9), zetac_y = matrix(rep(0.1, 45), ncol = 9),
                               zetaa_t = c(1985, 2013, 2015, 2016, 2020), zetab_t = c(1985, 2013, 2015, 2016, 2020), zetac_t = c(1985, 2013, 2015, 2016, 2020),
                               eP0 = rep(0, 9), eP1a = rep(0.1, 9), eP1b = rep(0.1, 9), eP1c = rep(0.1, 9), eP1d = rep(0, 9), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result1 = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  newpars <- lhs_parameters(1, set_null = "eP1a", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$eP1a

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP1a = newpars, eP1b = newpars, eP1c = newpars)))

  result2 = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  expect_true(sum(unlist(result1[c(grep("I[0-9]", names(result1)))])) < sum(unlist(result2[c(grep("I[0-9]", names(result2)))])))

})



# FORCE OF INFECTION SET TO ZERO
###################################################################################################################################
###################################################################################################################################


# BETA
# only FSW infected, do males get infected?
test_that("beta 1", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default[which(rownames(ranges_default) != "RR_beta_FtM"),], forced_pars = list(RR_beta_FtM = 0, time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})

# only GPF infected, do males get infected?
test_that("beta 2", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default[which(rownames(ranges_default) != "RR_beta_FtM"),], forced_pars = list(RR_beta_FtM = 0, time = time_default, I11_init = c(0, 0, 1000, 0, 0, 0, 0, 0, 0), I01_init = c(0, 0, 1000, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})

# only clients infected, do females get infected?
test_that("beta 3", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default[which(rownames(ranges_default) != "betaMtoF_baseline"),],
                               forced_pars = list(betaMtoF_comm = 0, betaMtoF_noncomm = 0, betaMtoF_baseline = 0, movement = 0,
                                                  time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4]))==0)
})


# R
# if R is 0, no infections
test_that("R 1", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), R = 0, infect_ART = 0, infect_acute = 0, infect_AIDS = 0))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})




# n
# only group 1 infected, does group 2 get infected if n = 0?
test_that("n 1", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(betaMtoF_noncomm = 0.001,

                                                  # n_comm = matrix(c(0), nrow = 9, ncol = 9), n_noncomm = matrix(c(0), nrow = 9, ncol = 9),
                                                  n_y_noncomm_1985 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2002 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2015 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2016 = matrix(0, ncol=9, nrow=9),
                                                  n_y_comm_1985 = matrix(0, ncol=9, nrow=9), n_y_comm_2002 = matrix(0, ncol=9, nrow=9), n_y_comm_2015 = matrix(0, ncol=9, nrow=9), n_y_comm_2016 = matrix(0, ncol=9, nrow=9),


                                                  time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})

# n
# only group 1 infected, does group 2 get infected if n > 0?
test_that("n 2", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(betaMtoF_noncomm = 0.001, n_y_comm_2002 = matrix(c(0.1), nrow = 9, ncol = 9), n_y_noncomm_2002 = matrix(c(0.1), nrow = 9, ncol = 9),
                                                  time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)
})


# only group 2 infected, does group 1 get infected n=0?
test_that("n 3", {

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(betaMtoF_noncomm = 0.001,
                                                  # n_y_comm_2002 = matrix(c(0), nrow = 9, ncol = 9), n_y_noncomm_2002 = matrix(c(0), nrow = 9, ncol = 9),
                                                  n_y_noncomm_1985 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2002 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2015 = matrix(0, ncol=9, nrow=9), n_y_noncomm_2016 = matrix(0, ncol=9, nrow=9),
                                                  n_y_comm_1985 = matrix(0, ncol=9, nrow=9), n_y_comm_2002 = matrix(0, ncol=9, nrow=9), n_y_comm_2015 = matrix(0, ncol=9, nrow=9), n_y_comm_2016 = matrix(0, ncol=9, nrow=9),

                                                  time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4]))==0)
})

# only group 2 infected, does group 1 get infected n>0?
test_that("n 3", {

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(betaMtoF_noncomm = 0.001, n_y_comm_2002 = matrix(c(0.1), nrow = 9, ncol = 9), n_y_noncomm_2002 = matrix(c(0), nrow = 9, ncol = 9),
                                                  time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4])) > 0)
})



# spent ages trying to do condom but doesn't work...



# fc and ec
# only group 1 infected, does group 2 get infected?
# condom efficacy is 1 and frequency of condom use is 1 - no infections
test_that("fc ec 1", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 fc_y_comm_1985 = matrix(1, ncol=9,nrow=9), fc_y_comm_1993 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_1995 = matrix(1, ncol=9,nrow=9), fc_y_comm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2002 = matrix(1, ncol=9,nrow=9), fc_y_comm_2005 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2008 = matrix(1, ncol=9,nrow=9), fc_y_comm_2012 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2015 = matrix(1, ncol=9,nrow=9), fc_y_comm_2016 = matrix(1, ncol=9,nrow=9),
                                 fc_t_comm = c(1985, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015, 2016, 2020),
                                 fc_y_noncomm_1985 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1993 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2002 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2008 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2011 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2015 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2016 = matrix(1, ncol=9,nrow=9), fc_t_noncomm = c(1985, 1993, 1998, 2002, 2008, 2011, 2015, 2016, 2020),
                                 ec = rep(1, 9), ignore_ranges_fc_c = 1,
                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})

# condom efficacy is NOT 1 and frequency of condom use is 1 - some infections
test_that("fc ec 1b", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 fc_y_comm_1985 = matrix(1, ncol=9,nrow=9), fc_y_comm_1993 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_1995 = matrix(1, ncol=9,nrow=9), fc_y_comm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2002 = matrix(1, ncol=9,nrow=9), fc_y_comm_2005 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2008 = matrix(1, ncol=9,nrow=9), fc_y_comm_2012 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2015 = matrix(1, ncol=9,nrow=9), fc_y_comm_2016 = matrix(1, ncol=9,nrow=9),
                                 fc_t_comm = c(1985, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015, 2016, 2020),
                                 fc_y_noncomm_1985 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1993 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2002 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2008 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2011 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2015 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2016 = matrix(1, ncol=9,nrow=9), fc_t_noncomm = c(1985, 1993, 1998, 2002, 2008, 2011, 2015, 2016, 2020),
                                 ec = rep(0.99, 9), ignore_ranges_fc_c = 1,
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),
                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)
})

# condom efficacy is 1 and frequency of condom use is NOT 1 - some infections
test_that("fc ec 1c", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(betaMtF_noncomm = 0.001,
                                 fc_y_comm_1985 = matrix(1, ncol=9,nrow=9), fc_y_comm_1993 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_1995 = matrix(1, ncol=9,nrow=9), fc_y_comm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2002 = matrix(1, ncol=9,nrow=9), fc_y_comm_2005 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2008 = matrix(1, ncol=9,nrow=9), fc_y_comm_2012 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2015 = matrix(1, ncol=9,nrow=9), fc_y_comm_2016 = matrix(1, ncol=9,nrow=9),
                                 fc_t_comm = c(1985, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015, 2016, 2020),
                                 fc_y_noncomm_1985 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1993 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2008 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2011 = matrix(0.9, ncol=9,nrow=9), fc_y_noncomm_2015 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2016 = matrix(1, ncol=9,nrow=9), fc_t_noncomm = c(1985, 1993, 1998, 2002, 2008, 2011, 2015, 2016, 2020),
                                 ec = rep(1, 9), ignore_ranges_fc_c = 1,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)
})


# only group 2 infected, does group 1 get infected?
test_that("fc ec 2", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 fc_y_comm_1985 = matrix(1, ncol=9,nrow=9), fc_y_comm_1993 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_1995 = matrix(1, ncol=9,nrow=9), fc_y_comm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2002 = matrix(1, ncol=9,nrow=9), fc_y_comm_2005 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2008 = matrix(1, ncol=9,nrow=9), fc_y_comm_2012 = matrix(1, ncol=9,nrow=9),
                                 fc_y_comm_2015 = matrix(1, ncol=9,nrow=9), fc_y_comm_2016 = matrix(1, ncol=9,nrow=9),
                                 fc_t_comm = c(1985, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015, 2016, 2020),
                                 fc_y_noncomm_1985 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1993 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_1998 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2002 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2008 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2011 = matrix(1, ncol=9,nrow=9), fc_y_noncomm_2015 = matrix(1, ncol=9,nrow=9),
                                 fc_y_noncomm_2016 = matrix(1, ncol=9,nrow=9), fc_t_noncomm = c(1985, 1993, 1998, 2002, 2008, 2011, 2015, 2016, 2020),
                                 ec = rep(1, 9), ignore_ranges_fc_c = 1,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),

                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4]))==0)
})







# not sure if need... was trying to do fP and eP and had copied and pasted above for the condom.
#
# fP and eP
# only group 1 infected, does group 2 get infected?
# prep efficacy is 1 and frequency of prep use is 1 - no infections
# NOTE KAPPA AND ZETA MUST BE 0!
test_that("fP eP 1", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6]))==0)
})



#
# PREP efficacy is NOT 1 and frequency of PREP use is 1 - some infections
test_that("fP eP 1b", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1, betaMtF_noncomm = 0.001,
                                 eP0 = rep(0.9, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1, betaMtF_noncomm = 0.001,
                                 eP0 = rep(1, 9), eP1a = rep(0.9, 9), eP1b = rep(1, 9), eP1c = rep(1, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1, betaMtF_noncomm = 0.001,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(0.9, 9), eP1c = rep(1, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1, betaMtF_noncomm = 0.001,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(0.9, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)


})

# PrEP efficacy is 1 and frequency of prep use is NOT 1 - some infections
test_that("fP eP 1c", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(1, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(0.99, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,


                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0), I01_init = c(1000, 0, 0, 0, 0, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,5], xx$I11[,5], xx$I01[,6], xx$I11[,6])) > 0)
})
#
#
# only group 2 infected, does group 1 get infected?
test_that("fP eP 2", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(1, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4]))==0)
})

test_that("fP eP 2b", {
  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(0.9, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(1, 9), eP1d = rep(1, 9),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))



  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]

  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(0.9, 9), eP1b = rep(1, 9), eP1c = rep(1, 9), eP1d = rep(1, 9), zetaa_y = matrix(0.1, ncol = 9, nrow = 5),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),

                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(0.9, 9), eP1c = rep(1, 9), eP1d = rep(1, 9), zetab_y = matrix(0.1, ncol = 9, nrow = 5),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),

                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4])) > 0)

  parameters <- lhs_parameters(1, Ncat = 9, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(PrEPOnOff = 1,
                                 eP0 = rep(1, 9), eP1a = rep(1, 9), eP1b = rep(1, 9), eP1c = rep(0.9, 9), eP1d = rep(1, 9), zetac_y = matrix(0.1, ncol = 9, nrow = 5),
                                 fP_y_comm = matrix(1, nrow = 5, ncol = 9), fP_y_noncomm = matrix(1, nrow = 5, ncol = 9),
                                 fP_t_comm = c(1985, 2014, 2015, 2016, 2030),  fP_t_noncomm = c(1985, 2014, 2015, 2016, 2030),
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),


                                 time = time_default, I11_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0), I01_init = c(0, 0, 0, 0, 1000, 0, 0, 0, 0)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  xx <- result[c(grep("I01", names(result)), grep("I11", names(result)))]
  expect_true(sum(c(xx$I01[,1], xx$I11[,1], xx$I01[,2], xx$I11[,2], xx$I01[,3], xx$I11[,3], xx$I01[,4], xx$I11[,4])) > 0)
})



# ALL FORCES OF INFECTION POSITIVE?!



# DISEASE PROGRESSION
###################################################################################################################################
###################################################################################################################################

# Setting progression rate from acute to CD4>500 to zero
# done by setting tau[0-9]1 and gamma[0-9]1 to 0

test_that("acute to CD4>500 zero", {
  relevant_parameters = parameter_names[c(grep("gamma[0-9]1", parameter_names), grep("testing_prob_y", parameter_names))]


  parameters <- lhs_parameters(1, I11_init = rep(100, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  all_infected = result[c(grep("I[0-9]2|I[0-9]3|I[0-9]4|I[0-9]5", names(result)))]
  expect_true(all(unlist(all_infected) == 0))
})

# Setting progression rate from CD4>500 to CD4 350-500 to zero
# done by setting gamma[0-9]2 to 0

test_that("CD4>500 to CD4 350-500 zero", {
  relevant_parameters = parameter_names[c(grep("gamma[0-9]2", parameter_names))]


  parameters <- lhs_parameters(1, I11_init = rep(100, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  all_infected = result[c(grep("I[0-9]3|I[0-9]4|I[0-9]5", names(result)))]
  expect_true(all(unlist(all_infected) == 0))
})

# Setting progression rate from CD4 350-500 to CD4 200-349 to zero
# done by setting gamma[0-9]3 to 0

test_that("CD4 350-500 to CD4 200-349 zero", {
  relevant_parameters = parameter_names[c(grep("gamma[0-9]3", parameter_names))]
  parameters <- lhs_parameters(1, I11_init = rep(100, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  all_infected = result[c(grep("I[0-9]4|I[0-9]5", names(result)))]
  expect_true(all(unlist(all_infected) == 0))
})


# Setting progression rate from CD4 200-349 to CD4 <200 to zero
# done by setting gamma[0-9]4 to 0

test_that("CD4 200-349 to CD4 <200 to zero", {
  relevant_parameters = parameter_names[c(grep("gamma[0-9]4", parameter_names))]
  parameters <- lhs_parameters(1, I11_init = rep(100, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  all_infected = result[c(grep("I[0-9]5", names(result)))]
  expect_true(all(unlist(all_infected) == 0))
})

# CARE CASCADE PROGRESSION
###################################################################################################################################
###################################################################################################################################


# PREP

test_that("prep", {
  relevant_parameters = parameter_names[c(grep("zeta", parameter_names))]
  parameters <- lhs_parameters(1,



                               I11_init = rep(0, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(time = time_default,
                               prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                               sigma = c(1,1,1,1,1,1,1,1,1),
                               test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                               PrEPOnOff = 1


                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I11", names(result)), grep("S1", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(all(diff(N) != 0))
})


test_that("no prep", {
  relevant_parameters = parameter_names[c(grep("zeta", parameter_names))]
  parameters <- lhs_parameters(1,



                               I11_init = rep(0, 9), set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(time = time_default,
                                                  prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                                  sigma = c(1,1,1,1,1,1,1,1,1),
                                                  test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                                  PrEPOnOff = 0


                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I11", names(result)), grep("S1", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(all(diff(N) == 0))
})

test_that("prep increases", {

  parameters <- lhs_parameters(1, I11_init = rep(0, 9), zetaa_y = matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0), byrow = F, ncol=9, nrow = 5), zetab_y = matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0), byrow = F, ncol=9, nrow = 5), zetac_y = matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0), byrow = F, ncol=9, nrow = 5),
                               par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(time = time_default,
                                                  prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                                  sigma = c(1,1,1,1,1,1,1,1,1),
                                                  test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                                  PrEPOnOff = 1


                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I11", names(result)), grep("S1", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(!all(diff(N) == 0))
})

# NO TESTING

test_that("no testing", {
  relevant_parameters = parameter_names[c(grep("testing_prob_y", parameter_names))]

  parameters <- lhs_parameters(1, set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I2[0-9]|I3[0-9]|I4[0-9]", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(all(diff(N) == 0))
})


# NO ART

test_that("no ART", {
  relevant_parameters = parameter_names[c(grep("ART_prob_y", parameter_names))]
  parameters <- lhs_parameters(1, set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I3[0-9]|I4[0-9]", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(all(diff(N) == 0))
})


# NO DROP OUT

test_that("no drop out", {
  relevant_parameters = parameter_names[c(grep("phi", parameter_names))]
  parameters <- lhs_parameters(1, set_null = relevant_parameters, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I4[0-9]", names(result)))]
  N <- rowSums(do.call(cbind, xx))
  expect_true(all(diff(N) == 0))
})

# BALANCING
###################################################################################################################################
###################################################################################################################################

# DUNNO!
test_that("B check 0", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  expect_equal(result$B_check_comm, result$B_check_noncomm)
  expect_equal(result$B_check_comm, rep(0, length(time_default)))
})

# c_comm_balanced? contained in the above test tbh


# p makes sense?

test_that("p makes sense", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
  true_vec = c()
  for(i in 1:length(time_default))
    true_vec[i] = all(round(rowSums(result$p_comm[i,,]), 4) %in% seq(0, 1))
  expect_true(all(true_vec))

  for(i in 1:length(time_default))
    true_vec[i] = all(round(rowSums(result$p_noncomm[i,,]), 4) %in% seq(0, 1))
  expect_true(all(true_vec))

})



# CALCULATING PREVALENCE
###################################################################################################################################
###################################################################################################################################
test_that("prevalence", {


  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  all_infected = result[c(grep("I[0-9]", names(result)))]
  all = result[c(grep("I[0-9]", names(result)), grep("^S[0-9]", names(result)))]
  for(i in 1:9)
    expect_equal(result$prev[,i], 100 * rowSums(do.call(cbind, lapply(all_infected, function(x) x <- x[,i]))) / rowSums(do.call(cbind, lapply(all, function(x) x <- x[,i]))), tolerance = 1e-6)

  # this will need to be tested against overall prevalence
  #over_prevalence = rowSums(do.call(cbind, all_infected)) / rowSums(do.call(cbind, all))

  # result$prev}

})

#  OVERALL PREVALENCE IS EQUAL TO WEIGHTED AVERAGE OF ALL PREVALENCES

# CALCULATING INCIDENCE
###################################################################################################################################
###################################################################################################################################

# set all mortality to zero, set births to zero
# incidence can be calculated by:
# lambda * S

#dont understand this test

# test_that("comparing incidence", {
#   relevant_parameters = parameter_names[c(grep("gamma[0-9]4", parameter_names))]
#   parameters <- lhs_parameters(1, I11_init = c(100,100), set_null = relevant_parameters)[[1]]
#   result = run_model(parameters, main_model, time)
#
#   parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
#   result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]
#
#   all_infected = result[c(grep("I[0-9]5", names(result)))]
#   expect_true(all(unlist(all_infected) == 0))
# })

# CALCULATING FORCE OF INFECTION
###################################################################################################################################
###################################################################################################################################



# LOGICAL CHECKS FOR MODEL
###################################################################################################################################
###################################################################################################################################

# increase beta, increase overall prevalence

test_that("beta vs prevalence", {

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 betaMtoF_noncomm = 0.001, time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(beta_noncomm =  x$beta_noncomm * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

# increase R, increase overall prevalence

test_that("R vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default, R = 0.1,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(R =  x$R * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

# increase n, increase overall prevalence

test_that("n vs prevalence", {

  parameters <- lhs_parameters(1,par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,  forced_pars = list(
    time = time_default,
    n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
    n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
  ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(n_y_comm =  x$n_y_comm * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))

})

test_that("n vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(n_y_noncomm =  x$n_y_noncomm * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

# increase c, increase overall prevalence

test_that("c_comm vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(c_y_comm =  x$c_y_comm * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

test_that("c_noncomm vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(c_y_noncomm =  x$c_y_noncomm * 1.01)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})


# increase fc, decrease overall prevalence

test_that("fc vs prevalence comm", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  newpars <- lhs_parameters(1, set_null = "fc_y_comm", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$fc_y_comm

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(fc_y_comm = newpars)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

test_that("fc vs prevalence noncomm", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  newpars <- lhs_parameters(1, set_null = "fc_y_noncomm", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$fc_y_noncomm

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(fc_y_noncomm = newpars)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

# increase ec, decrease overall prevalence

test_that("ec vs prevalence", {

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,
                                 ec = rep(0.9, 9)
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(ec =  x$ec * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})



# increase fP, decrease overall prevalence

test_that("fP vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,
                                 zetaa_y = array(data=c(0.2), dim = c(5, 9)),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1

                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  newpars <- lhs_parameters(1, set_null = "fP_y_comm", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$fP_y_comm

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(fP_y_comm = newpars)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

test_that("fP vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,
                                 zetaa_y = array(data=c(0.2), dim = c(5, 9)),
                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  newpars <- lhs_parameters(1, set_null = "fP_y_noncomm", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$fP_y_noncomm

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(fP_y_noncomm = newpars)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

# increase eP, decrease overall prevalence

test_that("eP vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1,
                                 betaMtF_noncomm = 0.001, eP0 = rep(0.9, 9)
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP0 =  x$eP0 * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


  parameters <- lhs_parameters(1,zetaa_y = matrix(0.1, ncol = 9, nrow = 5), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001, eP1a = rep(0.9, 9),
                                 zetaa_y = array(data=c(0.2), dim = c(5, 9)),
                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP1a =  x$eP1a * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


  parameters <- lhs_parameters(1, zetaa_y = matrix(0.1, ncol = 9, nrow = 5), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1,betaMtF_noncomm = 0.001, eP1b = rep(0.9, 9),
                                 zetaa_y = array(data=c(0.2), dim = c(5, 9)),
                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP1b =  x$eP1b * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))

  parameters <- lhs_parameters(1, zetaa_y = matrix(0.1, ncol = 9, nrow = 5), par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001, eP1c = rep(0.9, 9),
                                 zetaa_y = array(data=c(0.2), dim = c(5, 9)),
                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP1c =  x$eP1c * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))

  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1d = rep(0.9, 9), eP0 = rep(0.9, 9),
                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  p1 <- parameters

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(eP1d =  x$eP1d * 0.5)))

  p2 <- parameters

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))

})

# unlist(p1[c(unlist(p1) == unlist(p2))])
# which( == T)






# increase prep uptake, decrease overall prevalence

test_that("zeta vs prevalence", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,

                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  newpars <- lhs_parameters(1, set_null = "sigma", par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]]$sigma

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(sigma = newpars)))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))




  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(sigma =  x$sigma * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))





  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(test_rate_prep =  x$test_rate_prep * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))



  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0, 9), kappab = rep(0, 9), kappac = rep(0, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(PrEPOnOff =  x$PrEPOnOff * 0)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


  # below is no longer true but maybe still useful remembering that prep may be disadvantageous if the testing rate of those on prep is slower and ART is good
  #   # in this test, i made gammas and taus 0 to make sure prep doesn't advtange by going to ART quicker
  #   parameters <- lhs_parameters(1, zetaa = c(0.1, 0.1), eP0 = c(0, 0), eP1a = c(0, 0), eP1b = c(0, 0), eP1c = c(0, 0), gamma01 = c(0, 0), gamma11 = c(0, 0), tau01 = c(0, 0), tau11 = c(0, 0))[[1]]
  #   result1 = run_model(parameters, main_model, time)
  #   parameters <- lhs_parameters(1, zetaa = c(0.09, 0.09), eP0 = c(0, 0), eP1a = c(0, 0), eP1b = c(0, 0), eP1c = c(0, 0), gamma01 = c(0, 0), gamma11 = c(0, 0), tau01 = c(0, 0), tau11 = c(0, 0))[[1]]
  #   result2 = run_model(parameters, main_model, time)
  #   expect_equal(result1$cumuInf[length(time)], result2$cumuInf[length(time)])

})





# increase ART uptake, decrease overall prevalence

test_that("ART vs prevalence", {

  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0, 9), kappab = rep(0, 9), kappac = rep(0, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(ART_prob_y =  x$ART_prob_y * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})


# increase testing, decrease overall cumuinf

test_that("testing vs prevalence", {

  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0, 9), kappab = rep(0, 9), kappac = rep(0, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("I[0-9][0-9]", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(testing_prob_y =  x$testing_prob_y * 0.5)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})






# increase prep adherence movement, increase overall infections


test_that("adherence movements vs infections", {

  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0, 9), kappab = rep(0, 9), kappac = rep(0, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))


  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(psia =  x$psia * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0, 9), kappab = rep(0, 9), kappac = rep(0, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(psib =  x$psib * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})


# increase prep drop out, increase overall infections


test_that("prep dropout vs infections", {
  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0.1, 9), kappab = rep(0.1, 9), kappac = rep(0.1, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(kappaa =  x$kappaa * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))



  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0.1, 9), kappab = rep(0.1, 9), kappac = rep(0.1, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(kappab =  x$kappab * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))

  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0.1, 9), kappab = rep(0.1, 9), kappac = rep(0.1, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(kappac =  x$kappac * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0.1, 9), kappab = rep(0.1, 9), kappac = rep(0.1, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(kappa1 =  x$kappa1 * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))


})


# increase ART drop out, increase overall infections

test_that("ART dropout vs prevalence", {

  parameters <- lhs_parameters(1,  par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default,
                               forced_pars = list(
                                 time = time_default,
                                 n_y_noncomm = array(data = c(1), dim=c(5, 9, 9)),kappaa = rep(0.1, 9), kappab = rep(0.1, 9), kappac = rep(0.1, 9),
                                 n_y_comm = array(data = c(1), dim=c(5, 9, 9)),ignore_ranges_fc_c = 1, betaMtF_noncomm = 0.001,
                                 eP1a = rep(0.6, 9), eP1b = rep(0.4, 9), eP1c = rep(0, 9),
                                 eP1d = rep(0, 9), eP0 = rep(0, 9),

                                 prep_intervention_y = matrix(c(1), ncol=9, nrow=4),
                                 sigma = c(1,1,1,1,1,1,1,1,1),
                                 test_rate_prep = c(4,1,1,1,1,1,1,1,1),
                                 PrEPOnOff = 1
                               ))

  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]

  xx <- result[c(grep("cumuInf", names(result)))]
  N1 <- rowSums(do.call(cbind, xx))

  parameters <- lapply(parameters, function(x) modifyList(as.list(x), list(phi2 =  x$phi2 * 2)))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  xx <- result[c(grep("cumuInf", names(result)))]
  N2 <- rowSums(do.call(cbind, xx))

  expect_true(sum(N2) > sum(N1))
})

test_that("movement in = out", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  for (i in 1:9)
    expect_equal(sum(result$rate_move_in[1,,i]), -result$rate_move_out[1,i])
})

test_that("if rate_leave_client is 0, then client out should equal zero", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default[which(rownames(ranges_default) != "rate_leave_client"),],
                               forced_pars = list(time = time_default, rate_leave_client = 0))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  expect_true(all(result$rate_move_out[,5] == 0))

})

test_that("if rate_leave_pro_FSW is 0, then FSW out should equal zero", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default[which(rownames(ranges_default) != "rate_leave_pro_FSW"),],
                               forced_pars = list(time = time_default, rate_leave_pro_FSW = 0))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  expect_true(all(result$rate_move_out[,1] == 0))

})

test_that("if rate_leave_low_FSW is 0, then FSW low out should equal zero", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default[which(rownames(ranges_default) != "rate_leave_low_FSW"),],
                               forced_pars = list(time = time_default, rate_leave_low_FSW = 0))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]


  expect_true(all(result$rate_move_out[,2] == 0))

})



test_that("There should be no one but pro FSW on PrEP", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default,
                               forced_pars = list(time = time_default
                                                  ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]




  expect_true(sum(c(result$S1a[,2:9], result$S1b[,2:9], result$S1c[,2:9])) == 0)

})

test_that("fPs sum to 1", {
  parameters <- lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default,
                               ranges = ranges_default,
                               forced_pars = list(time = time_default
                               ))
  result = run_model_for_tests(number_simulations = 1, time = time_default, parameters = parameters)[[1]]




  expect_true(all(colSums(do.call(rbind, result[c("fPa", "fPb", "fPc")])) == 1))

})
