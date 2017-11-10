# helper files have functions that will be used a lot in the model tests


par_gridplot2 = function(result, parm) {
  require(plyr)
  fc_df = aperm(result[parm][[1]], c(2, 3, 1))
  fc_df_list = alply(fc_df, 3)

  fc_df_list_applied = lapply(fc_df_list, function(x) {colnames(x) = rownames(x) = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients",
                                                                             "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou")
  return(x)})

  dat = data.frame(row =
                     rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients",
                           "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou"), 1, each = 9),
                   col =
                     rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients",
                           "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou"), 9),
  value = unlist(lapply(fc_df_list_applied, c)),
  year = unlist(sort(rep(time, 81))))
  dat$row = factor(dat$row, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients",
                                       "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou"))
  dat$col = factor(dat$col, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients",
                                       "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou"))


  return(ggplot(dat, aes(x = year, y = value, color = value)) + geom_line(size = 2) + facet_grid(row~col) + theme_bw())
}



# best_set_default --------------------------------------------------------
best_set_default = list(
  init_clientN_from_PCR=0,
  initial_Ntot = 286114,

  frac_women_ProFSW = 0.0024,
  frac_women_LowFSW = 0.0027,
  frac_women_exFSW = 0.0024,

  frac_men_client = 0.2,
  frac_women_virgin = 0.1,
  frac_men_virgin = 0.1,

  prev_init_FSW = 0.0326,
  prev_init_rest = 0.0012,
  # N_init = c(672, 757, 130895, 672, 27124, 100305, 14544, 11145, 0),
  # fraction_F = 0.5,
  fraction_F = 0.515666224,

  epsilon_1985 = 0.059346131 * 1.5,
  epsilon_1992 = 0.053594832 * 1.5,
  epsilon_2002 = 0.026936907 * 1.5,
  epsilon_2013 = 0.026936907 * 1.5,
  epsilon_2016 = 0.026936907 * 1.5,
  # mu = c(0.02597403, 0.02597403, 0.02597403, 0.02597403, 0.02739726, 0.02739726, 0.02597403, 0.02739726, 0.02597403), # women 1/((27 + 50)/2) # men 1/((25 +  48)/2)
  #   c_comm = c(750, 52, 0, 0, 13.5, 0, 0, 0, 0),
  #   c_noncomm = c(0.38, 0.38, 0.88, 0.88, 4, 1.065, 0, 0, 0), # partner change rate lowlevel FSW same as pro, others are approximations from various surveys
  #
  muF = 0.02597403,
  muM = 0.02739726,
  # PARTNER CHANGE RATE
  c_comm_1985 = c(1229.5, 52, 0, 0, 10.15873, 0, 0, 0, 0), # (1020 + 1439)/2
  c_comm_1993 = c(1229.5, 52, 0, 0, 10.15873, 0, 0, 0, 0), # (1020 + 1439)/2
  c_comm_1995 = c(1280, 52, 0, 0, 10.15873, 0, 0, 0, 0), # (1135 + 1425)/2
  c_comm_1998 = c(881, 52, 0, 0, 10.15873, 0, 0, 0, 0), # (757 + 1005)/2
  c_comm_2002 = c(598.5, 52, 0, 0, 11.08109, 0, 0, 0, 0), # (498 + 699)/2, (13.387-10.15873)/14 * 4 + 10.15873
  c_comm_2005 = c(424, 52, 0, 0, 11.77286, 0, 0, 0, 0), # (366 + 482)/2, (13.387-10.15873)/14 * 7 + 10.15873
  c_comm_2008 = c(371.5, 52, 0, 0, 12.46464, 0, 0, 0, 0), # (272 + 471)/2, (13.387-10.15873)/14 * 10 + 10.15873
  c_comm_2012 = c(541, 52, 0, 0, 13.387, 0, 0, 0, 0), # (459 + 623)/2
  c_comm_2015 = c(400, 52, 0, 0, 17.15294, 0, 0, 0, 0), # (309 + 491)/2
  c_comm_2016 = c(400, 52, 0, 0, 17.15294, 0, 0, 0, 0), # (309 + 491)/2

  c_noncomm_1985 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0), # (0.4682779 + 0.3886719 + 0.2729358)/3
  c_noncomm_1993 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0),
  c_noncomm_1995 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0),
  c_noncomm_1998 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0),
  c_noncomm_2002 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0),
  c_noncomm_2005 = c(0.3766285, 0.3766285, 0.9610526, 0.9610526, 2.028986, 1.337444, 0, 0, 0),
  c_noncomm_2008 = c(0.3766285, 0.3766285, 0.7943578, 0.7943578, 2.028986, 0.7878543, 0, 0, 0),
  c_noncomm_2012 = c(0.3766285, 0.3766285, 0.7943578, 0.7943578, 8.086957, 0.7878543, 0, 0, 0),
  c_noncomm_2015 = c(0.3766285, 0.3766285, 0.7943578, 0.7943578, 6.258258, 0.7878543, 0, 0, 0),
  c_noncomm_2016 = c(0.3766285, 0.3766285, 0.7943578, 0.7943578, 6.258258, 0.7878543, 0, 0, 0),


  #think about transforming to matrix
  betaMtoF_comm = 0.00051, # RR circumcision = 0.44
  betaFtoM_comm = 0.02442*0.44,
  betaMtoF_noncomm = 0.003,
  betaFtoM_noncomm = 0.0038*0.44,

  infect_acute = 9, # RR for acute phase
  infect_AIDS = 2, #7.27, # RR for AIDS phase
  infect_ART =  c(0, rep_len(0, 8)),
  ec = rep_len(0.8, 9), # from kate's paper on nigeria SD couples
  eP0 = c(0, rep_len(0, 8)), # assumptions!
  eP1a = c(0.9, rep_len(0, 8)),
  eP1b = c(0.45, rep_len(0, 8)),
  eP1c = c(0, rep_len(0, 8)),
  eP1d = c(0, rep_len(0, 8)),
  gamma01 = 0.4166667, #years
  SC_to_200_349 = 3.4,
  gamma04 = 4.45, #years

  kappaa = rep(0.2, 9),
  kappab = rep(0.2, 9),
  kappac = rep(0.2, 9),
  kappa1 = rep(0.2, 9),


  alpha01 = rep_len(0, 9),
  alpha02 = rep_len(0, 9),
  alpha03 = rep_len(0.05, 9),
  alpha04 = rep_len(0.08, 9),
  alpha05 = rep_len(0.27, 9), #1/2.9
  alpha11 = rep_len(0, 9),
  alpha22 = rep_len(0, 9),
  alpha23 = rep_len(0.05, 9),
  alpha24 = rep_len(0.08, 9),
  alpha25 = rep_len(0.27, 9),
  alpha32 = rep_len(0, 9),
  alpha33 = rep_len(0.05, 9),
  alpha34 = rep_len(0.08, 9),
  alpha35 = rep_len(0.27, 9),
  alpha42 = rep_len(0, 9),
  alpha43 = rep_len(0.05, 9),
  alpha44 = rep_len(0.08, 9),
  alpha45 = rep_len(0.27, 9),


  #PREP
  zetaa_t = c(1985, 2013, 2015, 2016),
  zetaa_y = matrix(c(rep(0, 9), 0, rep(0, 9-1), rep(0, 9), rep(0, 9)), ncol = 9, byrow = T),
  zetab_t = c(1985, 2013, 2015, 2016),
  zetab_y = matrix(c(rep(0, 9), 0, rep(0, 9-1), rep(0, 9), rep(0, 9)), ncol = 9, byrow = T),
  zetac_t = c(1985, 2013, 2015, 2016),
  zetac_y = matrix(c(rep(0, 9), 0, rep(0, 9-1), rep(0, 9), rep(0, 9)), ncol = 9, byrow = T),
  # zetac_y = matrix(c(rep(0, 9), 0.0075, rep(0, 9-1), rep(0, 9), rep(0, 9)), ncol = 9, byrow = T),

  psia = rep_len(0.1,9),
  psib = rep_len(0.1,9),

  #TESTING

  test_rate_prep = c(4, 0, 0, 0, 0, 0, 0, 0, 0),
  sigma = c(1, 0, 0, 0, 0, 0, 0, 0, 0),
  prep_intervention_t = c(1985, 2013, 2015, 2016),
  prep_intervention_y = matrix(c(rep(0, 9), 1, rep(0, 9-1), rep(0, 9), rep(0, 9)), ncol = 9, byrow = T),

  testing_prob_t = c(1985, 2001, 2005, 2006, 2008, 2012, 2013, 2015, 2016),
  # testing_prob_y = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, # 1985 columns are the risk groups
  #                           0, 0, 0, 0, 0, 0, 0, 0, 0, # 2001
  #                           0, 0, 0, 0, 0, 0, 0, 0, 0, # 2005
  #                           0.142, 0.142, 0.142, 0.142, 0.142, 0.142, 0, 0, 0, # 2006 0.653/8 slope
  #                           0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0, 0, 0, # 2008 3*0.653/8
  #                           0.331, 0.331, 0.331, 0.331, 0.331, 0.331, 0, 0, 0, # 2012 7*0.653/8
  #                           0.331, 0.331, 0.331, 0.331, 0.331, 0.331, 0, 0, 0, # 2013
  #                           0.331, 0.331, 0.331, 0.331, 0.331, 0.331, 0, 0, 0, # 2015
  #                           0.331, 0.331, 0.331, 0.331, 0.331, 0.331, 0, 0, 0), # 2016
  # nrow = 9, ncol = 9, byrow = T),
  testing_prob_y = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, # 1985 columns are the risk groups
                            0, 0, 0, 0, 0, 0, 0, 0, 0, # 2001
                            0, 0, 0, 0, 0, 0, 0, 0, 0, # 2005
                            0.081625, 0.142, 0.142, 0.142, 0.0975, 0.0975, 0, 0, 0, # 2006 0.653/8 slope
                            0.244875, 0.21, 0.21, 0.21, 0.1, 0.1, 0, 0, 0, # 2008 3*0.653/8
                            0.571375, 0.331, 0.331, 0.331, 0.0582, 0.0582, 0, 0, 0, # 2012 7*0.653/8
                            0.653, 0.331, 0.331, 0.331, 0.0582, 0.0582, 0, 0, 0, # 2013
                            0.68, 0.331, 0.331, 0.331, 0.0582, 0.0582, 0, 0, 0, # 2015
                            0.68, 0.331, 0.331, 0.331, 0.0582, 0.0582, 0, 0, 0), # 2016
                          nrow = 9, ncol = 9, byrow = T),
  #ART
  ART_prob_t = c(1985, 2002, 2005, 2016),
  # ART_prob_y = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, # 1985
  #                       0, 0, 0, 0, 0, 0, 0, 0, 0, # 2002
  #                       0.1448571, 0.1448571, 0.1448571, 0.1448571, 0.1448571, 0.1448571, 0, 0, 0, # 2005 0.676/14 * 3
  #                       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0),
  #                     nrow = 4, ncol = 9, byrow = T), # 2016 GP: (0.8+0.552)/2
  ART_prob_y = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, # 1985
                        0, 0, 0, 0, 0, 0, 0, 0, 0, # 2002
                        0, 0.1448571, 0.1448571, 0.1448571, 0.1448571, 0.1448571, 0, 0, 0, # 2005 0.676/14 * 3
                        0.6739, 0.676, 0.676, 0.676, 0.676, 0.676, 0, 0, 0),
                      nrow = 4, ncol = 9, byrow = T), # 2016 GP: (0.8+0.552)/2
  RR_ART_CD4200 = 5.39,
  phi2 = c(0.105360516, rep_len(0.025,8)), # former sex workers drop out rate??!
  phi3 = c(0.105360516, rep_len(0.025,8)),
  phi4 = c(0.105360516, rep_len(0.025,8)),
  phi5 = c(0.105360516, rep_len(0.025,8)),
  ART_RR_prog = (1.3+3.45)/2,
  ART_RR_mort = (1.3+3.45)/2,

  #CONDOM

  fc_y_comm_1985 = matrix(
    c(0, 0, 0, 0, 0.145524, 0, 0, 0, 0, # 0.145524 is using John's FSW condom 1989 as prop of 1993, * our measure of 1993
      0, 0, 0, 0, 0.145524, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.145524, 0.145524, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_1993 = matrix(
    c(0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.536, 0.536, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_1995 = matrix(
    c(0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.536, 0.536, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_1998 = matrix(
    c(0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0.536, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.536, 0.536, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2002 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2005 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2008 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2012 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2015 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_comm_2015 = matrix(
    c(0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0.8, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.8, 0.8, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_1985 = matrix(
    c(0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_1993 = matrix(
    c(0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  # 1998
  # (0.33 + 0.2705314)/ 2 # average FSW client
  # (0.0326087 + 0.2705314)/ 2 # average client GPF
  # (0.0326087 + 0.04989035) / 2 # average gpm gpf

  fc_y_noncomm_1998 = matrix(
    c(0, 0, 0, 0, 0.3002657, 0, 0, 0, 0,
      0, 0, 0, 0, 0.3002657, 0, 0, 0, 0,
      0, 0, 0, 0, 0.15157, 0.04124952, 0, 0, 0,
      0, 0, 0, 0, 0.15157, 0.04124952, 0, 0, 0,
      0.3002657, 0.3002657, 0.15157, 0.15157, 0, 0, 0, 0, 0,
      0, 0, 0.04124952, 0.04124952, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  # 2008
  # (0.33 + 0.4)/ 2 # average FSW client (both approx)
  # ((0.05042017+0.241404781)/2 + 0.4)/ 2 # average client GPF (gpf averaged from 2 estimtes)
  # ((0.05042017+0.241404781)/2 + (0.07103825+0.34838295)/2) / 2 # average gpm gpf

  fc_y_noncomm_2002 = matrix(
    c(0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0.365, 0.365, 0.2729562, 0.2729562, 0, 0, 0, 0, 0,
      0, 0, 0.1778115, 0.1778115, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_2008 = matrix(
    c(0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0.365, 0.365, 0.2729562, 0.2729562, 0, 0, 0, 0, 0,
      0, 0, 0.1778115, 0.1778115, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_2011 = matrix(
    c(0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0.365, 0.365, 0.2729562, 0.2729562, 0, 0, 0, 0, 0,
      0, 0, 0.1778115, 0.1778115, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_2015 = matrix(
    c(0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0.365, 0.365, 0.2729562, 0.2729562, 0, 0, 0, 0, 0,
      0, 0, 0.1778115, 0.1778115, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

  fc_y_noncomm_2016 = matrix(
    c(0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.365, 0, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0, 0, 0, 0, 0.2729562, 0.1778115, 0, 0, 0,
      0.365, 0.365, 0.2729562, 0.2729562, 0, 0, 0, 0, 0,
      0, 0, 0.1778115, 0.1778115, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),



  fc_t_comm = c(1985, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015, 2016),

  fc_t_noncomm = c(1985, 1993, 1998, 2002, 2008, 2011, 2015, 2016),


  n_y_comm_1985 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_y_comm_2002 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_y_comm_2015 = matrix(
    c(0.02), ncol=9, nrow = 9),

  n_y_comm_2016 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_t_comm = c(1985, 2002, 2015, 2016),


  n_y_noncomm_1985 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_y_noncomm_2002 = matrix(
    c(0.02), ncol=9, nrow = 9),

  n_y_noncomm_2015 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_y_noncomm_2016 = matrix(
    c(0.01), ncol=9, nrow = 9),

  n_t_noncomm = c(1985, 2002, 2015, 2016),

  rate_leave_pro_FSW = 0.2,
  FSW_leave_Cotonou_fraction = 0.1,
  rate_leave_low_FSW = 0.1,
  rate_leave_client = 0.05,
  dropout_rate_not_FSW = 0.025,
  replaceDeaths = 0,
  movement = 1,

  ART_recruit_rate_rest = 0.2,
  ART_reinit_rate_FSW = 0.25

)




par_seq_default = c("c_comm", "c_noncomm")
condom_seq_default = c("fc_y_comm", "fc_y_noncomm")
groups_seq_default = c("ProFSW", "LowFSW", "GPF", "FormerFSW", "Client", "GPM", "VirginF", "VirginM", "FormerFSWoutside")
years_seq_default = seq(1985, 2016)
time_default <- seq(1986, 2020, length.out = 35)

# ranges_default ----------------------------------------------------------
ranges_default = rbind(

  # MISC
  init_clientN_from_PCR = c(0,0),
  who_believe_comm = c(0, 1),

  # DEMOGRAPHIC

  fraction_F = c(0.512, 0.52), # fraction of population born female
  frac_women_ProFSW = c(0.0024, 0.0143), # fraction of women that are professional FSW
  frac_men_client = c(0.151, 0.4), # fraction of men that are clients
  frac_women_virgin = c(0.079, 0.2), # fraction of women that are virgins
  frac_men_virgin = c(0.070, 0.17), # fraction of men that are virgins

  prev_init_FSW = c(0.0132, 0.0659), # initial prevalence of FSW
  prev_init_rest = c(0.000313, 0.00294), # initial prevalence of the other groups


  # growth rates
  epsilon_1985 = c(0.08, 0.08),
  epsilon_1992 = c(0.08, 0.08),
  epsilon_2002 = c(0.06, 0.07),
  epsilon_2013 = c(0.04, 0.06),
  epsilon_2016 = c(0.04, 0.06),

  muF = c(0.01851852, 0.025), # female mortality
  muM = c(0.01851852, 0.025), # male mortality

  rate_leave_pro_FSW = c(0, 1), # rate of exit of professional sex work
  rate_leave_low_FSW = c(0, 1), # rate of exit of low level sex work

  fraction_FSW_foreign = c(0.5, 0.9),

  rate_leave_client = c(0, 0.189), # rate of exit of clients

  rate_enter_sexual_pop_F = c(1/(20-15), 1/(17-15)), # rate of entering sexual population women
  rate_enter_sexual_pop_M = c(1/(20-15), 1/(17-15)), # rate of entering sexual population men

  fraction_sexually_active_15_F = c(0.119, 0.17), # fraction of 15 year old women sexually active
  fraction_sexually_active_15_M = c(0.18, 0.35), # fraction of 15 year old men sexually active


  # BEHAVIOURAL

  # commercial partnerships
  c_comm_1993_ProFSW = c(192, 1277),
  c_comm_2005_ProFSW = c(81, 562),
  # c_comm_2015_ProFSW = c(71, 501),

  c_comm_1998_Client = c(8.39, 11.9),
  c_comm_2012_Client = c(11.8, 15),
  c_comm_2015_Client = c(14.5, 19.8),



  #non commercial partnerships
  c_non_comm_1985_ProFSW = c(0.31, 0.86),
  c_non_comm_1985_LowFSW = c(0.41, 1.04),
  c_non_comm_1985_Client= c(1.6, 7.9),

  c_noncomm_1998_GPF = c(0.93, 0.99),
  c_noncomm_2008_GPF = c(0.77, 0.82),

  c_noncomm_1998_GPM = c(1.24, 1.43),
  c_noncomm_2008_GPM = c(0.73, 0.84),


  # sex acts per partnership comm
  n_y_comm_1985_ProFSW_Client = c(1, 10.23),
  n_y_comm_1985_Client_ProFSW = c(1.45, 11.45),

  n_y_comm_1985_LowFSW_Client = c(1, 1),
  n_y_comm_1985_Client_LowFSW = c(1, 1),

  # sex acts per partnership noncomm

  n_y_noncomm_2002_ProFSW_Client = c(13, 20),
  n_y_noncomm_2015_ProFSW_Client = c(38.2, 60),
  n_y_noncomm_1985_GPF_GPM = c(29, 43.7),
  n_y_noncomm_1985_GPM_GPF = c(19.4, 46.7),


  #BETA
  betaMtoF_baseline = c(0.0006, 0.00109), # baseline male to female transmission rate
  RR_beta_FtM = c(0.53, 2), # RR for transmission female to male
  RR_beta_HSV2_comm = c(1.4, 2.1), # RR for commercial sex acts where the susceptible individual is infected HSV2
  RR_beta_HSV2_noncomm = c(2.2, 3.4), # RR for non commercial sex acts where the susceptible individual is infected HSV2
  prev_HSV2_FSW = c(0.8687271, 0.9403027), # prevalence HSV2 in FSW
  prev_HSV2_Client = c(0.14, 0.8687271), # prevalence HSV2 in clients
  prev_HSV2_GPF = c(0.2666742, 0.3236852), # prevalence of HSV2 in GPF
  prev_HSV2_GPM = c(0.09843545, 0.14108970), # prevalence of HSV2 in GPM
  RR_beta_circum = c(0.34, 0.72), # RR for transmission if susceptible individual is circumcised


  # Progression parameters

  infect_acute = c(4.47, 18.81), # RR for transmission rate if infected is acute stage
  infect_AIDS = c(4.45, 11.88), # RR for transmission rate if infected is in AIDS stage

  eff_ART = c(0.96, 1), # infectiousness RR when on ART (efficacy ART assuimed 90% * % undetectable which is 52.3%)

  ec = c(0.58, 0.95), # condom efficacy
  eP1a = c(0.9, 0.9), # prep efficacy perfect adherence
  eP1b = c(0, 0.9), # prep efficacy intermediate adherence
  eP1c = c(0, 0), # prep efficacy poor adherence

  SC_to_200_349 = c(2.2, 4.6),
  gamma04 = c(3.9, 5),

  alpha03 = c(0.03, 0.07),
  alpha04 = c(0.05, 0.12),
  alpha05 = c(0.23, 0.33),

  ART_RR_prog = c(1.3, 3.45),
  ART_RR_mort = c(1.3, 3.45),

  dropout_rate_not_FSW = c(0.0233, 0.274),


  # condoms

  fc_y_comm_1985_ProFSW_Client = c(0, 0),
  fc_y_comm_1993_ProFSW_Client = c(0.535, 0.687),
  fc_y_comm_2002_ProFSW_Client = c(0.536, 0.992),

  fc_y_noncomm_1985_ProFSW_Client = c(0, 0),
  fc_y_noncomm_2002_ProFSW_Client = c(0.19, 0.62),

  fc_y_noncomm_1985_GPF_GPM = 0,
  fc_y_noncomm_1998_GPF_GPM = c(0.0326087, 0.05042017),
  fc_y_noncomm_2011_GPF_GPM = c(0.161, 0.255),

  viral_supp_y_2014_ProFSW = c(0.91, 0.92),
  viral_supp_y_1986_rest = c(0.1, 0.2),
  ART_eff = c(0.96, 1),

  ART_recruit_rate_FSW = c(0.5, 1.5),
  ART_recruit_rate_rest = c(0.5, 1.5),

  ART_reinit_rate_FSW = c(0.25, 1.5),
  ART_reinit_rate_rest = c(0.25, 1.5)






)


outputs_default = c("prev", "frac_N", "Ntot", "epsilon", "rate_leave_client", "alphaItot", "prev_FSW", "prev_LowFSW", "prev_client", "prev_men", "prev_women", "c_comm_balanced", "c_noncomm_balanced", "who_believe_comm")



parameter_names = names(lhs_parameters(1, par_seq = par_seq_default, condom_seq = condom_seq_default, groups_seq = groups_seq_default, years_seq = years_seq_default, set_pars = best_set_default, ranges = ranges_default, time = time_default)[[1]])




all_lambda_pars = c("c_comm", "p_comm", "N", "beta_comm", "R", "fc_comm",
                    "fP_comm", "n_comm", "eP0", "eP1a", "eP1b", "eP1c", "eP1d", "ec", "fc_noncomm", "fP_noncomm",
                    "n_noncomm", "c_noncomm", "p_noncomm", "infect_ART",
                    "infect_acute", "infect_AIDS",  "beta_noncomm")

