
# CREATE parTab


# ranges ------------------------------------------------------------------

# yup
ranges = rbind(



  testing_prob_men_2006 = c(0.0975, 0.21),
  testing_prob_men_2008 = c(0.1, 0.26),
  testing_prob_men_2012 = c(0.058, 0.26), # NOTE 0.26 is from MICHEL 2008

  testing_prob_women_2006 = c(0.142, 0.4),
  testing_prob_women_2008 = c(0.21, 0.54),
  testing_prob_women_2012 = c(0.331, 0.513),


  infected_FSW_incoming = c(1,1),
  n_y_noncomm_1998_GPF_GPM = c(34, 44),
  n_y_noncomm_2011_GPF_GPM = c(29, 38),

  prev_non_ben_fsw_1993 = c(0.027, 0.163),
  # prev_non_ben_fsw_1993 = c(0.12, 0.20),
  # prev_non_ben_fsw_2015 = c(0.03, 0.157),


  prev_non_ben_fsw_2015 = c(0, 0.046),

  # prev_non_ben_fsw_2015 = c(0.10, 0.157), # high
  # prev_non_ben_fsw_2015 = c(0.03, 0.046),

  # MISC
  # init_clientN_from_PCR = c(0,0),
  who_believe_comm = c(0, 0),

  # # growth rates
  # epsilon_1985 = c(0.08, 0.08),
  # epsilon_1992 = c(0.08, 0.08),
  # epsilon_2002 = c(0.06, 0.07),
  # epsilon_2013 = c(0.04, 0.06),
  # epsilon_2016 = c(0.04, 0.06),

  epsilon_1985 = c(0.059, 0.059),
  epsilon_1992 = c(0.048, 0.058),
  epsilon_2002 = c(0.027, 0.027),
  epsilon_2013 = c(0.027, 0.027),
  epsilon_2016 = c(0.027, 0.027),

  # DEMOGRAPHIC

  fraction_F = c(0.512, 0.52), # fraction of population born female
  # frac_women_ProFSW = c(0.0024, 0.0036), # fraction of women that are professional FSW
  frac_women_ProFSW = c(0.0024, 0.00715), # fraction of women that are professional FSW
  frac_women_LowFSW = c(1, 2), # relative abundance of low FSW relative to pro FSW

  frac_men_client = c(0.066, 0.3), # fraction of men that are clients
  frac_women_virgin = c(0.079, 0.2), # fraction of women that are virgins
  frac_men_virgin = c(0.070, 0.17), # fraction of men that are virgins

  # prev_init_FSW = c(0.0132, 0.1), # initial prevalence of FSW
  prev_init_FSW = c(0.0132, 0.0659), # initial prevalence of FSW
  prev_init_rest = c(0.000313, 0.00294), # initial prevalence of the other groups




  muF = c(0.0187, 0.02), # female mortality
  muM = c(0.0194, 0.022), # male mortality


  # rate_leave_pro_FSW = c(0, 0.125),
  rate_leave_pro_FSW = c(0, 0.33), # rate of exit of professional sex work
  # rate_leave_low_FSW = c(0, 1), # rate of exit of low level sex work

  fraction_FSW_foreign = c(0.5, 0.9),
  # fraction_FSW_foreign = c(0.1, 0.5),

  rate_leave_client = c(0, 0.295), # rate of exit of clients
  # rate_leave_client = c(0, 0.2), # rate of exit of clients

  rate_enter_sexual_pop_F = c(0.2, 0.5), # rate of entering sexual population women
  rate_enter_sexual_pop_M = c(0.2, 0.5), # rate of entering sexual population men

  fraction_sexually_active_15_F = c(0.119, 0.17), # fraction of 15 year old women sexually active
  fraction_sexually_active_15_M = c(0.18, 0.35), # fraction of 15 year old men sexually active


  # BEHAVIOURAL

  # commercial partnerships
  c_comm_1993_ProFSW = c(192, 1277),

  c_comm_1995_ProFSW = c(192, 1277),

  c_comm_2005_ProFSW = c(81, 562),
  # c_comm_2015_ProFSW = c(71, 501),

  c_comm_1993_LowFSW = c(26, 78),




  c_comm_1998_Client = c(8.4, 32),

  c_comm_2002_Client = c(11.1, 19.8),
  # c_comm_2012_Client = c(11.8, 15),
  # c_comm_2015_Client = c(14.5, 19.8),



  #non commercial partnerships
  c_noncomm_1985_ProFSW = c(0.31, 0.86),
  c_noncomm_1985_LowFSW = c(0.41, 1.04),
  c_noncomm_1985_Client = c(1.6, 3.3),



  c_noncomm_1998_GPF = c(0.93, 0.99),
  c_noncomm_2008_GPF = c(0.77, 0.82),

  c_noncomm_1998_GPM = c(1.25, 1.43),
  c_noncomm_2008_GPM = c(0.73, 0.84),


  # sex acts per partnership comm
  n_y_comm_1985_ProFSW_Client = c(1, 3.3),
  # n_y_comm_1985_Client_ProFSW = c(1.45, 11.45),
  # n_y_comm_1985_ProFSW_Client = c(1, 4),
  # n_y_comm_2002_ProFSW_Client = c(1, 3),

  n_y_comm_1985_LowFSW_Client = c(1, 1),
  n_y_comm_1985_Client_LowFSW = c(1, 1),

  # sex acts per partnership noncomm

  n_y_noncomm_2002_ProFSW_Client = c(13, 20),
  n_y_noncomm_2015_ProFSW_Client = c(38.2, 60),

  # n_y_noncomm_1985_GPF_GPM = c(39, 100),
  # n_y_noncomm_1985_GPM_GPF = c(19.4, 46.7),


  #BETA
  betaMtoF_baseline = c(0.0006, 0.00109), # baseline male to female transmission rate
  RR_beta_FtM = c(0.53, 2), # RR for transmission female to male
  # RR_beta_HSV2_comm_a = c(1.4, 2.1), # RR for commercial sex acts where the susceptible individual is infected HSV2
  # RR_beta_HSV2_noncomm_a = c(2.2, 3.4), # RR for non commercial sex acts where the susceptible individual is infected HSV2

  RR_beta_HSV2_a_FSW = c(0.9, 2.3),
  RR_beta_HSV2_a_client = c(1.5, 2.2),
  RR_beta_HSV2_a_GPF = c(1.8, 3.4),
  RR_beta_HSV2_a_GPM = c(2.2, 4.3),

  prev_HSV2_FSW = c(0.87, 0.94), # prevalence HSV2 in FSW
  prev_HSV2_Client = c(0.18, 0.28), # prevalence HSV2 in clients
  prev_HSV2_GPF = c(0.27, 0.32), # prevalence of HSV2 in GPF
  prev_HSV2_GPM = c(0.098, 0.14), # prevalence of HSV2 in GPM
  RR_beta_circum = c(0.34, 0.72), # RR for transmission if susceptible individual is circumcised


  # Progression parameters

  infect_acute = c(4.5, 18.8), # RR for transmission rate if infected is acute stage
  infect_AIDS = c(4.5, 11.9), # RR for transmission rate if infected is in AIDS stage

  ART_eff = c(0.96, 0.99), # infectiousness RR when on ART (efficacy ART assuimed 90% * % undetectable which is 52.3%)

  viral_supp_y_1986_rest = c(0.424, 0.85),
  viral_supp_y_2015_ProFSW = c(0.75, 0.85),

  ec = c(0.58, 0.95), # condom efficacy

  # eP1a = c(0.9, 0.9), # prep efficacy perfect adherence
  # eP1b = c(0, 0.9), # prep efficacy intermediate adherence
  # eP1c = c(0, 0), # prep efficacy poor adherence


  SC_to_death = c(8.7, 12.3),
  dur_primary_phase = c(0.25, 0.42),
  dur_200_349 = c(2.3, 4.4),
  dur_below_200 = c(0.58, 3.17),


  alpha03 = c(0.01, 0.05),
  alpha04 = c(0.03, 0.1),

  ART_RR_prog = c(4.82, 10.23),

  # intervention_testing_increase = c(1, 2),
  # intervention_testing_increase = c(0.5, 2), # keep
  intervention_testing_increase = c(0, 0),

  RR_test_CD4200 = c(1, 5.4),

  # ART_recruit_rate_FSW = c(0.5, 6),
  # ART_recruit_rate_FSW = c(0.5, 1.5),
  # ART_recruit_rate_FSW = c(0.5, 3),
  ART_recruit_rate_FSW = c(0, 5),

  # ART_recruit_rate_rest = c(0.5, 1.5),
  # ART_recruit_rate_rest = c(0.5, 6),

  ART_recruit_rate_rest = c(3, 12), # NEW
  ART_init_ratio_MF = c(1.76, 3.8), # NEW

  # intervention_ART_increase = c(0, 12),
  # intervention_ART_increase = c(0, 24),
  # intervention_ART_increase = c(0.5, 5), # keep
  intervention_ART_increase = c(0, 0),


  dropout_rate_not_FSW = c(0.0233, 0.11),
  dropout_rate_FSW = c(0.0233, 0.11),

  ART_reinit_rate_FSW = c(0.25, 1.5),
  ART_reinit_rate_rest = c(0.25, 1.5),


  # condoms

  # fc_y_comm_1985_ProFSW_Client = c(0, 0),
  # fc_y_comm_1985_ProFSW_Client = c(0.54, 0.69),
  # fc_y_comm_1998_ProFSW_Client = c(0.54, 0.99),


  fc_y_comm_1985_ProFSW_Client = c(0, 0.18),

  fc_y_comm_1993_ProFSW_Client = c(0.18, 0.33),

  fc_y_comm_1998_ProFSW_Client = c(0.4, 0.73),

  fc_y_comm_2002_ProFSW_Client = c(0.61, 0.99),

  fc_y_comm_2008_ProFSW_Client = c(0.86, 0.99),


  fc_y_comm_1985_LowFSW_Client = 0,
  fc_y_comm_2015_LowFSW_Client = c(0.25, 0.52),

  fc_y_noncomm_1985_ProFSW_Client = 0,

  fc_y_noncomm_2002_ProFSW_Client = c(0.19, 0.62),

  fc_y_noncomm_1985_LowFSW_Client = 0,
  fc_y_noncomm_2015_LowFSW_Client = c(0.138, 0.383),

  # fc_y_noncomm_1985_GPF_GPM = 0,
  fc_y_noncomm_1998_GPF_GPM = c(0.033, 0.05),
  fc_y_noncomm_2011_GPF_GPM = c(0.16, 0.26)





)







# ranges end ------------------------------------------------------------------

# who_believe_comm = c(0, 0), !!!!


parTab_create = data.frame(
  names = rownames(ranges),
  values = ranges[,1],
  fixed = 0,
  steps = 0.1,
  lower_bound = ranges[,1],
  upper_bound = ranges[,2]

)


good_previous_fit = get(load("best_pars_combined_0904_top.Rdata"))

good_previous_fit_variedpars = good_previous_fit[[1]][rownames(ranges)]

j=0
for(i in 1:length(rownames(ranges)))
{
  if(parTab_create[i,"names"] %in% names(good_previous_fit_variedpars))
  {
    parTab_create[i, "values"] = as.numeric(good_previous_fit_variedpars[as.character(parTab_create[i,"names"])][[1]][1]);
    j=j+1
  }


}
sum(names(good_previous_fit_variedpars) == rownames(ranges), na.rm = T)





parTab_create[parTab_create$lower_bound ==parTab_create$upper_bound, "fixed"] = 1


rownames(parTab_create) = NULL
write.csv(parTab_create, file = "parTab.csv")



###

create_lik <- function(parTab, data, PRIOR_FUNC,...){


  epi_start = 1986
  # epi_end = 2030
  epi_end = 2017

  # setup -------------------------------------------------------------------
  par_seq = c("c_comm", "c_noncomm")
  condom_seq = c("fc_y_comm", "fc_y_noncomm", "n_y_comm", "n_y_noncomm")
  groups_seq = c("ProFSW", "LowFSW", "GPF", "FormerFSW", "Client", "GPM", "VirginF", "VirginM", "FormerFSWoutside")
  years_seq = seq(1985, 2016)
  time <- seq(epi_start, epi_end, length.out = epi_end - epi_start + 1)
  time_with_mid <- seq(epi_start, epi_end, length.out = (epi_end - epi_start + 0.5)*2)

  #####################################################

  # this is the best set of parameters (the fixed ones)
  # best_set ----------------------------------------------------------------


  best_set = list(
    init_clientN_from_PCR = 0,
    initial_Ntot = 286114,

    frac_women_ProFSW = 0.0024,
    frac_women_LowFSW = 0.0027,
    frac_women_exFSW = 0.0024,

    frac_men_client = 0.2,
    frac_women_virgin = 0.1,
    frac_men_virgin = 0.1,

    prev_init_FSW = 0.0326,
    prev_init_rest = 0.0012,

    nu = 0.022,

    # N_init = c(672, 757, 130895, 672, 27124, 100305, 14544, 11145, 0),
    # fraction_F = 0.5,
    fraction_F = 0.515666224,

    epsilon_1985 = 0.08,
    epsilon_1992 = 0.08,
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
    infect_ART = 0.9 * 0.523, # infectiousness RR when on ART (efficacy ART assuimed 90% * % undetectable which is 52.3%)
    ec = rep_len(0.8, 9), # from kate's paper on nigeria SD couples
    eP0 = c(0, rep_len(0, 8)), # assumptions!
    eP1a = c(0.9, rep_len(0, 8)),
    eP1b = c(0.45, rep_len(0, 8)),
    eP1c = c(0, rep_len(0, 8)),
    eP1d = c(0, rep_len(0, 8)),





    # gamma01 = 0.4166667, #years
    # gamma04 = 4.45, #years
    #

    alpha01 = rep_len(0, 9),
    alpha11 = rep_len(0, 9),
    alpha02 = rep_len(0, 9),
    alpha03 = 0.03,
    alpha04 = 0.07,
    alpha05 = 2,
    alpha11 = rep_len(0, 9),
    alpha22 = rep_len(0, 9),
    # alpha23 = rep_len(0.05, 9),
    # alpha24 = rep_len(0.08, 9),
    # alpha25 = rep_len(0.27, 9),
    alpha32 = rep_len(0, 9),
    # alpha33 = rep_len(0.05, 9),
    # alpha34 = rep_len(0.08, 9),
    # alpha35 = rep_len(0.27, 9),
    alpha42 = rep_len(0, 9),
    # alpha43 = rep_len(0.05, 9),
    # alpha44 = rep_len(0.08, 9),
    # alpha45 = rep_len(0.27, 9),


    test_rate_prep = c(4, 0, 0, 0, 0, 0, 0, 0, 0),
    sigma = c(0.82, 0, 0, 0, 0, 0, 0, 0, 0),
    prep_intervention_t = c(1985, 2015, 2016, 2016.001),
    prep_intervention_y = matrix(c(rep(0, 9), 1, rep(0, 9-1), 1, rep(0, 9-1), rep(0, 9)), ncol = 9, byrow = T),
    PrEPOnOff = 0,

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
                              0, 0.118, 0.118, 0.118, 0.08125, 0.08125, 0, 0, 0, # 2005 0.142*5/6 0.0975*5/6
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
    # phi2 = c(0.105360516, rep_len(0.025,8)), # former sex workers drop out rate??!
    # phi3 = c(0.105360516, rep_len(0.025,8)),
    # phi4 = c(0.105360516, rep_len(0.025,8)),
    # phi5 = c(0.105360516, rep_len(0.025,8)),

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
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_comm_2002 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_comm_2015 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_comm_2016 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_t_comm = c(1985, 2002, 2015, 2016),


    n_y_noncomm_1985 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_noncomm_2002 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_noncomm_1998 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_noncomm_2011 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_noncomm_2015 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_y_noncomm_2016 = matrix(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 9),

    n_t_noncomm = c(1985, 1998, 2002, 2011, 2015, 2016),

    rate_leave_pro_FSW = 0.2,
    FSW_leave_Cotonou_fraction = 0.1,
    rate_leave_low_FSW = 0.1,
    rate_leave_client = 0.05,
    dropout_rate_not_FSW = 0.025,
    replaceDeaths = 0,
    movement = 1,

    ART_recruit_rate_rest = 0.25,
    ART_recruit_rate_FSW = 0.25,

    ART_reinit_rate_FSW = 0.2,
    ART_reinit_rate_rest = 0.2

  )



  #







  # ranges ------------------------------------------------------------------

  # yup
  ranges = rbind(



    testing_prob_men_2006 = c(0.0975, 0.21),
    testing_prob_men_2008 = c(0.1, 0.26),
    testing_prob_men_2012 = c(0.058, 0.26), # NOTE 0.26 is from MICHEL 2008

    testing_prob_women_2006 = c(0.142, 0.4),
    testing_prob_women_2008 = c(0.21, 0.54),
    testing_prob_women_2012 = c(0.331, 0.513),


    infected_FSW_incoming = c(1,1),
    n_y_noncomm_1998_GPF_GPM = c(34, 44),
    n_y_noncomm_2011_GPF_GPM = c(29, 38),

    prev_non_ben_fsw_1993 = c(0.027, 0.163),
    # prev_non_ben_fsw_1993 = c(0.12, 0.20),
    # prev_non_ben_fsw_2015 = c(0.03, 0.157),


    prev_non_ben_fsw_2015 = c(0, 0.046),

    # prev_non_ben_fsw_2015 = c(0.10, 0.157), # high
    # prev_non_ben_fsw_2015 = c(0.03, 0.046),

    # MISC
    # init_clientN_from_PCR = c(0,0),
    who_believe_comm = c(0, 1),

    # # growth rates
    # epsilon_1985 = c(0.08, 0.08),
    # epsilon_1992 = c(0.08, 0.08),
    # epsilon_2002 = c(0.06, 0.07),
    # epsilon_2013 = c(0.04, 0.06),
    # epsilon_2016 = c(0.04, 0.06),

    epsilon_1985 = c(0.059, 0.059),
    epsilon_1992 = c(0.048, 0.058),
    epsilon_2002 = c(0.027, 0.027),
    epsilon_2013 = c(0.027, 0.027),
    epsilon_2016 = c(0.027, 0.027),

    # DEMOGRAPHIC

    fraction_F = c(0.512, 0.52), # fraction of population born female
    # frac_women_ProFSW = c(0.0024, 0.0036), # fraction of women that are professional FSW
    frac_women_ProFSW = c(0.0024, 0.00715), # fraction of women that are professional FSW
    frac_women_LowFSW = c(1, 2), # relative abundance of low FSW relative to pro FSW

    frac_men_client = c(0.066, 0.3), # fraction of men that are clients
    frac_women_virgin = c(0.079, 0.2), # fraction of women that are virgins
    frac_men_virgin = c(0.070, 0.17), # fraction of men that are virgins

    # prev_init_FSW = c(0.0132, 0.1), # initial prevalence of FSW
    prev_init_FSW = c(0.0132, 0.0659), # initial prevalence of FSW
    prev_init_rest = c(0.000313, 0.00294), # initial prevalence of the other groups




    muF = c(0.0187, 0.02), # female mortality
    muM = c(0.0194, 0.022), # male mortality


    # rate_leave_pro_FSW = c(0, 0.125),
    rate_leave_pro_FSW = c(0, 0.33), # rate of exit of professional sex work
    # rate_leave_low_FSW = c(0, 1), # rate of exit of low level sex work

    fraction_FSW_foreign = c(0.5, 0.9),
    # fraction_FSW_foreign = c(0.1, 0.5),

    rate_leave_client = c(0, 0.295), # rate of exit of clients
    # rate_leave_client = c(0, 0.2), # rate of exit of clients

    rate_enter_sexual_pop_F = c(0.2, 0.5), # rate of entering sexual population women
    rate_enter_sexual_pop_M = c(0.2, 0.5), # rate of entering sexual population men

    fraction_sexually_active_15_F = c(0.119, 0.17), # fraction of 15 year old women sexually active
    fraction_sexually_active_15_M = c(0.18, 0.35), # fraction of 15 year old men sexually active


    # BEHAVIOURAL

    # commercial partnerships
    c_comm_1993_ProFSW = c(192, 1277),

    c_comm_1995_ProFSW = c(192, 1277),

    c_comm_2005_ProFSW = c(81, 562),
    # c_comm_2015_ProFSW = c(71, 501),

    c_comm_1993_LowFSW = c(26, 78),




    c_comm_1998_Client = c(8.4, 32),

    c_comm_2002_Client = c(11.1, 19.8),
    # c_comm_2012_Client = c(11.8, 15),
    # c_comm_2015_Client = c(14.5, 19.8),



    #non commercial partnerships
    c_noncomm_1985_ProFSW = c(0.31, 0.86),
    c_noncomm_1985_LowFSW = c(0.41, 1.04),
    c_noncomm_1985_Client = c(1.6, 3.3),



    c_noncomm_1998_GPF = c(0.93, 0.99),
    c_noncomm_2008_GPF = c(0.77, 0.82),

    c_noncomm_1998_GPM = c(1.25, 1.43),
    c_noncomm_2008_GPM = c(0.73, 0.84),


    # sex acts per partnership comm
    n_y_comm_1985_ProFSW_Client = c(1, 3.3),
    # n_y_comm_1985_Client_ProFSW = c(1.45, 11.45),
    # n_y_comm_1985_ProFSW_Client = c(1, 4),
    # n_y_comm_2002_ProFSW_Client = c(1, 3),

    n_y_comm_1985_LowFSW_Client = c(1, 1),
    n_y_comm_1985_Client_LowFSW = c(1, 1),

    # sex acts per partnership noncomm

    n_y_noncomm_2002_ProFSW_Client = c(13, 20),
    n_y_noncomm_2015_ProFSW_Client = c(38.2, 60),

    # n_y_noncomm_1985_GPF_GPM = c(39, 100),
    # n_y_noncomm_1985_GPM_GPF = c(19.4, 46.7),


    #BETA
    betaMtoF_baseline = c(0.0006, 0.00109), # baseline male to female transmission rate
    RR_beta_FtM = c(0.53, 2), # RR for transmission female to male
    # RR_beta_HSV2_comm_a = c(1.4, 2.1), # RR for commercial sex acts where the susceptible individual is infected HSV2
    # RR_beta_HSV2_noncomm_a = c(2.2, 3.4), # RR for non commercial sex acts where the susceptible individual is infected HSV2

    RR_beta_HSV2_a_FSW = c(0.9, 2.3),
    RR_beta_HSV2_a_client = c(1.5, 2.2),
    RR_beta_HSV2_a_GPF = c(1.8, 3.4),
    RR_beta_HSV2_a_GPM = c(2.2, 4.3),

    prev_HSV2_FSW = c(0.87, 0.94), # prevalence HSV2 in FSW
    prev_HSV2_Client = c(0.18, 0.28), # prevalence HSV2 in clients
    prev_HSV2_GPF = c(0.27, 0.32), # prevalence of HSV2 in GPF
    prev_HSV2_GPM = c(0.098, 0.14), # prevalence of HSV2 in GPM
    RR_beta_circum = c(0.34, 0.72), # RR for transmission if susceptible individual is circumcised


    # Progression parameters

    infect_acute = c(4.5, 18.8), # RR for transmission rate if infected is acute stage
    infect_AIDS = c(4.5, 11.9), # RR for transmission rate if infected is in AIDS stage

    ART_eff = c(0.96, 0.99), # infectiousness RR when on ART (efficacy ART assuimed 90% * % undetectable which is 52.3%)

    viral_supp_y_1986_rest = c(0.424, 0.85),
    viral_supp_y_2015_ProFSW = c(0.75, 0.85),

    ec = c(0.58, 0.95), # condom efficacy

    # eP1a = c(0.9, 0.9), # prep efficacy perfect adherence
    # eP1b = c(0, 0.9), # prep efficacy intermediate adherence
    # eP1c = c(0, 0), # prep efficacy poor adherence


    SC_to_death = c(8.7, 12.3),
    dur_primary_phase = c(0.25, 0.42),
    dur_200_349 = c(2.3, 4.4),
    dur_below_200 = c(0.58, 3.17),


    alpha03 = c(0.01, 0.05),
    alpha04 = c(0.03, 0.1),

    ART_RR_prog = c(4.82, 10.23),

    # intervention_testing_increase = c(1, 2),
    # intervention_testing_increase = c(0.5, 2), # keep
    intervention_testing_increase = c(0, 0),

    RR_test_CD4200 = c(1, 5.4),

    # ART_recruit_rate_FSW = c(0.5, 6),
    # ART_recruit_rate_FSW = c(0.5, 1.5),
    # ART_recruit_rate_FSW = c(0.5, 3),
    ART_recruit_rate_FSW = c(0, 5),

    # ART_recruit_rate_rest = c(0.5, 1.5),
    # ART_recruit_rate_rest = c(0.5, 6),

    ART_recruit_rate_rest = c(3, 12), # NEW
    ART_init_ratio_MF = c(1.76, 3.8), # NEW

    # intervention_ART_increase = c(0, 12),
    # intervention_ART_increase = c(0, 24),
    # intervention_ART_increase = c(0.5, 5), # keep
    intervention_ART_increase = c(0, 0),


    dropout_rate_not_FSW = c(0.0233, 0.11),
    dropout_rate_FSW = c(0.0233, 0.11),

    ART_reinit_rate_FSW = c(0.25, 1.5),
    ART_reinit_rate_rest = c(0.25, 1.5),


    # condoms

    # fc_y_comm_1985_ProFSW_Client = c(0, 0),
    # fc_y_comm_1985_ProFSW_Client = c(0.54, 0.69),
    # fc_y_comm_1998_ProFSW_Client = c(0.54, 0.99),


    fc_y_comm_1985_ProFSW_Client = c(0, 0.18),

    fc_y_comm_1993_ProFSW_Client = c(0.18, 0.33),

    fc_y_comm_1998_ProFSW_Client = c(0.4, 0.73),

    fc_y_comm_2002_ProFSW_Client = c(0.61, 0.99),

    fc_y_comm_2008_ProFSW_Client = c(0.86, 0.99),


    fc_y_comm_1985_LowFSW_Client = 0,
    fc_y_comm_2015_LowFSW_Client = c(0.25, 0.52),

    fc_y_noncomm_1985_ProFSW_Client = 0,

    fc_y_noncomm_2002_ProFSW_Client = c(0.19, 0.62),

    fc_y_noncomm_1985_LowFSW_Client = 0,
    fc_y_noncomm_2015_LowFSW_Client = c(0.138, 0.383),

    # fc_y_noncomm_1985_GPF_GPM = 0,
    fc_y_noncomm_1998_GPF_GPM = c(0.033, 0.05),
    fc_y_noncomm_2011_GPF_GPM = c(0.16, 0.26)





  )








  # outputs -----------------------------------------------------------------
  outputs = c("HIV_positive_On_ART","lambda_sum_0","HIV_positive","pc_of_FOI_on_clients_from_pro_FSW", "S0", "S1a", "S1b", "S1c", "S1d", "prev", "frac_N", "Ntot", "epsilon", "rate_leave_client", "alphaItot", "prev_FSW", "prev_LowFSW", "prev_client", "prev_men", "prev_women", "c_comm_balanced", "c_noncomm_balanced", "who_believe_comm", "ART_coverage_FSW", "ART_coverage_men", "ART_coverage_women", "ART_coverage_all", "rho", "n_comm", "n_noncomm", "fc_comm", "fc_noncomm", "N", "cumuHIVDeaths", "lambda_0", "lambda_1a", "lambda_1b", "lambda_1c", "lambda_1d")


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
                                     100*0.084, 100*0.043),
                           lower = c(3.3, 8.2, 19.2, 48.02, 43.02, 36.58, 31.97, 30.42, 24.93, 23.01, 15.71,
                                     100*0.05898524, 100*0.068218538, 100*0.04293149, 100*0.034772131, 100*0.012660836, 100*0.006039259,
                                     100*0.024181624, 100*0.030073668, 100*0.012980254,
                                     100*0.022857312, 100*0.012427931, 100*0.007517563,
                                     # 100*0.091838441, 100*0.026704897),
                                     100*0.055700329, 100*0.024043597),

                           upper = c(3.3, 8.2, 19.2, 58.48, 54.42, 44.67, 46.27, 39.38, 33.88, 32.23, 22.01,
                                     100*0.11561791, 100*0.115608811, 100*0.105215792, 100*0.090216628, 100*0.051602442, 100*0.035338436,
                                     100*0.047726245, 100*0.052817187, 100*0.035296286,
                                     100*0.047183668, 100*0.029774338, 100*0.028546718,
                                     100*0.120857355, 100*0.069311506))
  prev_points_all = prev_points
  prev_points = prev_points[-c(1,2,3),]



  prev_points_low_fsw_2015 = prev_points
  prev_points_low_fsw_2015[prev_points_low_fsw_2015$time == 2015 & prev_points_low_fsw_2015$variable == "Pro FSW", "lower"] = 13.79


  prev_points_extended_low = data.frame(time = c(1986, 1987, 1988, 1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015,
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
                                                  100*0.084, 100*0.043),
                                        lower = c(3.3, 8.2, 19.2, 48.02, 43.02, 36.58, 31.97, 30.42, 24.93, 23.01, 15.71,
                                                  100*0.05898524, 100*0.068218538, 100*0.04293149, 100*0.034772131, 100*0.012660836, 100*0.006039259,
                                                  100*0.024181624, 100*0.030073668, 100*0.012980254,
                                                  100*0.022857312, 100*0.012427931, 100*0.007517563,
                                                  100*0, 100*0),

                                        upper = c(3.3, 8.2, 19.2, 58.48, 54.42, 44.67, 46.27, 39.38, 33.88, 32.23, 22.01,
                                                  100*0.11561791, 100*0.115608811, 100*0.105215792, 100*0.090216628, 100*0.051602442, 100*0.035338436,
                                                  100*0.047726245, 100*0.052817187, 100*0.035296286,
                                                  100*0.047183668, 100*0.029774338, 100*0.028546718,
                                                  100*0.120857355, 100*0.069311506))






  prev_points_FSW_only = data.frame(time = c(1993, 1995, 1998, 2002, 2005, 2008, 2012, 2015
  ),
  variable = c(rep("Pro FSW", 8)
  ),
  value = c(53.3, 48.7, 40.6, 38.9, 34.8, 29.3, 27.4, 18.7
  ),
  lower = c(48.02, 43.02, 36.58, 31.97, 30.42, 24.93, 23.01, 15.71
  ),
  upper = c(58.48, 54.42, 44.67, 46.27, 39.38, 33.88, 32.23, 22.01))


  prev_points_FSW_only_all_8_LB = prev_points_FSW_only
  prev_points_FSW_only_all_8_LB[prev_points_FSW_only_all_8_LB$time == 2015,"lower"] = 13.79

  prev_points_FSW_only_no_2012 <- prev_points_FSW_only_all_8_LB[-7,]



  prev_points_FSW_only_even_less_2 = prev_points_FSW_only[c(1, 4, 6, 8),]
  prev_points_FSW_Cotonou_centrale_lower_bound = prev_points_FSW_only_even_less_2
  prev_points_FSW_Cotonou_centrale_lower_bound[prev_points_FSW_Cotonou_centrale_lower_bound$time == 2015,"lower"] = 13.79



  prev_points_FSW_Cotonou_centrale_lower_bound_mid_year = prev_points_FSW_Cotonou_centrale_lower_bound
  prev_points_FSW_Cotonou_centrale_lower_bound_mid_year$time = prev_points_FSW_Cotonou_centrale_lower_bound_mid_year$time + 0.5
  prev_points_FSW_Cotonou_centrale_lower_bound$x = c(195, 74, 122, 116)
  prev_points_FSW_Cotonou_centrale_lower_bound$N = c(366, 190, 417, 620)

  # frac N data points ------------------------------------------------------
  frac_N_data_points = data.frame(time = c(1998, 2014,
                                           1998, 1998,
                                           1998, 2008, 2011,
                                           1998, 2008, 2011),
                                  point = c(1.43*0.515666224, 0.24*0.515666224,
                                            100*0.195738802*(1-0.515666224), 40*(1-0.515666224),
                                            100*0.1292392*0.515666224, 100*0.0972973*0.515666224, 100*0.18*0.515666224,
                                            100*0.124632*(1-0.515666224), 100*0.08840413*(1-0.515666224), 100*0.1175*(1-0.515666224)),
                                  variable = c("Pro FSW", "Pro FSW",
                                               "Clients", "Clients",
                                               "Virgin female", "Virgin female", "Virgin female",
                                               "Virgin male", "Virgin male", "Virgin male"))

  # frac_N_discard_points = data.frame(variable = c("Pro FSW", "Clients", "Virgin female", "Virgin male", "Low-level FSW"),
  #                                    min = c(0.001237599, 0.1509*(1-0.515666224), 0.07896475*0.515666224, 0.07039551*(1-0.515666224), 2*0.001237599),
  #                                    max = c(0.007374027, 0.40 * (1-0.515666224), 0.2*0.515666224, 0.17*(1-0.515666224), 5*0.007374027))


  frac_N_discard_points = data.frame(variable = c("Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"),
                                     min = c(0.001237599, 0.074*(1-0.515666224), 0.07896475*0.515666224, 0.07039551*(1-0.515666224), 0.0048*0.516, 1),
                                     max = c(0.0143*0.515666224/2, 0.3 * (1-0.515666224), 0.2*0.515666224, 0.17*(1-0.515666224),  0.0143*0.516, 5))

  frac_N_discard_points_no_FSW_LB = data.frame(variable = c("Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"),
                                               min = c(0, 0.074*(1-0.515666224), 0.07896475*0.515666224, 0.07039551*(1-0.515666224), 0.0048*0.516, 1),
                                               max = c(0.0143*0.515666224/2, 0.3 * (1-0.515666224), 0.2*0.515666224, 0.17*(1-0.515666224),  0.0143*0.516, 5))

  frac_N_discard_points_no_FSW_LB

  # Ntot data points ------------------------------------------------------

  Ntot_data_points = data.frame(time = c(1992, 2002, 2013, 2020, 2030),
                                point = c(404359.0418, 681559.032, 913029.606, 1128727.062, 1423887.65),
                                lower = c(343705.15, 579325.15, 776075.5, 959417.95, 1210304.8),
                                upper = c(465012.85, 783792.85, 1049984.5, 1298036.05, 1637471.2),
                                colour = c("data", "data", "data", "predicted", "predicted"))

  # ART coverage data points ------------------------------------------------







  ART_data_points = data.frame(time = c(2011, 2012, 2013, 2014, 2015, 2016, 2017,

                                        2012, 2014, 2015, 2016, 2017
  ),
  Lower = c(0.33, 0.42, 0.44, 0.39, 0.44, 0.51, 0.57,

            0.09,	0.14,	0.2,	0.22,	0.3

  ),
  Upper = c(0.52, 0.66, 0.69, 0.61, 0.69, 0.8, 0.8,
            0.13, 0.2, 0.28, 0.3, 0.42

  ),
  variable = c("All", "All", "All", "All", "All", "All", "All",
               "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))



  ART_data_points_extra = data.frame(time = c(2011, 2012, 2013, 2014, 2015, 2016, 2017,

                                              2012, 2014, 2015, 2016, 2017, 2017.5
  ),
  Lower = c(0.33, 0.42, 0.44, 0.39, 0.44, 0.51, 0.57,

            0.08,	0.14,	0.2,	0.50,	0.56, 0.8

  ),
  Upper = c(0.52, 0.66, 0.69, 0.61, 0.69, 0.8, 0.8,
            0.11, 0.2, 0.28, 0.7, 0.79, 0.9

  ),
  variable = c("All", "All", "All", "All", "All", "All", "All",
               "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))



  ART_data_points_FSW = ART_data_points[c(8, 9, 10, 11, 12),]

  ART_data_points_FSW_last_2 = ART_data_points[c(11, 12),]

  ART_data_points_first_and_last_FSW = ART_data_points[c(8, 12),]


  # first and last GP, first FSW and 2016, 2017 FSW
  ART_data_points_1611 = ART_data_points[c(1, 7, 8, 11, 12),]

  # first and last GP, all fsw before intervention
  ART_data_points_1911 = ART_data_points[c(1, 7, 8, 9, 10),]

  # NUMBERS OF FSW ON ART

  ART_data_points_with_numbers = data.frame(time = c(2011, 2017,

                                                     2012, 2014, 2015#, 2016, 2017
  ),
  Lower = c(0.33, 0.57,

            27,	34,	42#,	45,	62

  ),
  Upper = c(0.52, 0.8,
            37, 46, 56#, 83, 107

  ),
  variable = c("All", "All",
               "Numbers FSW", "Numbers FSW", "Numbers FSW"#, "Numbers FSW", "Numbers FSW"
  ))

  ART_data_points_with_numbers_reduced = ART_data_points_with_numbers[c(1,2,5),]


  ART_data_points_NULL = data.frame(time = c(2011
  ),
  Lower = c(-1

  ),
  Upper = c(1

  ),
  variable = c("All"))

  # PrEP_fitting ------------------------------------------------

  # PrEP_fitting = data.frame(time = c(2016, 2017, 2016, 2017, 2016, 2017),
  #                           group = c("S1a", "S1a", "S1b", "S1b", "S1c", "S1c"),
  #                           lower = c(50, 40, 50, 40, 50, 40),
  #                           point = c(55, 45, 50, 40, 50, 40),
  #                           upper = c(61, 66, 61, 66, 61, 66)
  #
  #
  # )
  #
  # PrEP_fitting = data.frame(time = c(2016, 2017),
  #                           group = c("S1a", "S1a"),
  #                           lower = c(50, 40),
  #                           point = c(55, 45),
  #                           upper = c(61, 66)
  #
  #
  # )

  PrEP_fitting = NULL

  #
  #
  # ART_data_points_FSW = data.frame(time = c(2010, 2013, 2014, 2015, 2016, 2017
  # ),
  # Lower = c(0.07,	0.13,	0.17,	0.50,	0.56, 0.80
  #
  # ),
  # Upper = c(0.10, 0.18, 0.24, 0.7, 0.79, 0.9
  #
  # ),
  # variable = c("Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))


  # ART_data_points_FSW_last_3 = data.frame(time = c(2015, 2016, 2017
  # ),
  # Lower = c(0.50,	0.56, 0.80
  #
  # ),
  # Upper = c(0.7, 0.79, 0.9
  #
  # ),
  # variable = c("Pro FSW", "Pro FSW", "Pro FSW"))

  # ART_data_points_first_and_last = data.frame(time = c(2010, 2017
  # ),
  # Lower = c(0.07, 0.80
  #
  # ),
  # Upper = c(0.10, 0.9
  #
  # ),
  # variable = c("Pro FSW", "Pro FSW"))




  # ART_data_points_all = data.frame(time = c(2010, 2011, 2012, 2013, 2014, 2015, 2016,
  #                                           2010, 2011, 2012, 2013, 2014, 2015, 2016,
  #                                           2010, 2011, 2012, 2013, 2014, 2015, 2016,
  #                                           2010, 2013, 2014, 2015, 2016
  # ),
  # Lower = c(0.32, 0.4, 0.43, 0.38, 0.43, 0.49, 0.552,
  #           0.32, 0.4, 0.43, 0.38, 0.43, 0.49, 0.552,
  #           0.32, 0.4, 0.43, 0.38, 0.43, 0.49, 0.552,
  #           7.5,	12.5,	17.4,	50.0,	56.1
  #
  # ),
  # Upper = c(0.522449, 0.653061, 0.702041, 0.620408, 0.702041, 0.8, 0.8,
  #           0.522449, 0.653061, 0.702041, 0.620408, 0.702041, 0.8, 0.8,
  #           0.522449, 0.653061, 0.702041, 0.620408, 0.702041, 0.8, 0.8,
  #           0.522449, 0.620408, 0.702041, 0.8, 0.8
  #
  # ),
  # variable = c("Women", "Women", "Women", "Women", "Women", "Women", "Women",
  #              "Men", "Men", "Men", "Men", "Men", "Men", "Men",
  #              "All", "All", "All", "All", "All", "All", "All",
  #              "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))

  #####################################################

  ART_data_points_lazymcmc = data.frame(time = c(2017, 2017, 2015),
                                        variable = c("Men", "Women","Pro FSW"),
                                        x = c(3908, 9784, 49))


  number_simulations=1





  par_names <- rownames(ranges)


  # pars = parTab_create$values

  ## using the `...` bit.
  likelihood_func <- function(pars){


    # put pars into way in which the model can read it

    # names(pars) = par_names


    ranges_new = cbind(pars, pars)
    rownames(ranges_new) = par_names

    parameters = cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                         ranges = ranges_new, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)



    #ranges = NULL needs to be changed to incorporate those that vary!

    # run model
    # res = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)


    res = cotonou::return_outputs(parameters[[1]], cotonou::main_model, time = time, outputs = outputs)

    # calculate likelihood
    # likelihood_list = lapply(res, cotonou::likelihood_lazymcmc, time = time,
    #                          prev_points = prev_points_FSW_Cotonou_centrale_lower_bound,
    #                          frac_N_discard_points = frac_N_discard_points_no_FSW_LB,
    #                          Ntot_data_points = Ntot_data_points,
    #                          ART_data_points = ART_data_points_lazymcmc,
    #                          PrEP_fitting = PrEP_fitting)

    lik = cotonou::likelihood_lazymcmc(res, time = time,
                        prev_points = prev_points_FSW_Cotonou_centrale_lower_bound,
                        frac_N_discard_points = frac_N_discard_points_no_FSW_LB,
                        Ntot_data_points = Ntot_data_points,
                        ART_data_points = ART_data_points_lazymcmc,
                        PrEP_fitting = PrEP_fitting)



    # lik = likelihood_list[[1]]

    return(lik)
    # return(par_names)
  }
  return(likelihood_func)
}


parTab <- read.csv("parTab.csv")

# f <- create_lik(parTab, NULL, NULL)
# f(unlist(pars))
f <- create_lik(parTab)
f(parTab$values)


devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)






before = Sys.time()
mcmcPars <- c("iterations"=4000,popt=0.44,opt_freq=1000,thin=10,burnin=0,adaptive_period=4000,save_block=100)
Rprof(tmp <- tempfile())
res <- run_MCMC(parTab,data,mcmcPars,"test",create_lik,NULL,NULL,0.1)
Rprof()
summaryRprof(tmp)
after = Sys.time()

after - before


chain <- read.csv(res$file)
plot(coda::as.mcmc(chain))




mcmc_results = read.csv("test_univariate_chain_2505.csv")


plot(mcmc_results$lnlike)

#from 100 onwards?

mcmc_results_without_burnin = mcmc_results[-c(1:100),]


test_set = mcmc_results_without_burnin[45,]


names(test_set)

ranges_post = ranges
for(i in 1:length(ranges[,1])) {
  ranges_post[i,1] = as.numeric(test_set[rownames(ranges)[i]])
  ranges_post[i,2] = ranges_post[i,1]
}


parameters = cotonou::lhs_parameters(number_simulations, set_pars = best_set, Ncat = 9, time = time,
                                     ranges = ranges_post, par_seq = par_seq, condom_seq = condom_seq, groups_seq = groups_seq, years_seq = years_seq)



#ranges = NULL needs to be changed to incorporate those that vary!

# run model
# res = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)


# res = cotonou::return_outputs(parameters[[1]], cotonou::main_model, time = time, outputs = CEA_outputs)

res_best_runs = lapply(parameters, cotonou::return_outputs, cotonou::main_model, time = time, outputs = CEA_outputs)

# res_best_runs[[1]] = res

# ignore these ######################################
# frac_ProFSW_F = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N), function(x) (x[,1]/(x[,1] + x[,2] + x[,3] + x[,4] + x[,7])))), 2, cotonou::quantile_95)))
# frac_ProFSW_F = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N/x$frac_F), function(x) x[,1])), 2, cotonou::quantile_95)))



annual_client_volume_pro_FSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x)  {x$c_comm[,1]}))))
colnames(annual_client_volume_pro_FSW) = c("time", as.character(seq(1, (length(annual_client_volume_pro_FSW[1,])-1))))
annual_client_volume_pro_FSW_melted = reshape2::melt(annual_client_volume_pro_FSW, id.vars = c("time"))



N_Pro_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1)])}))))
colnames(N_Pro_FSW) = c("time", as.character(seq(1, (length(N_Pro_FSW[1,])-1))))
N_Pro_FSW_melted = reshape2::melt(N_Pro_FSW, id.vars = c("time"))
colnames(N_Pro_FSW_melted) = c("time", "run", "value")
N_Low_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(2)])}))))
colnames(N_Low_FSW) = c("time", as.character(seq(1, (length(N_Low_FSW[1,])-1))))
N_Low_FSW_melted = reshape2::melt(N_Low_FSW, id.vars = c("time"))
colnames(N_Low_FSW_melted) = c("time", "run", "value")
N_Client = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(5)])}))))
colnames(N_Client) = c("time", as.character(seq(1, (length(N_Client[1,])-1))))
N_Client_melted = reshape2::melt(N_Client, id.vars = c("time"))
colnames(N_Client_melted) = c("time", "run", "value")



Fraction_F =   do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1, 2, 3, 4, 7)])}), rowSums)) /
  (do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(5, 6)])}), rowSums)) + do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1, 2, 3, 4, 7)])}), rowSums)))


lambda_sum_0_ProFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,1]))))
lambda_sum_0_LowFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,2]))))
lambda_sum_0_GPF = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,3]))))
lambda_sum_0_FormerFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,4]))))
lambda_sum_0_Client = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,5]))))
lambda_sum_0_GPM = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,6]))))
lambda_sum_0_Virgin_Female = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,7]))))
lambda_sum_0_Virgin_Male = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,8]))))
lambda_sum_0_Former_FSW_Outside = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,9]))))




lambda_sum_0_indiv = rbind(lambda_sum_0_ProFSW, lambda_sum_0_LowFSW, lambda_sum_0_GPF,
                           lambda_sum_0_FormerFSW, lambda_sum_0_Client, lambda_sum_0_GPM)


lambda_sum_0_indiv = data.frame(time, rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM"), each = length(time)), lambda_sum_0_indiv)


colnames(lambda_sum_0_indiv) = c("time", "variable", as.character(seq(1, length(lambda_sum_0_ProFSW[1,]))))
lambda_sum_0_indiv_melted = reshape2::melt(lambda_sum_0_indiv, id.vars = c("time", "variable"))
colnames(lambda_sum_0_indiv_melted) = c("time", "variable", "run", "value")
lambda_sum_0_indiv_melted$variable = factor(lambda_sum_0_indiv_melted$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM"))






frac_ProFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,1])), 2, cotonou::quantile_95)))
frac_LowFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,2])), 2, cotonou::quantile_95)))
frac_GPF = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,3])), 2, cotonou::quantile_95)))
frac_FormerFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,4])), 2, cotonou::quantile_95)))
frac_Client = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,5])), 2, cotonou::quantile_95)))
frac_GPM = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,6])), 2, cotonou::quantile_95)))
frac_Virgin_Female = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,7])), 2, cotonou::quantile_95)))
frac_Virgin_Male = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,8])), 2, cotonou::quantile_95)))
frac_Former_FSW_Outside = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,9])), 2, cotonou::quantile_95)))
frac_Active_FSW = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) {100*(x$frac_N[,1] + x$frac_N[,2])})), 2, cotonou::quantile_95)))
Ratio_Low_Pro = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) {x$frac_N[,2]/ x$frac_N[,1]})), 2, cotonou::quantile_95)))

frac = rbind(frac_ProFSW, frac_LowFSW, frac_GPF, frac_FormerFSW, frac_Client, frac_GPM, frac_Virgin_Female, frac_Virgin_Male, frac_Former_FSW_Outside, frac_Active_FSW, Ratio_Low_Pro)
frac = data.frame(frac, group = rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"), each = length(time)))
colnames(frac) = c("time", "Lower", "Median", "Upper", "variable")
frac$variable = factor(frac$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"))

prev_FSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_FSW)), 2, cotonou::quantile_95))
prev_LowFSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_LowFSW)), 2, cotonou::quantile_95))
prev_client = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_client)), 2, cotonou::quantile_95))
prev_women = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_women)), 2, cotonou::quantile_95))
prev_men = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_men)), 2, cotonou::quantile_95))
prev = rbind(prev_FSW, prev_LowFSW, prev_client, prev_women, prev_men)
prev = data.frame(time, prev, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)))
colnames(prev) = c("time", "Lower", "Median", "Upper", "variable")
prev$variable = factor(prev$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))


prev_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_FSW)))
prev_LowFSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_LowFSW)))
prev_client_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_client)))
prev_women_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_women)))
prev_men_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_men)))
prev_indiv = rbind(prev_FSW_indiv, prev_LowFSW_indiv, prev_client_indiv, prev_women_indiv, prev_men_indiv)

prev_indiv = data.frame(time, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)), prev_indiv)


colnames(prev_indiv) = c("time", "variable", as.character(seq(1, length(prev_FSW_indiv[1,]))))

prev_indiv_melted = reshape2::melt(prev_indiv, id.vars = c("time", "variable"))

colnames(prev_indiv_melted) = c("time", "variable", "run", "value")

prev_indiv_melted$variable = factor(prev_indiv_melted$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))




Ntot = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$Ntot)), 2, cotonou::quantile_95)))
colnames(Ntot) = c("time", "Lower", "Median", "Upper")


ART_coverage_women = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_women)), 2, cotonou::quantile_95))
ART_coverage_men = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_men)), 2, cotonou::quantile_95))
ART_coverage_FSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_FSW)), 2, cotonou::quantile_95))
ART_coverage_all = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_all)), 2, cotonou::quantile_95))
ART_coverage = rbind(ART_coverage_women, ART_coverage_men, ART_coverage_FSW, ART_coverage_all)
ART_coverage = data.frame(time, ART_coverage, rep(c("Women", "Men", "Pro FSW", "All"), each = length(time)))
colnames(ART_coverage) = c("time", "Lower", "Median", "Upper", "variable")
ART_coverage$variable = factor(ART_coverage$variable, levels = c("Pro FSW", "Women", "Men", "All"))
ART_coverage = ART_coverage[ART_coverage$variable == "All" | ART_coverage$variable == "Pro FSW",]



ART_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_FSW)))
ART_women_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_women)))
ART_men_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_men)))

ART_all_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_all)))


ART_indiv = rbind(ART_FSW_indiv, ART_all_indiv)

ART_indiv = data.frame(time, rep(c("Pro FSW", "All"), each = length(time)), ART_indiv)


colnames(ART_indiv) = c("time", "variable", as.character(seq(1, length(ART_FSW_indiv[1,]))))

ART_indiv_melted = reshape2::melt(ART_indiv, id.vars = c("time", "variable"))

colnames(ART_indiv_melted) = c("time", "variable", "run", "value")

ART_indiv_melted$variable = factor(ART_indiv_melted$variable, levels = c("Pro FSW", "All"))


# N of FSW on ART


N_ART_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$HIV_positive_On_ART[,1])))


N_ART_FSW_indiv = data.frame(time, rep(c("Pro FSW"), each = length(time)), N_ART_FSW_indiv)


colnames(N_ART_FSW_indiv) = c("time", "variable", as.character(seq(1, length(N_ART_FSW_indiv[1,])-2)))

N_ART_FSW_indiv_melted = reshape2::melt(N_ART_FSW_indiv, id.vars = c("time", "variable"))

colnames(N_ART_FSW_indiv_melted) = c("time", "variable", "run", "value")



# N on ART and N off ART diagnosed
Diagnosed_Off_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I22 + x$I23 + x$I24 + x$I25}), function(x) x[,1]))))
Diagnosed_On_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I32 + x$I33 + x$I34 + x$I35}), function(x) x[,1]))))
Diagnosed_Dropout_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I42 + x$I43 + x$I44 + x$I45}), function(x) x[,1]))))
Diagnosed_ALL_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I22 + x$I23 + x$I24 + x$I25 + x$I32 + x$I33 + x$I34 + x$I35 + x$I42 + x$I43 + x$I44 + x$I45}), function(x) x[,1]))))

Diagnosed_FSW = data.frame(rbind(Diagnosed_Off_ART_FSW, Diagnosed_On_ART_FSW, Diagnosed_Dropout_ART_FSW, Diagnosed_ALL_FSW),
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout", "All diagnosed"), each = length(time)))
colnames(Diagnosed_FSW) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_FSW[1,])-1))), "variable")

Diagnosed_FSW_melted = reshape2::melt(Diagnosed_FSW, id.vars = c("time", "variable"))

Diagnosed_FSW_melted$group = "FSW"

colnames(Diagnosed_FSW_melted) = c("time", "variable", "run", "value", "group")


Diagnosed_Off_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))
Diagnosed_On_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))
Diagnosed_Dropout_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))

Diagnosed_All = data.frame(rbind(Diagnosed_Off_ART_All, Diagnosed_On_ART_All, Diagnosed_Dropout_ART_All),
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_All) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_All[1,])-1))), "variable")

Diagnosed_All_melted = reshape2::melt(Diagnosed_All, id.vars = c("time", "variable"))

Diagnosed_All_melted$group = "All"


colnames(Diagnosed_All_melted) = c("time", "variable", "run", "value", "group")


Diagnosed = rbind(Diagnosed_FSW_melted, Diagnosed_All_melted)



Diagnosed_Off_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(5, 6)])}), rowSums))))
Diagnosed_On_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(5, 6)])}), rowSums))))
Diagnosed_Dropout_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(5, 6)])}), rowSums))))

Diagnosed_Men = data.frame(rbind(Diagnosed_Off_ART_Men, Diagnosed_On_ART_Men, Diagnosed_Dropout_ART_Men),
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_Men) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_Men[1,])-1))), "variable")

Diagnosed_Men_melted = reshape2::melt(Diagnosed_Men, id.vars = c("time", "variable"))

Diagnosed_Men_melted$group = "Men"


colnames(Diagnosed_Men_melted) = c("time", "variable", "run", "value", "group")



Diagnosed_Off_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))
Diagnosed_On_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))
Diagnosed_Dropout_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))

Diagnosed_Women = data.frame(rbind(Diagnosed_Off_ART_Women, Diagnosed_On_ART_Women, Diagnosed_Dropout_ART_Women),
                             rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_Women) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_Women[1,])-1))), "variable")

Diagnosed_Women_melted = reshape2::melt(Diagnosed_Women, id.vars = c("time", "variable"))

Diagnosed_Women_melted$group = "Women"


colnames(Diagnosed_Women_melted) = c("time", "variable", "run", "value", "group")

Diagnosed_Women_Men = rbind(Diagnosed_Women_melted, Diagnosed_Men_melted)

Diagnosed_Women_Men_ratio = data.frame(Diagnosed_Women_melted[,c("time", "variable", "run")],
                                       value = Diagnosed_Women_melted$value/Diagnosed_Men_melted$value)


HIV_deaths = data.frame(time[-1], t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x) x$cumuHIVDeaths), rowSums), diff) )))

colnames(HIV_deaths) = c("time", as.character(seq(1, (length(HIV_deaths[1,])-1))))

HIV_deaths_melted = reshape2::melt(HIV_deaths, id.vars = c("time"))







pc_OfWomen_ProFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1]))/
                                         (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7]))


                                         )))

pc_OfWomen_LowFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2]))/
                                         (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))



pc_OfWomen_Active_FSW = data.frame(time, t(100*(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) + do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])))/
                                             (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))


pc_OfWomen_GPF = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3]))/
                                      (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))

pc_OfWomen_FormerFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))/
                                            (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))



pc_OfMen_Client = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5]))/
                                       (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                          do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                          do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))


pc_OfMen_GPM = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))/
                                    (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                       do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                       do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))


pc_OfWomen_VF = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7]))/
                                     (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) +
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))


pc_OfMen_VM = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8]))/
                                   (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                      do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                      do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))

Ratio_Low_Pro = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {x$frac_N[,2]/ x$frac_N[,1]}))))


frac_by_gender = rbind(pc_OfWomen_ProFSW,
                       pc_OfWomen_LowFSW,
                       pc_OfWomen_GPF,
                       pc_OfWomen_FormerFSW,
                       pc_OfMen_Client,
                       pc_OfMen_GPM,
                       pc_OfWomen_VF,
                       pc_OfMen_VM,
                       pc_OfWomen_Active_FSW,
                       Ratio_Low_Pro
)
frac_by_gender = data.frame(frac_by_gender, group = rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"), each = length(time)))
colnames(frac_by_gender) = c("time",as.character(seq(1, (length(pc_OfWomen_ProFSW[1,])-1))),  "variable")
frac_by_gender$variable = factor(frac_by_gender$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"))


frac_by_gender_melted = reshape2::melt(frac_by_gender, id.vars = c("time", "variable"))

colnames(frac_by_gender_melted) = c("time", "variable", "run", "value")

colnames(pc_OfWomen_ProFSW) = c("time", as.character(seq(1, (length(pc_OfWomen_ProFSW[1,])-1))))
colnames(pc_OfWomen_LowFSW) = c("time", as.character(seq(1, (length(pc_OfWomen_LowFSW[1,])-1))))
colnames(pc_OfWomen_Active_FSW) = c("time", as.character(seq(1, (length(pc_OfWomen_Active_FSW[1,])-1))))
colnames(pc_OfMen_Client) = c("time", as.character(seq(1, (length(pc_OfMen_Client[1,])-1))))


colnames(pc_OfWomen_VF) = c("time", as.character(seq(1, (length(pc_OfWomen_VF[1,])-1))))
colnames(pc_OfMen_VM) = c("time", as.character(seq(1, (length(pc_OfMen_VM[1,])-1))))

pc_OfWomen_ProFSW_melted = reshape2::melt(pc_OfWomen_ProFSW, id.vars = c("time"))
pc_OfWomen_LowFSW_melted = reshape2::melt(pc_OfWomen_LowFSW, id.vars = c("time"))
pc_OfWomen_Active_FSW_melted = reshape2::melt(pc_OfWomen_Active_FSW, id.vars = c("time"))
pc_OfMen_Client_melted = reshape2::melt(pc_OfMen_Client, id.vars = c("time"))



pc_OfWomen_VF_melted = reshape2::melt(pc_OfWomen_VF, id.vars = c("time"))
pc_OfMen_VM_melted = reshape2::melt(pc_OfMen_VM, id.vars = c("time"))







condom_Pro_FSW_comm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_comm[,1,][,5]}))))
condom_Pro_FSW_noncomm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,1,][,5]}))))

condom_Pro_FSW = data.frame(time, rbind(condom_Pro_FSW_comm, condom_Pro_FSW_noncomm), rep(c("Commercial", "Non commercial"), each = length(time)))
colnames(condom_Pro_FSW) = c("time", as.character(seq(1, (length(condom_Pro_FSW[1,])-2))), "variable")

condom_Pro_FSW_melted = reshape2::melt(condom_Pro_FSW, id.vars = c("time", "variable"))
colnames(condom_Pro_FSW_melted) = c("time", "variable", "run", "value")


condom_Low_FSW_comm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_comm[,2,][,5]}))))
condom_Low_FSW_noncomm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,2,][,5]}))))

condom_Low_FSW = data.frame(time, rbind(condom_Low_FSW_comm, condom_Low_FSW_noncomm), rep(c("Commercial", "Non commercial"), each = length(time)))
colnames(condom_Low_FSW) = c("time", as.character(seq(1, (length(condom_Low_FSW[1,])-2))), "variable")

condom_Low_FSW_melted = reshape2::melt(condom_Low_FSW, id.vars = c("time", "variable"))
colnames(condom_Low_FSW_melted) = c("time", "variable", "run", "value")

condom_GPF_noncomm = data.frame(time, t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,3,][,6]})))))
colnames(condom_GPF_noncomm) = c("time", as.character(seq(1, (length(condom_GPF_noncomm[1,])-1))))
condom_GPF_noncomm_melted = reshape2::melt(condom_GPF_noncomm, id.vars = c("time"))

condom_GPM_noncomm = data.frame(time, t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,6,][,3]})))))
colnames(condom_GPM_noncomm) = c("time", as.character(seq(1, (length(condom_GPM_noncomm[1,])-1))))
condom_GPM_noncomm_melted = reshape2::melt(condom_GPM_noncomm, id.vars = c("time"))






testing_rate_ratio_F_M = data.frame(time, t((do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$testing_prob}), function(x) {return(x[,c(3)])}))/do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$testing_prob}), function(x) {return(x[,c(6)])})))))

colnames(testing_rate_ratio_F_M) = c("time", as.character(seq(1, (length(testing_rate_ratio_F_M[1,])-1))))

testing_rate_ratio_F_M_melted = reshape2::melt(testing_rate_ratio_F_M, id.vars = c("time"))

colnames(testing_rate_ratio_F_M_melted) = c("time", "run", "value")

Diagnosed_Women_Men_ratio_On_ART = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed On ART",]


testing_rate_ratio_F_M_Women_Men_ratio_On_ART = data.frame(x=testing_rate_ratio_F_M_melted$value, y=Diagnosed_Women_Men_ratio_On_ART$value)


pfFSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x)  {x$pfFSW[,1]}))))
colnames(pfFSW) = c("time", as.character(seq(1, (length(pfFSW[1,])-1))))
pfFSW_melted = reshape2::melt(pfFSW, id.vars = c("time"))
colnames(pfFSW_melted) = c("time", "run", "value")



ART_inits = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTinitiations))}))))
colnames(ART_inits) = c("time", as.character(seq(1, (length(ART_inits[1,])-1))))
ART_inits_melted = reshape2::melt(ART_inits, id.vars = c("time"))
colnames(ART_inits_melted) = c("time", "run", "value")



# Ratio female:male of rate of initiation per infected = 8.0/4.68 = 1.7
# to compare to the data, I am taking the number of initiations that years as a fraction of the total nubmer of PLHIV that year



# DENOMINTOR = ALL POSTIIVES
female_ART_init_rate_from_all_HIV_pos = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)])[-length(time)] # number of HIV positive females at the beginning of the year
  return(a/b)
}))))


male_ART_init_rate_from_all_HIV_pos = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  c = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of male ART initiations per year
  d = rowSums(x$HIV_positive[,c(5,6)])[-length(time)] # number of HIV positive males at the beginning of the year
  return(c/d)
}))))







# DENOMINTOR = ALL POSTIIVES NOT ON ART

female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)] - x$HIV_positive_On_ART[,c(1,2,3,4)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return(a/b)
}))))
# female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART



male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(5,6)] - x$HIV_positive_On_ART[,c(5,6)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return(a/b)
}))))
# male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART





# FEMALE TO MALE RATIO OF INITIATION RATES WHERE DENOMINATOR IS ART POSITIVE NOT ON ART

FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART = data.frame(time = time[-length(time)], t(do.call(rbind, lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)] - x$HIV_positive_On_ART[,c(1,2,3,4)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  c = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of female ART initiations per year
  d = rowSums(x$HIV_positive[,c(5,6)] - x$HIV_positive_On_ART[,c(5,6)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return((a/b)/(c/d))
}))))

colnames(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART) = c("time", as.character(seq(1, (length(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART[1,])-1))))

FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART_melted = reshape2::melt(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART, id.vars = c("time"))



ART_inits_ratio_F_over_M =  data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) /
    diff(rowSums(x$cumuARTinitiations[,c(5,6)]))}))))
ART_inits_ratio_F_over_M_melted = reshape2::melt(ART_inits_ratio_F_over_M, id.vars = c("time"))




ART_inits_cumu = data.frame(time = time, t(do.call(rbind,lapply(res_best_runs, function(x)  {rowSums(x$cumuARTinitiations)}))))
colnames(ART_inits_cumu) = c("time", as.character(seq(1, (length(ART_inits_cumu[1,])-1))))
ART_inits_cumu_melted = reshape2::melt(ART_inits_cumu, id.vars = c("time"))
colnames(ART_inits_cumu_melted) = c("time", "run", "value")



ART_REinits_cumu = data.frame(time = time, t(do.call(rbind,lapply(res_best_runs, function(x)  {rowSums(x$cumuARTREinitiations)}))))
colnames(ART_REinits_cumu) = c("time", as.character(seq(1, (length(ART_REinits_cumu[1,])-1))))
ART_REinits_cumu_melted = reshape2::melt(ART_REinits_cumu, id.vars = c("time"))
colnames(ART_REinits_cumu_melted) = c("time", "run", "value")



ART_REinits = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTREinitiations))}))))
colnames(ART_REinits) = c("time", as.character(seq(1, (length(ART_REinits[1,])-1))))
ART_REinits_melted = reshape2::melt(ART_REinits, id.vars = c("time"))
colnames(ART_REinits_melted) = c("time", "run", "value")




frac_ART_inits_REinits = data.frame(time = time[-length(time)])
for(i in 2:(length(res_best_runs)+1))
{
  frac_ART_inits_REinits[,i] = ART_REinits[,i]/(ART_inits[,i] + ART_REinits[,i]) * 100
}
colnames(frac_ART_inits_REinits) = c("time", as.character(seq(1, (length(frac_ART_inits_REinits[1,])-1))))
frac_ART_inits_REinits_melted = reshape2::melt(frac_ART_inits_REinits, id.vars = c("time"))
colnames(frac_ART_inits_REinits_melted) = c("time", "run", "value")




ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex = data.frame(
  Women_2014 = unlist(female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2014),c(2:(length(res_best_runs)+1))]),
  Women_2015 = unlist(female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2015),c(2:(length(res_best_runs)+1))]),

  Men_2014 = unlist(male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2014),c(2:(length(res_best_runs)+1))]),
  Men_2015 = unlist(male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2015),c(2:(length(res_best_runs)+1))])
)


ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted = reshape2::melt(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex)

levels(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted) = c("Women_2014","Women_2015", "Men_2014","Men_2015")


###################






require(ggplot2)




# want plots on number of initiations of ART, dropout, reinitiations! compare to data points from the reports!

# need to do difference between years for cumulative art inits etc


prev_axes = data.frame(variable = c(rep("Pro FSW", 2),
                                    rep("Clients", 2),
                                    rep("Women", 2),
                                    rep("Men", 2),
                                    rep("Low-level FSW", 2)),
                       time = c(rep(c(1986, 2025), 5)),
                       value = c(0, 70, 0, 70, 0, 15, 0, 15, 0, 70)
)
prev_points_80s = prev_points_all[c(1,2,3),]


prev_points[prev_points$time == "2015", "lower"][1] = 13.79

# plot prevalence in each group indiv runs
g1=ggplot() + geom_line(data = prev_indiv_melted, aes(x = time, y = value, factor = variable, factor = run), alpha = 0.3) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "prevalence (%)") +
  geom_point(data = prev_points, aes(x = time, y = value))+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper))+
  geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "red")+
  geom_blank(data = prev_axes, aes(x = time, y = value))+
  theme(text = element_text(size=20))

g1


ART_data_points_with_numbers_FSW = data.frame(time = c(2012, 2014, 2015, 2016, 2017),
                                              Lower = c(27,	34,	42, 53,	73),
                                              Upper = c(37, 46, 56, 83, 107),
                                              variable = rep("Pro FSW", 5),
                                              datatype = c("data", "data", "data", "counterfactual", "counterfactual"))

ART_data_points_with_numbers_FSW_points = data.frame(time = c(2012, 2014, 2015, 2016, 2017),
                                                     value = c(32,	40,	49, 53,	73),
                                                     variable = rep("Pro FSW", 5),
                                                     datatype = c("data", "data", "data", "counterfactual", "counterfactual"))


g7=ggplot() + geom_line(data = N_ART_FSW_indiv_melted, aes(x = time, y = value, factor = variable, factor = run), alpha = 1) + ggtitle("Number of professional FSW on ART") +theme_bw() +
  geom_errorbar(data = ART_data_points_with_numbers_FSW, aes(x = time, ymin = Lower, ymax = Upper, col = datatype)) + labs(x = "Year") +
  scale_x_continuous(breaks = seq(2000, 2018, 1), limits = c(2000,2018))+
  labs(y = "")+
  geom_point(data = ART_data_points_with_numbers_FSW_points, aes(x = time, y = value, col = datatype))+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2008.8, y = 70, size = 7,col = "black", label = "Data on numbers of FSW file active \n from Luc +- 15% \n Fitting to 2015 estimate only")

g7




frac_by_gender_discard_points_no_FSW_LB = data.frame(variable = c("Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"),
                                                     min = c(0, 0.074, 0.07896475, 0.07039551, 0.0048, .01),
                                                     max = c(0.0143/2, 0.3 , 0.2, 0.17,  0.0143, .05))



# plot fraction in each group
g2=ggplot(frac_by_gender_melted) + geom_line(aes(x = time, y = value, factor = run))  +
  theme_bw() + labs(y = "Percent in each group in their respective gender (%)") +  facet_wrap(~variable, scales = "free") +
  # geom_point(data = frac_N_data_points, aes(x = time, y = point), size = I(2), color = "red", shape = 15) +
  geom_hline(data = frac_by_gender_discard_points_no_FSW_LB, aes(yintercept = 100*min), size = I(0.5), color = "red", linetype = 1, alpha = 0.7) +
  geom_hline(data = frac_by_gender_discard_points_no_FSW_LB, aes(yintercept = 100*max), size = I(0.5), color = "red", linetype = 1, alpha = 0.7)+
  theme(text = element_text(size=20))
g2






g3=ggplot() + geom_line(data = pc_OfWomen_ProFSW_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + labs(y = "% of all women 15-59 that are pro FSW")+
  theme(text = element_text(size=20)) + geom_hline(yintercept = 0.24,col = "orange", linetype="dashed") + geom_hline(yintercept = 0.72, linetype="dashed",col = "red")+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2002.8, y = 0.7, size = 7,col = "red", label = "Fitted to below this line")+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2003.8, y = 0.3, size = 7,col = "orange", label = "Original estimate lower bound ")
g3






# g3b=ggplot() + geom_line(data = pc_OfWomen_ProFSW_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + labs(y = "% of all women 15-59 that are pro FSW")+
#   theme(text = element_text(size=20)) + geom_hline(yintercept = 0.24,col = "orange", linetype="dashed") + geom_hline(yintercept = 0.72, linetype="dashed",col = "red")+
#   theme(text = element_text(size=24),
#         legend.text=element_text(size=18),
#         legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2002.8, y = 0.7, size = 7,col = "red", label = "Fitted to below this line")+
#   theme(text = element_text(size=24),
#         legend.text=element_text(size=18),
#         legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2003.8, y = 0.3, size = 7,col = "orange", label = "Original estimate lower bound ")
# g3b


g4=ggplot() + geom_line(data = pc_OfWomen_LowFSW_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all women 15-59 that are low level FSW")+
  theme(text = element_text(size=20))


g4b = ggplot() + geom_line(data = pc_OfWomen_VF_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all women 15-59 that are virgins")+
  theme(text = element_text(size=20))+ geom_hline(yintercept = c(7.9, 20), col = "red") +
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2010, y = 19, size = 7,col = "red", label = "Fitting")

g4b


g5=ggplot() + geom_line(data = pc_OfMen_Client_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + labs(y = "% of all men 15-59 that are clients")+
  theme(text = element_text(size=20)) + geom_hline(yintercept = c(7.4, 30), col = "red") +
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2003.8, y = 28, size = 7,col = "red", label = "Fitting")

g5

g5b = ggplot() + geom_line(data = pc_OfMen_VM_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all men 15-59 that are virgins")+
  theme(text = element_text(size=20))+ geom_hline(yintercept = c(7, 17), col = "red")+
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2010, y = 15, size = 7,col = "red", label = "Fitting")

g5b
# # plot prevalence in each group
# ggplot() + geom_line(data = prev, aes(x = time, y = Median))+ geom_ribbon(data = prev, aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "prevalence (%)") +
#   geom_point(data = prev_points, aes(x = time, y = value))+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper)) +
#   geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "red")+
#   geom_blank(data = prev_axes, aes(x = time, y = value))
#




# plot total population size
g6=ggplot(Ntot) + geom_line(aes(x = time, y = Median)) + geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) +
  theme_bw() + labs(y = "Total population size of Grand Cotonou") +
  geom_point(data = Ntot_data_points, aes(x = time, y = point, color = colour), size = I(2), shape = 15) + geom_errorbar(data = Ntot_data_points, aes(x = time, ymax = upper, ymin = lower, color = colour), width = 2)+
  theme(text = element_text(size=20))


# ggplot(ART_coverage) +
#   geom_line(aes(x = time, y = Median))+ geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() +
#   facet_wrap(~variable) + labs(y = "ART coverage ") +
#   geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower, ymax = Upper), colour = "darkred")





ART_data_points = data.frame(time = c(2011, 2012, 2013, 2014, 2015, 2016, 2017,

                                      2012, 2014, 2015, 2016, 2017
),
Lower = c(0.33, 0.42, 0.44, 0.39, 0.44, 0.51, 0.57,

          0.09,	0.14,	0.2,	0.22,	0.3

),
Upper = c(0.52, 0.66, 0.69, 0.61, 0.69, 0.8, 0.8,
          0.13, 0.2, 0.28, 0.3, 0.42

),
variable = c("All", "All", "All", "All", "All", "All", "All",
             "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))


ART_data_points$Data = c(rep("UNAIDS (fitted 2011 & 2017 only)", 7), rep("File active FSW (X-validation)", 5))

art_cov = data.frame(time = 2017, variable = "All", Upper = 0.91, Lower = 0.6, Data = "Audit report 2017 (x-validation)")

ART_data_points = rbind(ART_data_points, art_cov)

ART_text = data.frame(time = 2010, variable = "Pro FSW", value = 50)




g8=ggplot() + geom_line(data = ART_indiv_melted, aes(x = time, y = value*100, factor = variable, factor = run), alpha = 0.3) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "ART coverage (%)") +
  geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower*100, ymax = Upper*100, col = Data), size = I(1))+
  theme(text = element_text(size=20)) + theme(legend.position = "top") + geom_text(data = ART_text, aes(x = time, y = value, factor = variable),
                                                                                   label = "2016 & 2017 are counterfactuals \n NOT including FSW on TasP/PrEP", size = 5)
g8


g9=ggplot() + geom_line(data = lambda_sum_0_indiv_melted[lambda_sum_0_indiv_melted$time > 2009,], aes(x = time, y = value*100, factor = variable, factor = run), alpha = 0.3) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "Force of infection on susceptibles off PrEP (%)")+
  theme(text = element_text(size=20))


Diagnosed_FSW_melted$col = Diagnosed_FSW_melted$variable
Diagnosed_FSW_melted$col = ifelse(Diagnosed_FSW_melted$col == "Diagnosed Off ART", "#7fc97f",
                                  ifelse(Diagnosed_FSW_melted$col == "Diagnosed On ART", "#beaed4",
                                         ifelse(Diagnosed_FSW_melted$col == "Dropout", "#fdc086", 0)))


# 2012 FSW IBBA ever tested 0.75 - 0.83
# EQS mapping 2012 N = 889 - 1391
# 2012 prevalence = 27.4% [0.23 - 0.322]
# lower N diagnosed = 0.75 * 889 * 0.23 =
# upper N diagnosed = 0.83 * 1391 * 0.322 =


# bounds of infected ever tested = 304 * 0.75, 304 * 0.83 = 228 - 252
pc_diagnosed = data.frame(time = 2012, Lower = 153, Upper = 372, variable = "Pro FSW")

g10=ggplot() + geom_line(data = Diagnosed_FSW_melted, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + theme_bw() + labs(y = "Numbers of pro FSW in Grand Cotonou")+
  scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), #max(Diagnosed_FSW_melted$value), 20) )+
                                  380, 40) )+
  theme(text = element_text(size=20), legend.position = "top")+
  geom_errorbar(data = ART_data_points_with_numbers_FSW, size=1,aes(x = time, ymin = Lower, ymax = Upper), col ="#4daf4a") +
  labs(x = "Year")+
  geom_errorbar(data = pc_diagnosed, size=1, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  labs(x = "Year")+ scale_colour_brewer(type = "qual", palette = 6) +
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 2003.8, y = 222, size = 7,col = "#e41a1c", label = "Number of FSW * \n prevalence * \n % ever tested \n (cross-validation)")

g10
# #
# ggplot() + geom_line(data = Diagnosed_All_melted, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
#   geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
#   geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
#   geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
#   geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2016, label="Everyone", y=-800))
# #


# I2x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Diagnosed Off ART",]
# I3x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Diagnosed On ART",]
# I4x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Dropout",]


I2xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Diagnosed Off ART",]
I3xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Diagnosed On ART",]
I4xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Dropout",]


I2xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Diagnosed Off ART",]
I3xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Diagnosed On ART",]
I4xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Dropout",]







diag_off_art_auditF = data.frame(time = 2017, Lower = 195, Upper = 209, variable = "Women")
diag_off_art_auditM = data.frame(time = 2017, Lower = 100, Upper = 109, variable = "Men")



g10a = ggplot() + geom_line(data = I2xF, aes(x = time, y = value, factor = variable, factor = run), colour = "black", alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
  geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
  geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
  geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
  geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5))+
  geom_errorbar(data = diag_off_art_auditF, size=1, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  ggtitle("Numbers of all women in Grand Cotonou diagnosed NOT on ART compared to data in red (Pre-ART audit report 2017)")
g10a




g10b = ggplot() + geom_line(data = I2xM, aes(x = time, y = value, factor = variable, factor = run), colour = "black", alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
  geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
  geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
  geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
  geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5))+
  geom_errorbar(data = diag_off_art_auditM, size=1, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  ggtitle("Numbers of all men in Grand Cotonou diagnosed NOT on ART compared to data in red (Pre-ART audit report 2017)")
g10b




I2xF_2017 = data.frame(
  # Women_2016=I2xF[I2xF$time == 2016, "value"],
  Women_2017=I2xF[I2xF$time == 2017, "value"],
  # Men_2016=I2xM[I2xM$time == 2016, "value"],
  Men_2017=I2xM[I2xM$time == 2017, "value"]
)
I2xF_2017_melted = reshape2::melt(I2xF_2017)

g10c = ggplot(I2xF_2017_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  geom_errorbar(aes(x = "Women_2017", ymax = 201, ymin = 215), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "Men_2017", ymax = 108, ymin = 117), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sex", y = "")+ ggtitle("Across the difference runs, the number of women and men in Grand Cotonou \n diagnosed but NOT on ART in 2016 and 2017 vs the data in red for 2017")
g10c

I3xF_2017 = data.frame(Women=I3xF[I3xF$time == 2017, "value"], Men=I3xM[I3xM$time == 2017, "value"])
I3xF_2017_melted = reshape2::melt(I3xF_2017)

g10d = ggplot(I3xF_2017_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  geom_errorbar(aes(x = "Women", ymax = 6155, ymin = 13050), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "Men", ymax = 2369, ymin = 5223), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sex", y = "")+
  ggtitle("Number of women and men in Grand Cotonou diagnosed \n ON ART in 2017 vs the data in red")
g10d

I4xF_2015 = data.frame(Women=I4xF[I4xF$time == 2015, "value"], Men=I4xM[I4xM$time == 2015, "value"])
I4xF_2017_melted = data.frame(variable = "All", value = rowSums(I4xF_2015))




g12=ggplot() + geom_line(data = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed On ART",], aes(x = time, y = value, factor = variable, factor = run), alpha = 1) +
  theme_bw() + labs(y = "")+
  ggtitle("Ratio of Women:Men on ART compared to data (2017) in red from audit 2017")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+1), 4) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2017, y = 2.5), col = "red", size = I(03))+
  # geom_hline(yintercept = c(1,2), col = "orange") +
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))#+ annotate("text", x = 1994, y = 2.5, size = 7,col = "orange")

g12


g10e = ggplot(I4xF_2017_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  # geom_errorbar(aes(x = "Women", ymax = 6155, ymin = 13050), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "All", ymax = 56, ymin = 56), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sex", y = "")+
  ggtitle("Number of women and men in Grand Cotonou in dropout \n from ART in 2017 vs the data in red (Lower bound)")
g10e


g10f = ggplot() + geom_line(data = ART_inits_cumu_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ggtitle("Cumulative ART initiations of all groups combined")

g10f

g10g = ggplot() + geom_line(data = ART_inits_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_errorbar(aes(x = 2015, ymax = 2051, ymin = 1318), col = "red", width = I(0.5), size = I(02))+
  ggtitle("Numbers of ART initiations each year of all groups combined vs. the data in red")


g10g




g10h=ggplot() + geom_line(data = ART_inits_ratio_F_over_M_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+3), 2) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 2.4), col = "red", size = I(03))+
  ggtitle("Ratio of ART initiations women:Men, \n compared with data from audit 2017 in red")
g10h




g10i=ggplot() + geom_line(data = FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+3), 2) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 1.7), col = "red", size = I(03))+
  ggtitle("Ratio of ART initiation rate from HIV+ not on ART Women:Men, \n compared with data from audit 2017 in red")
g10i

g10j = ggplot(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  labs(y = "")+
  geom_point(aes(x = "Women_2015", y = 0.169), col = "red", size = 3)+
  geom_point(aes(x = "Men_2015", y = 0.1), col = "red", size = 3)+
  ggtitle("ART initiation rate from HIV+ not on ART Women & Men in 2014 and 2015, \n compared with data for 2015 in red")
g10j


# Ratio female:male of rate of initiation per infected = 8.0/4.68 = 1.7




g10k = ggplot() + geom_line(data = ART_REinits_cumu_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  # geom_errorbar(aes(x = 2015, ymax = 2051, ymin = 1318), col = "red", width = I(0.5), size = I(02))+
  ggtitle("Cumulative numbers of ART RE initiations each year of all groups combined")


g10k


g10l = ggplot() + geom_line(data = ART_REinits_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 213), size = 5, col = "red") +
  ggtitle("Numbers of ART RE initiations each year of all groups combined vs. the data in red") +
  geom_text(col = "red",aes(x=2002, label="Data from monitoring reports for 2015 \n for Littoral + Atlantique + Oueme", y=300), size = 8)



g10l




g10m = ggplot() + geom_line(data = frac_ART_inits_REinits_melted[frac_ART_inits_REinits_melted$time > 1999,], aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2013.5, y = 8), size = 5, col = "red") +
  geom_point(aes(x = 2014.5, y = 7), size = 5, col = "red") +
  geom_point(aes(x = 2015, y = 7), size = 5, col = "red") +
  geom_point(aes(x = 2015.5, y = 9), size = 5, col = "red") +
  geom_point(aes(x = 2016, y = 12), size = 5, col = "red") +
  geom_point(aes(x = 2017, y = 11), size = 5, col = "red") +
  ggtitle("Percent of ART initiations which are RE initiations (%) vs. the data in red")  +
  geom_text(col = "red",aes(x=2005, label="Data from monitoring reports by semester for BENIN", y=30), size = 8)

g10m
# ggplot() + geom_line(data = I2x, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
#   geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
#   geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
#   geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
#   geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed")




#
#
# ggplot() + geom_line(data = Diagnosed_Women_Men, aes(x = time, y = value, factor = variable, factor = run, colour = variable, linetype = group), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) )
#
# ggplot() + geom_line(data = Diagnosed_Women_Men_ratio, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) +
#   theme_bw() + labs(y = "Ratio of Women:Men")+
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(1, 5, 0.2) )

# ggplot() + geom_line(data = Diagnosed, aes(x = time, y = value, factor = variable, factor = run, colour = group), alpha = 1) + theme_bw() + facet_wrap(~variable, scales = "fixed") + labs(y = "Numbers of all people in Grand Cotonou") +
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) )
#
# #

g11=ggplot() + geom_line(data = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed Off ART",], aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) +
  theme_bw() + labs(y = "Ratio of Women:Men")+
  scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20))



# cross validation to ratio female:male on ART from audit report 2017


HIV_deaths_aidsinfo = data.frame(time = seq(1990, 2016),
                                 lower = c(100,100,100,100,100,100,100,100,200,200,200,500,500,500,500,500,500,500,500,200,200,200,200,200,200,200,100),
                                 upper = c(300,300,300,300,300,400,400,500,600,900,1200,1500,1500,1500,1500,2000,2000,2000,2000,900,900,800,1300,1300,1300,1300,700))

g13=ggplot() + geom_line(data = HIV_deaths_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() +
  labs(y = "HIV deaths") +
  geom_errorbar(data = HIV_deaths_aidsinfo, aes(x = time, ymin = lower, ymax = upper), col = "orange")+
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 1700, size = 7,col = "orange",
                                                     label = "Cross validation")

g13



# EQS 2012 mapping 889 – 1391
EQS_mapping = data.frame(time = 2012, lower = 889, upper = 1391)
g14=ggplot() + geom_line(data = N_Pro_FSW_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "Number of Pro FSW") +
  theme(text = element_text(size=20)) +
  geom_errorbar(data = EQS_mapping, aes(x = time, ymin = lower, ymax = upper), col = "red")+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 1200, size = 7,col = "red",
                                                     label = "Fitting to EQS mapping,\n  CI given by Michel")

g14

N_low_level = data.frame(time = 2012, lower = 2000, upper = 3000)

g15=ggplot() + geom_line(data = N_Low_FSW_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "Number of Low level FSW") +
  theme(text = element_text(size=20)) +
  geom_errorbar(data = N_low_level, aes(x = time, ymin = lower, ymax = upper), col = "orange")+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 2750, size = 7,col = "orange",
                                                     label = "Cross validating")


low_pro_ratio = N_Low_FSW_melted
low_pro_ratio$value = N_Low_FSW_melted$value / N_Pro_FSW_melted$value

g22 = ggplot() + geom_line(data = low_pro_ratio, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "Ratio of low level FSW to pro FSW") +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 4.5, size = 7,col = "red",
                                                     label = "Fitting") +
  geom_hline(yintercept = c(1, 5), col = "red")


# lowndes 2002 19970

g16=ggplot() + geom_line(data = N_Client_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "Number of Clients")








g17=ggplot() + geom_line(data = condom_Pro_FSW_melted, aes(x = time, y = value, colour = variable, factor = run), alpha = 1) +
  theme_bw() + labs(y = "Condom use of Pro FSW")+
  theme(text = element_text(size=20))

g18=ggplot() + geom_line(data = condom_Low_FSW_melted, aes(x = time, y = value, colour = variable, factor = run), alpha = 1) +
  theme_bw() + labs(y = "Condom use of Low level FSW")+
  theme(text = element_text(size=20))

g19=ggplot() + geom_line(data = condom_GPF_noncomm_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "Condom use of GPF")+
  theme(text = element_text(size=20))

# ggplot() + geom_line(data = condom_GPM_noncomm_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
#   theme_bw() + labs(y = "Condom use of GPM")


#



g20=ggplot() + geom_line(data = annual_client_volume_pro_FSW_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "Annual client volume for Pro FSW")+
  theme(text = element_text(size=20))

g21=ggplot() + geom_line(data = pfFSW_melted, aes(x = time, y = value*100, factor = run), alpha = 1) +
  theme_bw() + labs(y = "HIV prevalence of incoming FSWs from neighbouring countries(%)")+
  theme(text = element_text(size=20))






plotlist = list(g1, g3, g4,g4b, g5,g5b, g6, g7, g8, g9, g10,
                g10a, g10b, g10c, g10d,g12, g10e, g10f, g10g, g10h, g10i, g10j,
                g13, g14, g15, g16, g17, g18, g19, g20, g21, g22)

# outputfolder = "0103"


for(i in 1:length(plotlist))

{
  jpeg(filename=paste0("C:\\Users\\eg1012\\Google Drive\\Imperial\\Cotonou\\Models\\Transmission Model\\Outputs/",
                       batch_folder, "/","all", i,".jpg"), width = 1000, height = 600)
  print(plotlist[i][[1]])
  dev.off()
}


