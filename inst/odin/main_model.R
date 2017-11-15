
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
########################################################## MAIN MODEL ######################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# INDEXING
##############################################################################
#
# Risk group (i)
# 1. Professional FSW
# 2. Low-level FSW
# 3. General population female
# 4. Former FSW in Cotonou
# 5. Clients
# 6. General population male
# 7. Virgin Female
# 8. Virgin Male
# 9. Former FSW in Benin, outside Cotonou (not involved in epidemic, but tracked anyway)

# Age (x)
# 1. 15 - 24
# 2. 25 - 34
# 3. 35 - 59



# NOTES
##############################################################################

# hopefully not dividing by 0 for the balancing
# lambdas should be different for each compartment??!?!?!?!
# I05 no gamma5s anywhere!
# calculate E for all, then fraction for each pop
# if sum(omega) != 1, then problem! - how do I hard code this error?

# each variable must have a derivative part and an initial part

# a and k do not exist here yet

config(include) = "FOI.c"

# ORDINARY DIFFERENTIAL EQUATIONS
##############################################################################

replaceDeaths = user()

# births and prep movement
E0[] = if(replaceDeaths == 1) mu[i] * N[i] + nu * N[i] + alphaItot[i] + epsilon * Ntot * omega[i] - S0[i] * zeta[i] else new_people_in_group[i] - S0[i] * zeta[i] + new_people_in_group_FSW_only[i]
E1a[] = zeta[i] * fPa * S0[i] - psia[i] * S1a[i] - kappaa[i] * S1a[i]
E1b[] = zeta[i] * fPb * S0[i] + psia[i] * S1a[i] - psib[i] * S1b[i]  - kappab[i] * S1b[i]
E1c[] = zeta[i] * fPc * S0[i] + psib[i] * S1b[i] - kappac[i] * S1c[i]
E1d[] = kappaa[i] * S1a[i] + kappab[i] * S1b[i] + kappac[i] * S1c[i]


deriv(S0[]) = E0[i] - S0[i] * lambda_sum_0[i] - S0[i] * mu[i] - S0[i] * nu + rate_move_out[i] * S0[i] + sum(in_S0[i, ])
deriv(S1a[]) = E1a[i] - S1a[i] * lambda_sum_1a[i] - S1a[i] * mu[i] - S1a[i] * nu + rate_move_out[i] * S1a[i] + sum(in_S1a[i, ])
deriv(S1b[]) = E1b[i] - S1b[i] * lambda_sum_1b[i] - S1b[i] * mu[i] - S1b[i] * nu + rate_move_out[i] * S1b[i] + sum(in_S1b[i, ])
deriv(S1c[]) = E1c[i] - S1c[i] * lambda_sum_1c[i] - S1c[i] * mu[i] - S1c[i] * nu + rate_move_out[i] * S1c[i] + sum(in_S1c[i, ])
deriv(S1d[]) = E1d[i] - S1d[i] * lambda_sum_1d[i] - S1d[i] * mu[i] - S1d[i] * nu + rate_move_out[i] * S1d[i] + sum(in_S1d[i, ])

#primary infection
deriv(I01[]) = S0[i] * lambda_sum_0[i] + S1d[i] * lambda_sum_1d[i] - I01[i] * (gamma01[i] + tau[i] + tau_intervention[i] + alpha01[i] + mu[i] + nu) + rate_move_out[i] * I01[i] + sum(in_I01[i, ]) +
  kappa1[i] * I11[i]

deriv(I11[]) = S1a[i] * lambda_sum_1a[i] + S1b[i] * lambda_sum_1b[i] + S1c[i] * lambda_sum_1c[i] -
  I11[i] * (gamma11[i] + test_rate_prep[i] + kappa1[i] + alpha11[i] + mu[i] + nu) + rate_move_out[i] * I11[i] + sum(in_I11[i, ]) # test_rate_prep is tau^1,1


#chronic
deriv(I02[]) = gamma01[i] * I01[i] + gamma11[i] * I11[i] - I02[i] * (gamma02[i] + tau[i] + tau_intervention[i] + alpha02[i] + mu[i] + nu) + rate_move_out[i] * I02[i] + sum(in_I02[i, ])
deriv(I03[]) = gamma02[i] * I02[i] - I03[i] * (gamma03[i] + tau[i] + tau_intervention[i] + alpha03[i] + mu[i] + nu) + rate_move_out[i] * I03[i] + sum(in_I03[i, ])
deriv(I04[]) = gamma03[i] * I03[i] - I04[i] * (gamma04[i] + tau[i] + tau_intervention[i] + alpha04[i] + mu[i] + nu) + rate_move_out[i] * I04[i] + sum(in_I04[i, ])
deriv(I05[]) = gamma04[i] * I04[i] - I05[i] * (RR_test_CD4200*tau[i] + tau_intervention[i] + alpha05[i] + mu[i] + nu) + rate_move_out[i] * I05[i] + sum(in_I05[i, ])

deriv(I22[]) = (tau[i] + tau_intervention[i]) * I01[i] + test_rate_prep[i] * I11[i] + (tau[i] + tau_intervention[i]) * I02[i] - I22[i] * (gamma22[i] + rho_intervention[i] + rho[i]*ART_eligible_CD4_above_500 + alpha22[i] + mu[i] + nu) + rate_move_out[i] * I22[i] + sum(in_I22[i, ])
deriv(I23[]) = gamma22[i] * I22[i] + (tau[i] + tau_intervention[i]) * I03[i] - I23[i] * (gamma23[i] + rho[i]*ART_eligible_CD4_350_500 + rho_intervention[i] + alpha23[i] + mu[i] + nu) + rate_move_out[i] * I23[i] + sum(in_I23[i, ])
deriv(I24[]) = gamma23[i] * I23[i] + (tau[i] + tau_intervention[i]) * I04[i] - I24[i] * (gamma24[i] + rho[i]*ART_eligible_CD4_200_349 + rho_intervention[i] + alpha24[i] + mu[i] + nu) + rate_move_out[i] * I24[i] + sum(in_I24[i, ])
deriv(I25[]) = gamma24[i] * I24[i] + (RR_test_CD4200*tau[i] + tau_intervention[i]) * I05[i] - I25[i] * (rho[i]*ART_eligible_CD4_below_200 + rho_intervention[i] + alpha25[i] + mu[i] + nu) + rate_move_out[i] * I25[i] + sum(in_I25[i, ])

deriv(I32[]) = (rho[i]*ART_eligible_CD4_above_500 + rho_intervention[i]) * I22[i] + iota[i] * I42[i] - I32[i] * (gamma32[i] + phi2[i] + alpha32[i] + mu[i] + nu) + rate_move_out[i] * I32[i] + sum(in_I32[i, ])
deriv(I33[]) = gamma32[i] * I32[i] + (rho[i]*ART_eligible_CD4_350_500 + rho_intervention[i]) * I23[i] + iota[i] * I43[i] - I33[i] * (gamma33[i] + phi3[i] + alpha33[i] + mu[i] + nu) + rate_move_out[i] * I33[i] + sum(in_I33[i, ])
deriv(I34[]) = gamma33[i] * I33[i] + (rho[i]*ART_eligible_CD4_200_349 + rho_intervention[i]) * I24[i] + iota[i] * I44[i] - I34[i] * (gamma34[i] + phi4[i] + alpha34[i] + mu[i] + nu) + rate_move_out[i] * I34[i] + sum(in_I34[i, ])
deriv(I35[]) = gamma34[i] * I34[i] + (rho[i]*ART_eligible_CD4_below_200 + rho_intervention[i]) * I25[i] + iota[i] * I45[i] - I35[i] * (phi5[i] + alpha35[i] + mu[i] + nu) + rate_move_out[i] * I35[i] + sum(in_I35[i, ])

deriv(I42[]) = phi2[i] * I32[i] - I42[i] * (gamma42[i] + iota[i] + alpha42[i] + mu[i] + nu) + rate_move_out[i] * I42[i] + sum(in_I42[i, ])
deriv(I43[]) = gamma42[i] * I42[i] + phi3[i] * I33[i] - I43[i] * (gamma43[i] + iota[i] + alpha43[i] + mu[i] + nu) + rate_move_out[i] * I43[i] + sum(in_I43[i, ])
deriv(I44[]) = gamma43[i] * I43[i] + phi4[i] * I34[i] - I44[i] * (gamma44[i] + iota[i] + alpha44[i] + mu[i] + nu) + rate_move_out[i] * I44[i] + sum(in_I44[i, ])
deriv(I45[]) = gamma44[i] * I44[i] + phi5[i] * I35[i] - I45[i] * (iota[i] + alpha45[i] + mu[i] + nu) + rate_move_out[i] * I45[i] + sum(in_I45[i, ])

# output(ART_eligible_CD4_below_200) = ART_eligible_CD4_below_200

# sum of all compartments
N[] = S0[i] + S1a[i] + S1b[i] + S1c[i] + S1d[i] + I01[i] + I11[i] + I02[i] + I03[i] + I04[i] + I05[i] +
  I22[i] + I23[i] + I24[i] + I25[i] + I32[i] + I33[i] + I34[i] + I35[i] +
  I42[i] + I43[i] + I44[i] + I45[i]

Ntot = if (Ncat == 9) (N[1] + N[2] + N[3] + N[4] + N[5] + N[6] + N[7] + N[8]) else sum(N)
Ntot_inc_former_FSW_nonCot = sum(N)

output(Ntot_inc_former_FSW_nonCot) = Ntot_inc_former_FSW_nonCot

# births due to population growth
# epsilon = interpolate(epsilon_t, epsilon_y, "linear")
epsilon = interpolate(epsilon_t, epsilon_y, "constant")

new_people = epsilon * Ntot
output(new_people) = new_people

# new entrants into each group
new_people_in_group[] = (epsilon + mu[i] + nu) * Ntot * omega[i]


new_people_in_group_FSW_only[] =  rate_leave_pro_FSW * N[i] * fraction_FSW_foreign * FSW_ONLY[i]

rate_leave_pro_FSW = user()
fraction_FSW_foreign = user()
FSW_ONLY[] = user()
dim(FSW_ONLY) = Ncat

dim(new_people_in_group) = Ncat
output(new_people_in_group[]) = new_people_in_group

dim(new_people_in_group_FSW_only) = Ncat
output(new_people_in_group_FSW_only[]) = new_people_in_group_FSW_only

# MOVEMENT
##############################################################################


rate_move_in[,] = user()
rate_move_out[] = user()



#moving in
in_S0[,] <- if (i == j) 0 else rate_move_in[i, j] * S0[j]
in_S1a[,] <- if (i == j) 0 else rate_move_in[i, j] * S1a[j]
in_S1b[,] <- if (i == j) 0 else rate_move_in[i, j] * S1b[j]
in_S1c[,] <- if (i == j) 0 else rate_move_in[i, j] * S1c[j]
in_S1d[,] <- if (i == j) 0 else rate_move_in[i, j] * S1d[j]
in_I01[,] <- if (i == j) 0 else rate_move_in[i, j] * I01[j]
in_I11[,] <- if (i == j) 0 else rate_move_in[i, j] * I11[j]
in_I02[,] <- if (i == j) 0 else rate_move_in[i, j] * I02[j]
in_I03[,] <- if (i == j) 0 else rate_move_in[i, j] * I03[j]
in_I04[,] <- if (i == j) 0 else rate_move_in[i, j] * I04[j]
in_I05[,] <- if (i == j) 0 else rate_move_in[i, j] * I05[j]
in_I22[,] <- if (i == j) 0 else rate_move_in[i, j] * I22[j]
in_I23[,] <- if (i == j) 0 else rate_move_in[i, j] * I23[j]
in_I24[,] <- if (i == j) 0 else rate_move_in[i, j] * I24[j]
in_I25[,] <- if (i == j) 0 else rate_move_in[i, j] * I25[j]
in_I32[,] <- if (i == j) 0 else rate_move_in[i, j] * I32[j]
in_I33[,] <- if (i == j) 0 else rate_move_in[i, j] * I33[j]
in_I34[,] <- if (i == j) 0 else rate_move_in[i, j] * I34[j]
in_I35[,] <- if (i == j) 0 else rate_move_in[i, j] * I35[j]
in_I42[,] <- if (i == j) 0 else rate_move_in[i, j] * I42[j]
in_I43[,] <- if (i == j) 0 else rate_move_in[i, j] * I43[j]
in_I44[,] <- if (i == j) 0 else rate_move_in[i, j] * I44[j]
in_I45[,] <- if (i == j) 0 else rate_move_in[i, j] * I45[j]



##


############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################





# mortality due to HIV infection
alphaItot[] =
  alpha01[i] * I01[i] + alpha11[i] * I11[i] + alpha02[i] * I02[i] + alpha03[i] * I03[i] + alpha04[i] * I04[i] +
  alpha05[i] * I05[i] + alpha22[i] * I22[i] + alpha23[i] * I23[i] + alpha24[i] * I24[i] + alpha25[i] * I25[i] +
  alpha32[i] * I32[i] + alpha33[i] * I33[i] + alpha34[i] * I34[i] + alpha35[i] * I35[i] +
  alpha42[i] * I42[i] + alpha43[i] * I43[i] + alpha44[i] * I44[i] + alpha45[i] * I45[i]

# BALANCING OF PARTNERSHIPS
##############################################################################

c_t_comm[] = user()
dim(c_t_comm) = user()
c_t_noncomm[] = user()
dim(c_t_noncomm) = user()

c_y_comm[,] = user()
dim(c_y_comm) = user()
c_y_noncomm[,] = user()
dim(c_y_noncomm) = user()

c_comm[] = interpolate(c_t_comm, c_y_comm, "linear")
c_noncomm[] = interpolate(c_t_noncomm, c_y_noncomm, "linear")

c_comm_balanced[] <- c_comm[i]
c_noncomm_balanced[] <- c_noncomm[i]

# BALANCING COMMERCIAL PARTNERSHIPS
##############################################################################


##############
# BALANCING BY CHANGING THE NUMBER OF CLIENTS : DO I WANT TO DO THIS ?
# not the code below won't work: need to change the inits before passing parameters
# this needs to be only at t=0
# N[5] = if(Ncat == 9 && t == 1986) (c_comm[1] * N[1] + c_comm[2] * N[2]) / c_comm[5] else N[5]
# S0_init[5] = if(Ncat == 9 && t == 1986) (S0[5]/N[5])*(c_comm[1] * N[1] + c_comm[2] * N[2]) / c_comm[5] else S0_init[5]
# I0_init[5] = if(Ncat == 9 && t == 1986) (I0[5]/N[5])*(c_comm[1] * N[1] + c_comm[2] * N[2]) / c_comm[5] else I0_init[5]
##############

who_believe_comm = user()

# 1 = believe clients, 0 = believe FSW
c_comm_balanced[5] = if(Ncat == 9 && who_believe_comm == 0) (c_comm[1] * N[1] + c_comm[2] * N[2])/N[5] else c_comm_balanced[5]
c_comm_balanced[1] = if(Ncat == 9 && who_believe_comm == 1) (c_comm[5] * N[5] - c_comm[2] * N[2])/N[1] else c_comm_balanced[1]

output(who_believe_comm) = who_believe_comm

##############
# BALANCING BY CHANGING THE PARTNER CHANGE RATE OF CLIENTS
# c_comm_balanced[5] = if(Ncat == 9) (c_comm[1] * N[1] + c_comm[2] * N[2])/N[5] else c_comm_balanced[5]
##############

##############
# BALANCING BY CHANGING THE PARTNER CHANGE RATE OF PRO FSW
# c_comm_balanced[1] = if(Ncat == 9) (c_comm[5] * N[5] - c_comm[2] * N[2])/N[1] else c_comm_balanced[1]
##############


# BALANCING NON-COMMERCIAL PARTNERSHIPS
##############################################################################


##############
# BALANCING BY CHANGING THE PARTNER CHANGE RATE OF GPF (AND FORMER FSW)
c_noncomm_balanced[3] = if(Ncat == 9) (N[5] * c_noncomm[5] + N[6] * c_noncomm[6] - N[1] * c_noncomm[1] - N[2] * c_noncomm[2]) / (N[3] + N[4]) else c_noncomm_balanced[3]
c_noncomm_balanced[4] = if(Ncat == 9) c_noncomm_balanced[3] else c_noncomm_balanced[4]
##############

# CHECKING THE BALANCING IS CORRECT
##############################################################################
B_check_comm = if(Ncat == 9) c_comm_balanced[1]*N[1] + c_comm_balanced[2]*N[2] + c_comm_balanced[3]*N[3] + c_comm_balanced[4]*N[4] - c_comm_balanced[5]*N[5] - c_comm_balanced[6]*N[6] else 1
B_check_noncomm = if(Ncat == 9) c_noncomm_balanced[1]*N[1] + c_noncomm_balanced[2]*N[2] + c_noncomm_balanced[3]*N[3] + c_noncomm_balanced[4]*N[4] - c_noncomm_balanced[5]*N[5] - c_noncomm_balanced[6]*N[6] else 1



# PROBABILITY OF SEXUAL CONTACT (MIXING)
##############################################################################

# p_comm_balanced_old[,] = if(sum(M_comm[i,]) == 1 && M_comm[i,j] == 1) 1 else N[j] * c_comm_balanced[j] * M_comm[j, i] / (N[i] * c_comm_balanced[i])
# dim(p_comm_balanced_old) = c(Ncat, Ncat)
# output(p_comm_balanced_old[,]) = p_comm_balanced_old

comm_partnerships_requested[,] = N[j] * c_comm_balanced[j] * M_comm[i,j]
noncomm_partnerships_requested[,] = N[j] * c_noncomm_balanced[j] * M_noncomm[i,j]

dim(comm_partnerships_requested) = c(Ncat, Ncat)


dim(noncomm_partnerships_requested) = c(Ncat, Ncat)



p_comm[,] = if(M_comm[i, j] == 0) 0 else M_comm[i, j] * N[j] * c_comm_balanced[j] / sum(comm_partnerships_requested[i,])
p_noncomm[,] = if(M_noncomm[i, j] == 0) 0 else M_noncomm[i, j] * N[j] * c_noncomm_balanced[j] / sum(noncomm_partnerships_requested[i,])

# p_comm[,] = user()
# p_noncomm[,] = user()




# INTERPOLATING FUNCTIONS
##############################################################################


ART_eligible_CD4_above_500 = interpolate(ART_eligible_CD4_above_500_t, ART_eligible_CD4_above_500_y, "constant")
ART_eligible_CD4_350_500 = interpolate(ART_eligible_CD4_350_500_t, ART_eligible_CD4_350_500_y, "constant")
ART_eligible_CD4_200_349 = interpolate(ART_eligible_CD4_200_349_t, ART_eligible_CD4_200_349_y, "constant")
ART_eligible_CD4_below_200 = interpolate(ART_eligible_CD4_below_200_t, ART_eligible_CD4_below_200_y, "constant")


output(ART_eligible_CD4_above_500) = ART_eligible_CD4_above_500
output(ART_eligible_CD4_350_500) = ART_eligible_CD4_350_500
output(ART_eligible_CD4_200_349) = ART_eligible_CD4_200_349
output(ART_eligible_CD4_below_200) = ART_eligible_CD4_below_200


ART_eligible_CD4_above_500_t[] = user()
ART_eligible_CD4_350_500_t[] = user()
ART_eligible_CD4_200_349_t[] = user()
ART_eligible_CD4_below_200_t[] = user()
ART_eligible_CD4_above_500_y[] = user()
ART_eligible_CD4_350_500_y[] = user()
ART_eligible_CD4_200_349_y[] = user()
ART_eligible_CD4_below_200_y[] = user()

dim(ART_eligible_CD4_above_500_t) = user()
dim(ART_eligible_CD4_350_500_t) = user()
dim(ART_eligible_CD4_200_349_t) = user()
dim(ART_eligible_CD4_below_200_t) = user()
dim(ART_eligible_CD4_above_500_y) = user()
dim(ART_eligible_CD4_350_500_y) = user()
dim(ART_eligible_CD4_200_349_y) = user()
dim(ART_eligible_CD4_below_200_y) = user()


infect_ART_y[,] = user()
dim(infect_ART_y) = user()
output(infect_ART_y[,]) = infect_ART_y

infect_ART_t[] = user()
dim(infect_ART_t) = user()
output(infect_ART_t[]) = infect_ART_t

infect_ART[] = interpolate(infect_ART_t, infect_ART_y, "constant")



testing_prob_t[] = user()
testing_prob_y[,] = user()
dim(testing_prob_t) = user()
dim(testing_prob_y) = user()
dim(testing_prob) = Ncat
testing_prob[] = interpolate(testing_prob_t, testing_prob_y, "linear")
output(testing_prob[]) = testing_prob
#
# ART_prob_t[] = user()
# ART_prob_y[,] = user()
# dim(ART_prob_t) = user()
# dim(ART_prob_y) = user()
# dim(ART_prob) = Ncat

# ART_prob[] = interpolate(ART_prob_t, ART_prob_y, "constant")

# output(ART_prob[]) = ART_prob


prep_intervention_t[] = user()
prep_intervention_y[,] = user()

fc_comm[,] = interpolate(fc_t_comm, fc_y_comm, "linear")
fP_comm[] = interpolate(fP_t_comm, fP_y_comm, "linear")

fc_noncomm[,] = interpolate(fc_t_noncomm, fc_y_noncomm, "linear")
fP_noncomm[] = interpolate(fP_t_noncomm, fP_y_noncomm, "linear")

n_comm[,] = interpolate(n_t_comm, n_y_comm, "linear")
n_noncomm[,] = interpolate(n_t_noncomm, n_y_noncomm, "linear")

# prep_intervention is offering rate sigma is acceptance rate
zeta[] = (tau[i] + tau_intervention[i]) * sigma[i] * prep_offered[i] * PrEPOnOff

tau_intervention_t[] = user()
tau_intervention_y[,] = user()

tau_intervention[] = interpolate(tau_intervention_t, tau_intervention_y, "constant")

dim(tau_intervention) = Ncat

dim(tau_intervention_t) = user()
dim(tau_intervention_y) = user()

rho_intervention_t[] = user()
rho_intervention_y[,] = user()

rho_intervention[] = interpolate(rho_intervention_t, rho_intervention_y, "constant")
dim(rho_intervention) = Ncat

dim(rho_intervention_t) = user()
dim(rho_intervention_y) = user()



output(tau_intervention[]) = tau_intervention
output(rho_intervention[]) = rho_intervention


PrEPOnOff = user()
output(PrEPOnOff) = PrEPOnOff
fPa = user()
fPb = user()
fPc = user()

output(fPa) = fPa
output(fPb) = fPb
output(fPc) = fPc

dim(zeta) = Ncat

output(zeta[]) = zeta

dim(sigma) = Ncat
sigma[] = user()


prep_offered[] = interpolate(prep_intervention_t, prep_intervention_y, "constant")

#FOI of j on i
lambda[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                              I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                              I42[j], I43[j], I44[j], I45[j],
                                              N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP[i], ec[i],
                                              fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])

#FOI of j on i. PrEP adherence category 0 (off PrEP)
lambda_0[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                                I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                                I42[j], I43[j], I44[j], I45[j],
                                                N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP0[i], ec[i],
                                                fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])
#FOI of j on i. PrEP adherence category 1a (daily adherence)
lambda_1a[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                                 I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                                 I42[j], I43[j], I44[j], I45[j],
                                                 N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP1a[i], ec[i],
                                                 fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])
#FOI of j on i. PrEP adherence category 1b (intermittent adherence)
lambda_1b[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                                 I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                                 I42[j], I43[j], I44[j], I45[j],
                                                 N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP1b[i], ec[i],
                                                 fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])
#FOI of j on i. PrEP adherence category 1c (no adherence)
lambda_1c[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                                 I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                                 I42[j], I43[j], I44[j], I45[j],
                                                 N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP1c[i], ec[i],
                                                 fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])
#FOI of j on i. PrEP adherence category 1d (dropout)
lambda_1d[,] = if (i == j) 0 else compute_lambda(c_comm_balanced[i], p_comm[i,j], S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
                                                 I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
                                                 I42[j], I43[j], I44[j], I45[j],
                                                 N[j], beta_comm[i], R, fc_comm[i,j], fP_comm[i], n_comm[i,j], eP1d[i], ec[i],
                                                 fc_noncomm[i,j], fP_noncomm[i], n_noncomm[i,j], c_noncomm_balanced[i], p_noncomm[i,j], infect_ART[j], infect_acute, infect_AIDS, beta_noncomm[i])

lambda_sum_0[] = sum(lambda_0[i,])
lambda_sum_1a[] = sum(lambda_1a[i,])
lambda_sum_1b[] = sum(lambda_1b[i,])
lambda_sum_1c[] = sum(lambda_1c[i,])
lambda_sum_1d[] = sum(lambda_1d[i,])

# TESTING
##############################################################################
dim(tau) = Ncat
output(tau[]) = tau

RR_test_CD4200 = user()

tau[] = -log(1-testing_prob[i])


test_rate_prep[] = user()
dim(test_rate_prep) = Ncat


# ART
##############################################################################
dim(rho) = Ncat
output(rho[]) = rho

dim(iota) = Ncat
output(iota[]) = iota

iota[] = user()


# rho[] = -log(1-ART_prob[i])

rho[] = user()

output(infect_ART[]) = infect_ART




# OUTPUTS
##############################################################################


FSW_out = rate_move_out[1] * N[1]



output(FSW_out) = FSW_out


FSW_in = rate_move_in[1, 3] * N[3]

output(FSW_in) = FSW_in


output(rate_move_in[,]) = rate_move_in
output(rate_move_out[]) = rate_move_out


output(E0[]) = E0
output(E1a[]) = E1a
output(E1b[]) = E1b
output(E1c[]) = E1c
output(E1d[]) = E1d

deriv(cumuInf[]) = S0[i] * lambda_sum_0[i] + S1a[i] * lambda_sum_1a[i] + S1b[i] * lambda_sum_1b[i] + S1c[i] * lambda_sum_1c[i] + S1d[i] * lambda_sum_1d[i]
deriv(OnPrEP[]) = zeta[i] * S0[i]

deriv(OnPrEP1a[])  = zeta[i] * fPa * S0[i]
deriv(OnPrEP1b[])  = zeta[i] * fPb * S0[i]
deriv(OnPrEP1c[])  = zeta[i] * fPc * S0[i]


deriv(cumuHIVDeaths[]) = alpha01[i] * I01[i] + alpha11[i] * I11[i] + alpha02[i] * I02[i] + alpha03[i] * I03[i] + alpha04[i] * I04[i] +
  alpha05[i] * I05[i] + alpha22[i] * I22[i] + alpha23[i] * I23[i] + alpha24[i] * I24[i] + alpha25[i] * I25[i] +
  alpha32[i] * I32[i] + alpha33[i] * I33[i] + alpha34[i] * I34[i] + alpha35[i] * I35[i] +
  alpha42[i] * I42[i] + alpha43[i] * I43[i] + alpha44[i] * I44[i] + alpha45[i] * I45[i]

deriv(cumuAllDeaths[]) = (alpha01[i] + mu[i]) * I01[i] + (alpha11[i] + mu[i]) * I11[i] + (alpha02[i] + mu[i]) * I02[i] + (alpha03[i] + mu[i]) * I03[i] + (alpha04[i] + mu[i]) * I04[i] +
  (alpha05[i] + mu[i]) * I05[i] + (alpha22[i] + mu[i]) * I22[i] + (alpha23[i] + mu[i]) * I23[i] + (alpha24[i] + mu[i]) * I24[i] + (alpha25[i] + mu[i]) * I25[i] +
  (alpha32[i] + mu[i]) * I32[i] + (alpha33[i] + mu[i]) * I33[i] + (alpha34[i] + mu[i]) * I34[i] + (alpha35[i] + mu[i]) * I35[i] +
  (alpha42[i] + mu[i]) * I42[i] + (alpha43[i] + mu[i]) * I43[i] + (alpha44[i] + mu[i]) * I44[i] + (alpha45[i] + mu[i]) * I45[i] +
  mu[i] * (S0[i] + S1a[i] + S1b[i] + S1c[i] + S1d[i])

deriv(cumuARTinitiations[]) = (rho_intervention[i] + rho[i]*ART_eligible_CD4_above_500) * I22[i] + (rho[i]*ART_eligible_CD4_350_500 + rho_intervention[i]) * I23[i] +
  (rho[i]*ART_eligible_CD4_200_349 + rho_intervention[i]) * I24[i] + (rho[i]*ART_eligible_CD4_below_200 + rho_intervention[i]) * I25[i]



deriv(cumuARTREinitiations[]) = iota[i] * I42[i] + iota[i] * I43[i] + iota[i] * I44[i] + iota[i] * I45[i]




deriv(S0_infections[]) = S0[i] * lambda_sum_0[i]
deriv(S1a_infections[]) = S1a[i] * lambda_sum_1a[i]
deriv(S1b_infections[]) = S1b[i] * lambda_sum_1b[i]
deriv(S1c_infections[]) = S1c[i] * lambda_sum_1c[i]
deriv(S1d_infections[]) = S1d[i] * lambda_sum_1d[i]

initial(S0_infections[]) = 0
initial(S1a_infections[]) = 0
initial(S1b_infections[]) = 0
initial(S1c_infections[]) = 0
initial(S1d_infections[]) = 0

dim(S0_infections) = Ncat
dim(S1a_infections) = Ncat
dim(S1b_infections) = Ncat
dim(S1c_infections) = Ncat
dim(S1d_infections) = Ncat

# deriv(cumuTesting[]) = (tau[i] + tau_intervention[i] * PrEPOnOff) * S0[i] + (tau[i] + tau_intervention[i] * PrEPOnOff) * I01[i] +
#   test_rate_prep[i] * I11[i] + (tau[i] + tau_intervention[i] * PrEPOnOff) * I02[i] + (tau[i] + tau_intervention[i] * PrEPOnOff) * I03[i] +
#   (tau[i] + tau_intervention[i] * PrEPOnOff) * I04[i] + (RR_test_CD4200 * tau[i] + tau_intervention[i] * PrEPOnOff) * I05[i]



cumuInftot = sum(cumuInf)
output(cumuInftot) = cumuInftot

# fraction of group in each category INCLUDING FORMER FSW OUTSIDE BENIN
frac_N[] = N[i] / Ntot#_inc_former_FSW_nonCot
frac_F[] = if(Ncat == 9) (N[1] + N[2] + N[3] + N[4] + N[7])/ Ntot else 0
frac_N_sexualpop[] = if(Ncat == 9) N[i] / (N[1] + N[2] + N[3] + N[4] + N[5] + N[6]) else 0

frac_virgin = if(Ncat == 9) (N[7] + N[8])/(N[1] + N[2] + N[3] + N[4] + N[5] + N[6] + N[7] + N[8]) else 0
output(frac_virgin) = frac_virgin

dim(frac_F) = Ncat
output(frac_F[]) = frac_F
output(frac_N_sexualpop[]) = frac_N_sexualpop
dim (frac_N_sexualpop) = Ncat
# Calculations

output(alphaItot[]) = alphaItot


# PREVALENCE
# n.b. prevalence for all ages in each risk group will be useful, so maybe one prev array is too much info in one output
prev_FSW = 100 * (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
                    I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
                    I42[1] + I43[1] + I44[1] + I45[1]) / N[1]

prev_LowFSW = 100 * (I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
                       I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
                       I42[2] + I43[2] + I44[2] + I45[2]) / N[2]

output(prev_LowFSW) = prev_LowFSW

prev_client = 100 * (I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
                       I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
                       I42[5] + I43[5] + I44[5] + I45[5]) / N[5]

prev[] = 100 * (I01[i] + I11[i] + I02[i] + I03[i] + I04[i] + I05[i] +
                  I22[i] + I23[i] + I24[i] + I25[i] + I32[i] + I33[i] + I34[i] + I35[i] +
                  I42[i] + I43[i] + I44[i] + I45[i]) / N[i]


# prev_men = (prev[5] * N[5] + prev[6] * N[6] + prev[8] * N[8])/(N[5] + N[6] + N[8])
# prev_men_check = 100 * (I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
#                     I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
#                     I42[5] + I43[5] + I44[5] + I45[5] +
#                     I01[6] + I11[6] + I02[6] + I03[6] + I04[6] + I05[6] +
#                     I22[6] + I23[6] + I24[6] + I25[6] + I32[6] + I33[6] + I34[6] + I35[6] +
#                     I42[6] + I43[6] + I44[6] + I45[6] +
#                     I01[8] + I11[8] + I02[8] + I03[8] + I04[8] + I05[8] +
#                     I22[8] + I23[8] + I24[8] + I25[8] + I32[8] + I33[8] + I34[8] + I35[8] +
#                     I42[8] + I43[8] + I44[8] + I45[8]
#                     ) / (N[5] + N[6] + N[8])
prev_men = 100 * (I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
                    I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
                    I42[5] + I43[5] + I44[5] + I45[5] +
                    I01[6] + I11[6] + I02[6] + I03[6] + I04[6] + I05[6] +
                    I22[6] + I23[6] + I24[6] + I25[6] + I32[6] + I33[6] + I34[6] + I35[6] +
                    I42[6] + I43[6] + I44[6] + I45[6] +
                    I01[8] + I11[8] + I02[8] + I03[8] + I04[8] + I05[8] +
                    I22[8] + I23[8] + I24[8] + I25[8] + I32[8] + I33[8] + I34[8] + I35[8] +
                    I42[8] + I43[8] + I44[8] + I45[8]) / (N[5] + N[6] + N[8])

# prev_women = (prev[1] * N[1] + prev[2] * N[2] + prev[3] * N[3] + prev[4] * N[4] + prev[7] * N[7]) /
#   (N[1] + N[2] + N[3] + N[4] + N[7])
prev_women = 100 * (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
                      I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
                      I42[1] + I43[1] + I44[1] + I45[1] +
                      I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
                      I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
                      I42[2] + I43[2] + I44[2] + I45[2] +
                      I01[3] + I11[3] + I02[3] + I03[3] + I04[3] + I05[3] +
                      I22[3] + I23[3] + I24[3] + I25[3] + I32[3] + I33[3] + I34[3] + I35[3] +
                      I42[3] + I43[3] + I44[3] + I45[3] +
                      I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
                      I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
                      I42[4] + I43[4] + I44[4] + I45[4] +
                      I01[7] + I11[7] + I02[7] + I03[7] + I04[7] + I05[7] +
                      I22[7] + I23[7] + I24[7] + I25[7] + I32[7] + I33[7] + I34[7] + I35[7] +
                      I42[7] + I43[7] + I44[7] + I45[7]) / (N[1] + N[2] + N[3] + N[4] + N[7])

prev_cotonou = 100 * (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
                        I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
                        I42[1] + I43[1] + I44[1] + I45[1] +
                        I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
                        I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
                        I42[2] + I43[2] + I44[2] + I45[2] +
                        I01[3] + I11[3] + I02[3] + I03[3] + I04[3] + I05[3] +
                        I22[3] + I23[3] + I24[3] + I25[3] + I32[3] + I33[3] + I34[3] + I35[3] +
                        I42[3] + I43[3] + I44[3] + I45[3] +
                        I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
                        I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
                        I42[4] + I43[4] + I44[4] + I45[4] +
                        I01[7] + I11[7] + I02[7] + I03[7] + I04[7] + I05[7] +
                        I22[7] + I23[7] + I24[7] + I25[7] + I32[7] + I33[7] + I34[7] + I35[7] +
                        I42[7] + I43[7] + I44[7] + I45[7] +
                        I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
                        I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
                        I42[5] + I43[5] + I44[5] + I45[5] +
                        I01[6] + I11[6] + I02[6] + I03[6] + I04[6] + I05[6] +
                        I22[6] + I23[6] + I24[6] + I25[6] + I32[6] + I33[6] + I34[6] + I35[6] +
                        I42[6] + I43[6] + I44[6] + I45[6] +
                        I01[8] + I11[8] + I02[8] + I03[8] + I04[8] + I05[8] +
                        I22[8] + I23[8] + I24[8] + I25[8] + I32[8] + I33[8] + I34[8] + I35[8] +
                        I42[8] + I43[8] + I44[8] + I45[8]) / (N[1] + N[2] + N[3] + N[4] + N[7] + N[5] + N[6] + N[8])


# prev_women_check = 100 * (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
#                             I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
#                             I42[1] + I43[1] + I44[1] + I45[1] +
#                             I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
#                             I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
#                             I42[2] + I43[2] + I44[2] + I45[2] +
#                           I01[3] + I11[3] + I02[3] + I03[3] + I04[3] + I05[3] +
#                             I22[3] + I23[3] + I24[3] + I25[3] + I32[3] + I33[3] + I34[3] + I35[3] +
#                             I42[3] + I43[3] + I44[3] + I45[3] +
#                           I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
#                             I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
#                             I42[4] + I43[4] + I44[4] + I45[4] +
#                           I01[7] + I11[7] + I02[7] + I03[7] + I04[7] + I05[7] +
#                             I22[7] + I23[7] + I24[7] + I25[7] + I32[7] + I33[7] + I34[7] + I35[7] +
#                             I42[7] + I43[7] + I44[7] + I45[7]) / (N[1] + N[2] + N[3] + N[4] + N[7])

output(prev_men) = prev_men
# output(prev_men_check) = prev_men_check
output(prev_women) = prev_women
# output(prev_women_check) = prev_women_check

output(prev_cotonou) = prev_cotonou
# prev_test1 = 100*(I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
#   I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
#   I42[5] + I43[5] + I44[5] + I45[5] +
#   I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
#   I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
#   I42[4] + I43[4] + I44[4] + I45[4]) / (N[4] + N[5])
#
# prev_test2 = (prev[5] * N[5] + prev[4] * N[4])/(N[4] + N[5])
#
# output(prev_test1) = prev_test1
# output(prev_test2) = prev_test2

ART_coverage_men = (I32[5] + I33[5] + I34[5] + I35[5] +
                      I32[6] + I33[6] + I34[6] + I35[6]) /
    (I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
     I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
     I42[5] + I43[5] + I44[5] + I45[5] +
     I01[6] + I11[6] + I02[6] + I03[6] + I04[6] + I05[6] +
     I22[6] + I23[6] + I24[6] + I25[6] + I32[6] + I33[6] + I34[6] + I35[6] +
     I42[6] + I43[6] + I44[6] + I45[6])

ART_coverage_women = (I32[1] + I33[1] + I34[1] + I35[1] +
                        I32[2] + I33[2] + I34[2] + I35[2] +
                        I32[3] + I33[3] + I34[3] + I35[3] +
                        I32[4] + I33[4] + I34[4] + I35[4]) /
  (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
     I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
     I42[1] + I43[1] + I44[1] + I45[1] +
     I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
     I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
     I42[2] + I43[2] + I44[2] + I45[2] +
     I01[3] + I11[3] + I02[3] + I03[3] + I04[3] + I05[3] +
     I22[3] + I23[3] + I24[3] + I25[3] + I32[3] + I33[3] + I34[3] + I35[3] +
     I42[3] + I43[3] + I44[3] + I45[3] +
     I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
     I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
     I42[4] + I43[4] + I44[4] + I45[4])

ART_coverage_FSW = (I32[1] + I33[1] + I34[1] + I35[1]) / (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
     I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
     I42[1] + I43[1] + I44[1] + I45[1])

ART_coverage_all = (I32[1] + I33[1] + I34[1] + I35[1] +
                      I32[2] + I33[2] + I34[2] + I35[2] +
                      I32[3] + I33[3] + I34[3] + I35[3] +
                      I32[4] + I33[4] + I34[4] + I35[4] +
                      I32[5] + I33[5] + I34[5] + I35[5] +
                      I32[6] + I33[6] + I34[6] + I35[6]) /
  (I01[1] + I11[1] + I02[1] + I03[1] + I04[1] + I05[1] +
     I22[1] + I23[1] + I24[1] + I25[1] + I32[1] + I33[1] + I34[1] + I35[1] +
     I42[1] + I43[1] + I44[1] + I45[1] +
     I01[2] + I11[2] + I02[2] + I03[2] + I04[2] + I05[2] +
     I22[2] + I23[2] + I24[2] + I25[2] + I32[2] + I33[2] + I34[2] + I35[2] +
     I42[2] + I43[2] + I44[2] + I45[2] +
     I01[3] + I11[3] + I02[3] + I03[3] + I04[3] + I05[3] +
     I22[3] + I23[3] + I24[3] + I25[3] + I32[3] + I33[3] + I34[3] + I35[3] +
     I42[3] + I43[3] + I44[3] + I45[3] +
     I01[4] + I11[4] + I02[4] + I03[4] + I04[4] + I05[4] +
     I22[4] + I23[4] + I24[4] + I25[4] + I32[4] + I33[4] + I34[4] + I35[4] +
     I42[4] + I43[4] + I44[4] + I45[4] +
     I01[5] + I11[5] + I02[5] + I03[5] + I04[5] + I05[5] +
     I22[5] + I23[5] + I24[5] + I25[5] + I32[5] + I33[5] + I34[5] + I35[5] +
     I42[5] + I43[5] + I44[5] + I45[5] +
     I01[6] + I11[6] + I02[6] + I03[6] + I04[6] + I05[6] +
     I22[6] + I23[6] + I24[6] + I25[6] + I32[6] + I33[6] + I34[6] + I35[6] +
     I42[6] + I43[6] + I44[6] + I45[6])


ART_coverage[] = (I32[i] + I33[i] + I34[i] + I35[i]) /
  (I01[i] + I11[i] + I02[i] + I03[i] + I04[i] + I05[i] +
     I22[i] + I23[i] + I24[i] + I25[i] + I32[i] + I33[i] + I34[i] + I35[i] +
     I42[i] + I43[i] + I44[i] + I45[i])

output(ART_coverage[]) = ART_coverage

dim(ART_coverage) = Ncat

Men_on_ART = (I32[5] + I33[5] + I34[5] + I35[5] +
                I32[6] + I33[6] + I34[6] + I35[6])
Women_on_ART = (I32[1] + I33[1] + I34[1] + I35[1] +
                  I32[2] + I33[2] + I34[2] + I35[2] +
                  I32[3] + I33[3] + I34[3] + I35[3] +
                  I32[4] + I33[4] + I34[4] + I35[4])

ART_sex_ratio = Women_on_ART/Men_on_ART
output(ART_sex_ratio) = ART_sex_ratio

pc_S1a = (S1a[1]) * 100 / (S1a[1] + S1b[1] + S1c[1])
pc_S1b = (S1b[1]) * 100 / (S1a[1] + S1b[1] + S1c[1])
pc_S1c = (S1c[1]) * 100 / (S1a[1] + S1b[1] + S1c[1])

output(pc_S1a) = pc_S1a
output(pc_S1b) = pc_S1b
output(pc_S1c) = pc_S1c


output(Men_on_ART) = Men_on_ART
output(Women_on_ART) = Women_on_ART

output(ART_coverage_men) = ART_coverage_men
output(ART_coverage_women) = ART_coverage_women
output(ART_coverage_FSW) = ART_coverage_FSW
output(ART_coverage_all) = ART_coverage_all

output(Ncat) = Ncat
output(omega[]) = omega
output(Ntot) = Ntot
output(B_check_comm) = B_check_comm
output(B_check_noncomm) = B_check_noncomm

output(N[]) = N # is it worth outputting N? Once we have ages, it'll be better to have separate Ns for risk groups... but eugene ages can make N a matrix!
output(prev_FSW) = prev_FSW
output(prev_client) = prev_client
output(prev[]) = prev
output(frac_N[]) = frac_N
# output(lambda_sum[]) = lambda_sum

# FOI outputs
output(lambda[,]) = lambda
output(lambda_0[,]) = lambda_0
output(lambda_1a[,]) = lambda_1a
output(lambda_1b[,]) = lambda_1b
output(lambda_1c[,]) = lambda_1c
output(lambda_1d[,]) = lambda_1d

output(lambda_sum_0[]) = lambda_sum_0
output(lambda_sum_1a[]) = lambda_sum_1a
output(lambda_sum_1b[]) = lambda_sum_1b
output(lambda_sum_1c[]) = lambda_sum_1c
output(lambda_sum_1d[]) = lambda_sum_1d

output(fc_comm[,]) = fc_comm
output(fP_comm[]) = fP_comm
output(fc_noncomm[,]) = fc_noncomm
output(fP_noncomm[]) = fP_noncomm
output(c_comm[]) = c_comm
output(c_comm_balanced[]) = c_comm_balanced
output(c_noncomm[]) = c_noncomm
output(c_noncomm_balanced[]) = c_noncomm_balanced
# output(B[,]) = B
output(p_comm[,]) = p_comm
output(p_noncomm[,]) = p_noncomm

output(n_comm[,]) = n_comm
output(n_noncomm[,]) = n_noncomm

output(theta[,]) = theta

# output(rate_move_in[,]) = rate_move_in
# output(rate_move_out[]) = rate_move_out


output(epsilon) = epsilon
output(M_comm[,]) = M_comm
output(M_noncomm[,]) = M_noncomm

output(beta_comm[]) = beta_comm
output(beta_noncomm[]) = beta_noncomm

output(alpha05[]) = alpha05
output(alpha35[]) = alpha35

output(alpha02[]) = alpha02
output(alpha22[]) = alpha22
output(alpha01[]) = alpha01
output(alpha11[]) = alpha11

# INCIDENCE RATE

# = no. disease onsets / sum of "person-time" at risk
output(nu) = nu

output(mu[]) = mu
output(gamma01[]) = gamma01
output(gamma11[]) = gamma11

output(kappaa[]) = kappaa
output(kappab[]) = kappab
output(kappac[]) = kappac


output(kappa1[]) = kappa1

output(Nage) = Nage
output(dur_FSW) = dur_FSW
# output(rate_leave_FSW) = rate_leave_FSW
# output(rate_move_GPF_pFSW) = rate_move_GPF_pFSW
# output(rate_leave_client) = rate_leave_client
# output(prev[]) = prev

# in future nb eP eC constants
# lambda[,] = compute_lambda(S0[j], S1a[j], S1b[j], S1c[j], I01[j], I11[j], I02[j], I03[j], I04[j], I05[j],
#                              I22[j], I23[j], I24[j], I25[j], I32[j], I33[j], I34[j], I35[j],
#                              I42[j], I43[j], I44[j], I45[j],
#                              N[j], beta[i,j], R[i,j], fc[i,j], fP[i,j], n[i,j], eP[j], ec[j])

# parameters

Ncat = user()
Nage = user()


# do this for all vars!
initial(S0[]) = S0_init[i]
S0_init[] = user()
initial(S1a[]) = S1a_init[i]
S1a_init[] = user()
initial(S1b[]) = S1b_init[i]
S1b_init[] = user()
initial(S1c[]) = S1c_init[i]
S1c_init[] = user()
initial(S1d[]) = S1d_init[i]
S1d_init[] = user()

initial(I01[]) = I01_init[i]
I01_init[] = user()
initial(I11[]) = I11_init[i]
I11_init[] = user()

initial(I02[]) = I02_init[i]
I02_init[] = user()
initial(I03[]) = I03_init[i]
I03_init[] = user()
initial(I04[]) = I04_init[i]
I04_init[] = user()
initial(I05[]) = I05_init[i]
I05_init[] = user()

initial(I22[]) = I22_init[i]
I22_init[] = user()
initial(I23[]) = I23_init[i]
I23_init[] = user()
initial(I24[]) = I24_init[i]
I24_init[] = user()
initial(I25[]) = I25_init[i]
I25_init[] = user()

initial(I32[]) = I32_init[i]
I32_init[] = user()
initial(I33[]) = I33_init[i]
I33_init[] = user()
initial(I34[]) = I34_init[i]
I34_init[] = user()
initial(I35[]) = I35_init[i]
I35_init[] = user()

initial(I42[]) = I42_init[i]
I42_init[] = user()
initial(I43[]) = I43_init[i]
I43_init[] = user()
initial(I44[]) = I44_init[i]
I44_init[] = user()
initial(I45[]) = I45_init[i]
I45_init[] = user()

initial(cumuInf[]) = cumuInf_init[i]
cumuInf_init[] = user()

initial(cumuHIVDeaths[]) = 0
initial(cumuARTREinitiations[]) = 0
# initial(cumuTesting[]) = 0

initial(cumuARTinitiations[]) = 0

initial(cumuAllDeaths[]) = 0

initial(OnPrEP[]) = OnPrEP_init[i]
OnPrEP_init[] = user()

initial(OnPrEP1a[]) = 0
initial(OnPrEP1b[]) = 0
initial(OnPrEP1c[]) = 0




# initial(S0) = user()
# initial(S1a) = user()
# initial(S1b) = user()
# initial(S1c) = user()
#
# initial(I01) = user()
# initial(I11) = user()
#
# initial(I02) = user()
# initial(I03) = user()
# initial(I04) = user()
# initial(I05) = user()
#
# initial(I22) = user()
# initial(I23) = user()
# initial(I24) = user()
# initial(I25) = user()
#
# initial(I32) = user()
# initial(I33) = user()
# initial(I34) = user()
# initial(I35) = user()
#
# initial(I42) = user()
# initial(I43) = user()
# initial(I44) = user()
# initial(I45) = user()

dim(infect_ART) = Ncat

infect_acute = user()
infect_AIDS = user()

nu = user()

mu[] = user()
gamma01[] = user()
gamma02[] = user()
gamma03[] = user()
gamma04[] = user()

gamma11[] = user()

gamma22[] = user()
gamma23[] = user()
gamma24[] = user()

gamma32[] = user()
gamma33[] = user()
gamma34[] = user()

gamma42[] = user()
gamma43[] = user()
gamma44[] = user()


phi2[] = user()
phi3[] = user()
phi4[] = user()
phi5[] = user()

psia[] = user()
psib[] = user()




kappa1[] = user()
kappaa[] = user()
kappab[] = user()
kappac[] = user()

# note alpha is ordered differently...
alpha01[] = user()
alpha02[] = user()
alpha03[] = user()
alpha04[] = user()
alpha05[] = user()

alpha11[] = user()

alpha22[] = user()
alpha23[] = user()
alpha24[] = user()
alpha25[] = user()

alpha32[] = user()
alpha33[] = user()
alpha34[] = user()
alpha35[] = user()

alpha42[] = user()
alpha43[] = user()
alpha44[] = user()
alpha45[] = user()

# FOI parameters

beta_comm[] = user()
beta_noncomm[] = user()

# c_comm[] = user()
# c_noncomm[] = user()
# p_comm[,] = user()
# p_noncomm[,] = user()

ec[] = user()
eP[] = user()
eP0[] = user()
eP1a[] = user()
eP1b[] = user()
eP1c[] = user()
eP1d[] = user()



M_comm[,] = user()
M_noncomm[,] = user()

R = user()

# growth

omega[] = user()

# balancing
theta[,] = user()

epsilon_t[] = user()
epsilon_y[] = user()

dim(epsilon_t) = user()
dim(epsilon_y) = user()

fc_t_comm[] = user()
fc_y_comm[,,] = user()
dim(fc_t_comm) = user()
dim(fc_y_comm) = user()

fP_t_comm[] = user()
fP_y_comm[,] = user()
dim(fP_t_comm) = user()
dim(fP_y_comm) = user()

fc_t_noncomm[] = user()
fc_y_noncomm[,,] = user()
dim(fc_t_noncomm) = user()
dim(fc_y_noncomm) = user()

fP_t_noncomm[] = user()
fP_y_noncomm[,] = user()
dim(fP_t_noncomm) = user()
dim(fP_y_noncomm) = user()



n_t_comm[] = user()
n_y_comm[,,] = user()
dim(n_t_comm) = user()
dim(n_y_comm) = user()
n_t_noncomm[] = user()
n_y_noncomm[,,] = user()
dim(n_t_noncomm) = user()
dim(n_y_noncomm) = user()


dur_FSW = user()



# DIMMING



dim(prep_offered) = Ncat

#parameters

dim(omega) = Ncat
dim(frac_N) = Ncat

# care cascade
dim(mu) = Ncat

dim(gamma01) = Ncat
dim(gamma02) = Ncat
dim(gamma03) = Ncat
dim(gamma04) = Ncat
dim(gamma11) = Ncat
dim(gamma22) = Ncat
dim(gamma23) = Ncat
dim(gamma24) = Ncat
dim(gamma32) = Ncat
dim(gamma33) = Ncat
dim(gamma34) = Ncat
dim(gamma42) = Ncat
dim(gamma43) = Ncat
dim(gamma44) = Ncat


dim(phi2) = Ncat
dim(phi3) = Ncat
dim(phi4) = Ncat
dim(phi5) = Ncat

dim(psia) = Ncat
dim(psib) = Ncat


dim(prep_intervention_t) = user()
dim(prep_intervention_y) = user()

dim(kappaa) = Ncat
dim(kappab) = Ncat
dim(kappac) = Ncat
dim(kappa1) = Ncat

dim(alpha01) = Ncat
dim(alpha02) = Ncat
dim(alpha03) = Ncat
dim(alpha04) = Ncat
dim(alpha05) = Ncat
dim(alpha11) = Ncat
dim(alpha22) = Ncat
dim(alpha23) = Ncat
dim(alpha24) = Ncat
dim(alpha25) = Ncat
dim(alpha32) = Ncat
dim(alpha33) = Ncat
dim(alpha34) = Ncat
dim(alpha35) = Ncat
dim(alpha42) = Ncat
dim(alpha43) = Ncat
dim(alpha44) = Ncat
dim(alpha45) = Ncat

# FOI parameters
# dim(beta) = Ncat
dim(beta_comm) = Ncat
dim(beta_noncomm) = Ncat

dim(c_comm) = user()
dim(c_noncomm) = user()
dim(p_comm) = c(Ncat, Ncat)
dim(p_noncomm) = c(Ncat, Ncat)

dim(ec) = Ncat
dim(eP) = Ncat
dim(eP0) = Ncat
dim(eP1a) = Ncat
dim(eP1b) = Ncat
dim(eP1c) = Ncat
dim(eP1d) = Ncat

# dim(epsilon) = Ncat
dim(fc_comm) = c(Ncat, Ncat)
dim(fP_comm) = Ncat
dim(fc_noncomm) = c(Ncat, Ncat)
dim(fP_noncomm) = Ncat

dim(n_comm) = c(Ncat, Ncat)
dim(n_noncomm) = c(Ncat, Ncat)


dim(cumuHIVDeaths) = Ncat
dim(cumuARTinitiations) = Ncat
dim(cumuARTREinitiations) = Ncat
# dim(cumuTesting) = Ncat


dim(cumuAllDeaths) = Ncat


dim(cumuInf) = Ncat
dim(OnPrEP) = Ncat

dim(OnPrEP1a) = Ncat
dim(OnPrEP1b) = Ncat
dim(OnPrEP1c) = Ncat


# states and initial conditions
dim(S0) = Ncat
dim(S1a) = Ncat
dim(S1b) = Ncat
dim(S1c) = Ncat
dim(S1d) = Ncat

dim(I01) = Ncat
dim(I11) = Ncat

dim(I02) = Ncat
dim(I03) = Ncat
dim(I04) = Ncat
dim(I05) = Ncat

dim(I22) = Ncat
dim(I23) = Ncat
dim(I24) = Ncat
dim(I25) = Ncat

dim(I32) = Ncat
dim(I33) = Ncat
dim(I34) = Ncat
dim(I35) = Ncat

dim(I42) = Ncat
dim(I43) = Ncat
dim(I44) = Ncat
dim(I45) = Ncat

dim(S0_init) = Ncat
dim(S1a_init) = Ncat
dim(S1b_init) = Ncat
dim(S1c_init) = Ncat
dim(S1d_init) = Ncat

dim(I01_init) = Ncat
dim(I11_init) = Ncat

dim(I02_init) = Ncat
dim(I03_init) = Ncat
dim(I04_init) = Ncat
dim(I05_init) = Ncat

dim(I22_init) = Ncat
dim(I23_init) = Ncat
dim(I24_init) = Ncat
dim(I25_init) = Ncat

dim(I32_init) = Ncat
dim(I33_init) = Ncat
dim(I34_init) = Ncat
dim(I35_init) = Ncat

dim(I42_init) = Ncat
dim(I43_init) = Ncat
dim(I44_init) = Ncat
dim(I45_init) = Ncat




# other variables
dim(N) = Ncat
dim(E0) = Ncat
dim(E1a) = Ncat
dim(E1b) = Ncat
dim(E1c) = Ncat
dim(E1d) = Ncat

# other summary stats that are calculated
dim(cumuInf_init) = Ncat
dim(alphaItot) = Ncat
dim(prev) = Ncat
# dim(B) <- c(Ncat, Ncat)

# FOI parameters
dim(lambda) = c(Ncat, Ncat)
dim(lambda_0) = c(Ncat, Ncat)
dim(lambda_1a) = c(Ncat, Ncat)
dim(lambda_1b) = c(Ncat, Ncat)
dim(lambda_1c) = c(Ncat, Ncat)
dim(lambda_1d) = c(Ncat, Ncat)


dim(lambda_sum_0) = Ncat
dim(lambda_sum_1a) = Ncat
dim(lambda_sum_1b) = Ncat
dim(lambda_sum_1c) = Ncat
dim(lambda_sum_1d) = Ncat

dim(c_comm_balanced) <- Ncat
dim(c_noncomm_balanced) <- Ncat
dim(theta) <- c(Ncat, Ncat)
dim(OnPrEP_init) = Ncat

dim(M_comm) = c(Ncat, Ncat)
dim(M_noncomm) = c(Ncat, Ncat)




dim(in_S0) <- c(Ncat, Ncat)
dim(in_S1a) <- c(Ncat, Ncat)
dim(in_S1b) <- c(Ncat, Ncat)
dim(in_S1c) <- c(Ncat, Ncat)
dim(in_S1d) <- c(Ncat, Ncat)

dim(in_I01) <- c(Ncat, Ncat)
dim(in_I11) <- c(Ncat, Ncat)
dim(in_I02) <- c(Ncat, Ncat)
dim(in_I03) <- c(Ncat, Ncat)
dim(in_I04) <- c(Ncat, Ncat)
dim(in_I05) <- c(Ncat, Ncat)

dim(in_I22) <- c(Ncat, Ncat)
dim(in_I23) <- c(Ncat, Ncat)
dim(in_I24) <- c(Ncat, Ncat)
dim(in_I25) <- c(Ncat, Ncat)
dim(in_I32) <- c(Ncat, Ncat)
dim(in_I33) <- c(Ncat, Ncat)
dim(in_I34) <- c(Ncat, Ncat)
dim(in_I35) <- c(Ncat, Ncat)
dim(in_I42) <- c(Ncat, Ncat)
dim(in_I43) <- c(Ncat, Ncat)
dim(in_I44) <- c(Ncat, Ncat)
dim(in_I45) <- c(Ncat, Ncat)

dim(rate_move_in) <- c(Ncat, Ncat)
dim(rate_move_out) <- Ncat
