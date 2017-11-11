
res = lapply(result[[2]], cotonou::return_outputs, cotonou::main_model, time = time, outputs = outputs)


# ignore these ######################################
frac_ProFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,1])), 2, cotonou::quantile_95)))
frac_LowFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,2])), 2, cotonou::quantile_95)))
frac_GPF = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,3])), 2, cotonou::quantile_95)))
frac_FormerFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,4])), 2, cotonou::quantile_95)))
frac_Client = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,5])), 2, cotonou::quantile_95)))
frac_GPM = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,6])), 2, cotonou::quantile_95)))
frac_Virgin_Female = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,7])), 2, cotonou::quantile_95)))
frac_Virgin_Male = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,8])), 2, cotonou::quantile_95)))
frac_Former_FSW_Outside = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res, function(x) x$frac_N*100), function(x) x[,9])), 2, cotonou::quantile_95)))

frac_Active_FSW = data.frame(time, t(apply(do.call(rbind, lapply(res, function(x) {100*(x$frac_N[,1] + x$frac_N[,2])})), 2, cotonou::quantile_95)))
Ratio_Low_Pro = data.frame(time, t(apply(do.call(rbind, lapply(res, function(x) {x$frac_N[,2]/ x$frac_N[,1]})), 2, cotonou::quantile_95)))



frac = rbind(frac_ProFSW, frac_LowFSW, frac_GPF, frac_FormerFSW, frac_Client, frac_GPM, frac_Virgin_Female, frac_Virgin_Male, frac_Former_FSW_Outside, frac_Active_FSW, Ratio_Low_Pro)
frac = data.frame(frac, group = rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"), each = length(time)))
colnames(frac) = c("time", "Lower", "Median", "Upper", "variable")
frac$variable = factor(frac$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"))

prev_FSW = t(apply(do.call(rbind, lapply(res, function(x) x$prev_FSW)), 2, cotonou::quantile_95))
prev_LowFSW = t(apply(do.call(rbind, lapply(res, function(x) x$prev_LowFSW)), 2, cotonou::quantile_95))
prev_client = t(apply(do.call(rbind, lapply(res, function(x) x$prev_client)), 2, cotonou::quantile_95))
prev_women = t(apply(do.call(rbind, lapply(res, function(x) x$prev_women)), 2, cotonou::quantile_95))
prev_men = t(apply(do.call(rbind, lapply(res, function(x) x$prev_men)), 2, cotonou::quantile_95))
prev = rbind(prev_FSW, prev_LowFSW, prev_client, prev_women, prev_men)
prev = data.frame(time, prev, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)))
colnames(prev) = c("time", "Lower", "Median", "Upper", "variable")
prev$variable = factor(prev$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))


prev_FSW_indiv = t(do.call(rbind, lapply(res, function(x) x$prev_FSW)))
prev_LowFSW_indiv = t(do.call(rbind, lapply(res, function(x) x$prev_LowFSW)))
prev_client_indiv = t(do.call(rbind, lapply(res, function(x) x$prev_client)))
prev_women_indiv = t(do.call(rbind, lapply(res, function(x) x$prev_women)))
prev_men_indiv = t(do.call(rbind, lapply(res, function(x) x$prev_men)))
prev_indiv = rbind(prev_FSW_indiv, prev_LowFSW_indiv, prev_client_indiv, prev_women_indiv, prev_men_indiv)

prev_indiv = data.frame(time, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)), prev_indiv)


colnames(prev_indiv) = c("time", "variable", as.character(seq(1, length(prev_FSW_indiv[1,]))))

prev_indiv_melted = reshape2::melt(prev_indiv, id.vars = c("time", "variable"))

colnames(prev_indiv_melted) = c("time", "variable", "run", "value")

prev_indiv_melted$variable = factor(prev_indiv_melted$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))





Ntot = data.frame(time, t(apply(do.call(rbind, lapply(res, function(x) x$Ntot)), 2, cotonou::quantile_95)))
colnames(Ntot) = c("time", "Lower", "Median", "Upper")

ART_coverage_women = t(apply(do.call(rbind, lapply(res, function(x) x$ART_coverage_women)), 2, cotonou::quantile_95))
ART_coverage_men = t(apply(do.call(rbind, lapply(res, function(x) x$ART_coverage_men)), 2, cotonou::quantile_95))
ART_coverage_FSW = t(apply(do.call(rbind, lapply(res, function(x) x$ART_coverage_FSW)), 2, cotonou::quantile_95))
ART_coverage_all = t(apply(do.call(rbind, lapply(res, function(x) x$ART_coverage_all)), 2, cotonou::quantile_95))
ART_coverage = rbind(ART_coverage_women, ART_coverage_men, ART_coverage_FSW, ART_coverage_all)
ART_coverage = data.frame(time, ART_coverage, rep(c("Women", "Men", "Pro FSW", "All"), each = length(time)))
colnames(ART_coverage) = c("time", "Lower", "Median", "Upper", "variable")
ART_coverage$variable = factor(ART_coverage$variable, levels = c("Pro FSW", "Women", "Men", "All"))
ART_coverage = ART_coverage[ART_coverage$variable == "All" | ART_coverage$variable == "Pro FSW",]


frac_N_discard_points_graph = frac_N_discard_points
frac_N_discard_points_graph[frac_N_discard_points_graph$variable == "Low Pro Ratio", c(2,3)] = frac_N_discard_points_graph[frac_N_discard_points_graph$variable == "Low Pro Ratio",c(2,3)]/100

#####################################################

require(ggplot2)

# plot fraction in each group
ggplot(frac) + geom_line(aes(x = time, y = Median)) + geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) +
  theme_bw() + labs(y = "Percent in each group (%)") +  facet_wrap(~variable, scales = "free") +
  geom_point(data = frac_N_data_points, aes(x = time, y = point), size = I(2), color = "red", shape = 15) +
  geom_hline(data = frac_N_discard_points_graph, aes(yintercept = 100*min), size = I(0.5), color = "red", linetype = 1, alpha = 0.7) +
  geom_hline(data = frac_N_discard_points_graph, aes(yintercept = 100*max), size = I(0.5), color = "red", linetype = 1, alpha = 0.7)

prev_axes = data.frame(variable = c(rep("Pro FSW", 2),
                                    rep("Clients", 2),
                                    rep("Women", 2),
                                    rep("Men", 2),
                                    rep("Low-level FSW", 2)),
                       time = c(rep(c(1986, 2025), 5)),
                       value = c(0, 70, 0, 70, 0, 15, 0, 15, 0, 70)
)

prev_points_80s = prev_points_all[c(1,2,3),]

# plot prevalence in each group
ggplot() + geom_line(data = prev, aes(x = time, y = Median))+ geom_ribbon(data = prev, aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "prevalence (%)") +
  geom_point(data = prev_points, aes(x = time, y = value))+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper)) +
  geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "red")+
  geom_blank(data = prev_axes, aes(x = time, y = value))


# # plot prevalence in each group indiv runs
# ggplot() + geom_line(data = prev_indiv_melted, aes(x = time, y = value, factor = variable, factor = run), alpha = 0.3) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "prevalence (%)") +
#   geom_point(data = prev_points, aes(x = time, y = value))+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper))+
#   geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "red")+
#   geom_blank(data = prev_axes, aes(x = time, y = value))



# plot total population size
ggplot(Ntot) + geom_line(aes(x = time, y = Median)) + geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) +
  theme_bw() + labs(y = "Total population size of Grand Cotonou") +
  geom_point(data = Ntot_data_points, aes(x = time, y = point, color = colour), size = I(2), shape = 15) + geom_errorbar(data = Ntot_data_points, aes(x = time, ymax = upper, ymin = lower, color = colour), width = 2)


ggplot(ART_coverage) +
  geom_line(aes(x = time, y = Median))+ geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() +
  facet_wrap(~variable) + labs(y = "ART coverage ") +
  geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower, ymax = Upper), colour = "darkred")
