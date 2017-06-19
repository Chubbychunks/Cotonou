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

