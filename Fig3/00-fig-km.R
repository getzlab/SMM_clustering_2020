## km, forest plots


## write to file
write <- TRUE

out_folder <- '~/desktop/smm-out'
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

library('survival')
library('rawr')

## pts in dat but not in serial samples
ids <- c("SMM84", "SMM83", "SMM97")
ids <- ''

tte_factor <- c(Days = 1)
tte_factor <- c(Months = 30.437)
tte_factor <- c(Years = 365.242)


## cluster colors
cc <- c(
  C1 = '#0073C2FF',
  C2 = '#EFC000FF',
  C3 = '#868686FF',
  C4 = '#CD534CFF',
  C5 = '#7AA6DCFF',
  C6 = '#003C67FF'
)
cc2 <- c('blue4', 'goldenrod3', 'tomato3')

p <- './inst/manu/SMM_cluster_paper/Analysis files/'

ser <- data.table::fread(file.path(p, 'Clusters_serial_timepoints_updated_feb21.txt'))
ser <- as.data.frame(ser)

dat <- data.table::fread(file.path(p, 'SMM_progression_analysis_final_updated_August19_Mark_august20.txt'))
dat <- as.data.frame(dat)
names(dat) <- make.unique(tolower(names(dat)))


dat <- dat[!dat$`sample id` %in% ids, ]


dat <- within(dat, {
  id <- `sample id`
  
  cluster <- factor(clusters, c('C4', 'C6', 'C1', 'C2', 'C3', 'C5'))
  
  cluster1 <- factor(clusters, paste0('C', 1:6))
  
  cluster2 <- combine_levels(
    cluster,
    list(low = c('C4', 'C1', 'C6'), hi = c('C2', 'C3', 'C5'))
  )
  cluster2_cat <- combine_levels(
    cluster,
    list('C1, C4, C6' = c('C4', 'C1', 'C6'), 'C2, C3, C5' = c('C2', 'C3', 'C5'))
  )
  
  cluster3 <- combine_levels(
    cluster,
    list(low = 'C4', mid = c('C1', 'C6'), hi = c('C2', 'C3', 'C5'))
  )
  cluster3_cat <- combine_levels(
    cluster,
    list(C4 = 'C4', 'C1, C6' = c('C1', 'C6'), 'C2, C3, C5' = c('C2', 'C3', 'C5'))
  )
  
  ttp_time <- ttp / tte_factor
  ttp_ind <- progression_status
  
  mayo_cat <- factor(new_mayo_model)
  mayo_cat2 <- factor(new_mayo_model, 1:3, paste(c('Low', 'Intermediate', 'High'), 'Risk'))
})


cx <- coxph(Surv(ttp_time, ttp_ind) ~ cluster + mayo_cat, dat)
broom::tidy(cx, exponentiate = TRUE, conf.int = TRUE)

cx <- coxph(Surv(ttp_time, ttp_ind) ~ cluster, dat)
broom::tidy(cx, exponentiate = TRUE, conf.int = TRUE)

cx <- coxph(Surv(ttp_time, ttp_ind) ~ cluster, dat)
broom::tidy(cx, exponentiate = TRUE, conf.int = TRUE)


m <- c(11, 9, 2, 2)
kmby <- function(s, e, d, ...) {
  yat <- 0:2 / 2
  rawr::kmplot_by(
    s, e, d, ..., sub = '', lwd.surv = 2, strata_lab = FALSE, mar = m,
    yaxis.at = yat, yaxis.lab = yat * 100,
    bty = 'l', test_details = FALSE,
    xlab = 'Years',
    ylab = switch(
      e,
      os = 'Percent alive',
      pfs = 'Percent alive and\nfree from progression',
      ttp = 'Percent free from progression'
    )
  )
}

if (write)
  pdf(file.path(out_folder, 'fig3-ann.pdf'), height = 6, width = 8)
kmby('cluster1', 'ttp', dat, col.surv = cc)
kmby('cluster2_cat', 'ttp', dat, hr_text = TRUE, col.surv = cc2[-2])
kmby('cluster3_cat', 'ttp', dat, pw_test = TRUE, hr_text = TRUE, col.surv = cc2)
if (write)
  dev.off()

if (write)
  pdf(file.path(out_folder, 'fig3-clean.pdf'), height = 6, width = 8)
kmby('cluster1', 'ttp', dat, col.surv = cc)
kmby('cluster2_cat', 'ttp', dat, col.surv = cc2[-2])
kmby('cluster3_cat', 'ttp', dat, col.surv = cc2)
if (write)
  dev.off()


library('forest')
cx0 <- coxph(Surv(ttp_time, ttp_ind) ~ cluster, dat)
cx1 <- coxph(Surv(ttp_time, ttp_ind) ~ cluster + mayo_cat2, dat)
cx2 <- coxph(Surv(ttp_time, ttp_ind) ~ cluster2_cat + mayo_cat2, dat)

if (write)
  pdf(file.path(out_folder, 'fig3-forest0.pdf'), height = 4, width = 10)
f0 <- forest(
  cx0, header = 'Genetic subgroups',
  plotArgs = list(
    cex = 2, labels = levels(dat$cluster),
    show_conf = TRUE, layout = 'unified', xlim = c(0, 10), reset_par = FALSE,
    names = c('Characteristic', 'N (%)', 'HR (95% CI)', 'p-value')
  )
)
title(xlab = 'Hazard ratio')
if (write)
  dev.off()

if (write)
  pdf(file.path(out_folder, 'fig3-forest1.pdf'), height = 5, width = 10)
f1 <- forest(
  cx1, header = c('Genetic subgroups', '20-2-20 clinical model'),
  plotArgs = list(
    cex = 2, labels = c(levels(dat$cluster), levels(dat$mayo_cat2)),
    show_conf = TRUE, layout = 'unified', xlim = c(0, 10), reset_par = FALSE,
    names = c('Characteristic', 'N (%)', 'HR (95% CI)', 'p-value')
  )
)
title(xlab = 'Hazard ratio')
if (write)
  dev.off()

if (write)
  pdf(file.path(out_folder, 'fig3-forest2.pdf'), height = 4, width = 10)
f2 <- forest(
  cx2, header = c('Low- vs high-risk genetic subgroups', '20-2-20 clinical model'),
  plotArgs = list(
    cex = 2, labels = c(levels(dat$cluster2_cat), levels(dat$mayo_cat2)),
    show_conf = TRUE, layout = 'unified', xlim = c(0, 10), reset_par = FALSE,
    names = c('Characteristic', 'N (%)', 'HR (95% CI)', 'p-value')
  )
)
title(xlab = 'Hazard ratio')
if (write)
  dev.off()


summary.forest <- function(x, digits = 2L) {
  txt <- do.call('cbind', x$cleanfp_list$text)
  txt[is.na(txt)] <- ''
  txt[, 2L] <- apply(txt[, -1L], 1, function(x) paste(x, collapse = ' - '))
  txt[, 3L] <- x$cleanfp_list$`p-value`
  txt[is.na(txt) | txt == ' - '] <- ''
  colnames(txt) <- c('HR', '95% CI', 'p-value')
  term <- gsub('^ +', '&emsp;', x$cleanfp_list$Term)
  term[!grepl('^&em', term)] <- sprintf('<b>%s</b>', term[!grepl('^&em', term)])
  cbind(Term = term, txt)
}


s0 <- summary(f0)
s0 <- gsub('cluster', '', s0)
ht <- htmlTable::htmlTable(
  s0, align = 'lccc', caption = 'Univariable model: cluster.'
)
structure(ht, class = 'htmlTable')
if (write)
  rawr::write_htmlTable(ht, file = file.path(out_folder, '3.html'))


s1 <- summary(f1)
s1 <- gsub('cluster|mayo_cat2', '', s1)
ht <- htmlTable::htmlTable(
  s1, align = 'lccc', caption = 'Multivariable model: cluster, Mayo.'
)
structure(ht, class = 'htmlTable')
if (write)
  rawr::write_htmlTable(ht, file = file.path(out_folder, '3a.html'))

s2 <- summary(f2)
s2 <- gsub('cluster2_cat|mayo_cat2', '', s2)
ht <- htmlTable::htmlTable(
  s2, align = 'lccc', caption = 'Multivariable model: cluster risk, Mayo.'
)
structure(ht, class = 'htmlTable')
if (write)
  rawr::write_htmlTable(ht, file = file.path(out_folder, '3d.html'))



## repeat with 3 levels for cluster
cx2 <- coxph(Surv(ttp_time, ttp_ind) ~ cluster3_cat + mayo_cat2, dat)

if (write)
  pdf(file.path(out_folder, 'fig3-forest-low-high.pdf'), height = 4.25, width = 10)
f2 <- forest(
  cx2, header = c('Low- vs high-risk genetic subgroups', '20-2-20 clinical model'),
  plotArgs = list(
    cex = 2, labels = c(levels(dat$cluster3_cat), levels(dat$mayo_cat2)),
    show_conf = TRUE, layout = 'unified', xlim = c(0, 10), reset_par = FALSE,
    names = c('Characteristic', 'N (%)', 'HR (95% CI)', 'p-value')
  )
)
title(xlab = 'Hazard ratio')
if (write)
  dev.off()

s2 <- summary(f2)
s2 <- gsub('cluster3_cat|mayo_cat2', '', s2)
ht <- htmlTable::htmlTable(
  s2, align = 'lccc', caption = 'Multivariable model: cluster risk, Mayo.'
)
structure(ht, class = 'htmlTable')
if (write)
  rawr::write_htmlTable(ht, file = file.path(out_folder, '3d-low-high.html'))
