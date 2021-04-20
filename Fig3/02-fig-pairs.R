## eMP and eHgb, which are definitions of evolving M protein or Hgb over a
## one year period. Could you do a logistic regression for this on the cluster
## risks, either high vs Int vs low or just high and low as we want to make
## the point that these genetic subtypes could be predictive of evolving M protein
## or Hemoglobin


## run data steps in 00-fig-km.R thne 01-fig-traj scripts first

## write to file
write <- TRUE

out_folder <- '~/desktop/smm-out'
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

library('forest')

xx <- c(
  'Evolving M protein' = 'eMP',
  'Evolving hemoglobin' = 'eHb'
)
yy <- 'cluster2'

ser <- within(ser, {
  id <- `Sample ID`
  cluster <- factor(Cluster, c('C4', 'C6', 'C1', 'C2', 'C3', 'C5'))
  
  cluster2 <- combine_levels(
    cluster,
    list('Low risk' = c('C4', 'C1', 'C6'), 'High risk' = c('C2', 'C3', 'C5'))
  )
  
  cluster3 <- combine_levels(
    cluster,
    list('Low risk' = 'C4', Intermediate = c('C1', 'C6'), 'High risk' = c('C2', 'C3', 'C5'))
  )
  
  mayo_cat <- factor(New_mayo_model)
})

tmp <- merge(dat, ser[, c('id', xx)], all = TRUE)

models <- list(
  eMP = glm(eMP ~ cluster2, tmp, family = binomial('logit')),
  eMP = glm(eMP ~ cluster3, tmp, family = binomial('logit')),
  eHb = glm(eHb ~ cluster2, tmp, family = binomial('logit')),
  eHb = glm(eHb ~ cluster3, tmp, family = binomial('logit'))
)

prep_lists <- lapply(models, forest, plot = FALSE)
prep_lists <- lapply(prep_lists, function(x)
  `class<-`(x[[1L]], 'cleanfp_list'))
x <- Reduce(merge_forest, prep_lists)

group.col <- rep_len(c('grey95', 'none'), length(models))
group.col <- rep(group.col, sapply(prep_lists, function(x) length(x$Term)))

lbl <- rep(c(levels(tmp$cluster2), levels(tmp$cluster3)), 2)
lbl <- factor(lbl, c('low', 'mid', 'hi'), c('Low-risk', 'Intermediate', 'High-risk'))
lbl <- paste0(lbl, strrep(' ', seq_along(lbl)))

if (write)
  pdf(file.path(out_folder, 'logist-evolving.pdf'), height = 5, width = 10)
par(oma = c(0, 2, 0, 0))
plot(
  x, col.rows = group.col, reset_par = FALSE, cex = 2, layout = 'unified',
  labels = lbl, xlim = c(0, 30),
  show_conf = TRUE, names = c('Cluster', 'N (%)', 'OR (95% CI)', 'p-value')
)
rl <- rev(rle(group.col)$lengths)
text(grconvertX(0.025, 'ndc'), rev(cumsum(head(c(0, rl), -1)) + rl / 2) + 0.5,
     names(models), xpd = NA, srt = 90, adj = 0.5)
title(xlab = 'Odds ratio')
box('outer')
if (write)
  dev.off()

## boxplots
f <- function(data, x) {
  tmp <- droplevels(data[data$variable %in% x, ])
  sp <- split(tmp$value, tmp$time)
  res <- lapply(sp[-1L], function(x) {
    wil <- wilcox.test(sp[[1L]], x)
    con <- sort(sp[[1L]])
    com <- sort(x)
    
    data.frame(
      check.names = FALSE,
      'Median (range)' = sprintf('%.2f (%.2f, %.2f)', mean(con), min(con), max(con)),
      'Median (range)' = sprintf('%.2f (%.2f, %.2f)', mean(com), min(com), max(com)),
      'p-value' = wil$p.value
    )
  })
  
  do.call('rbind', res)
}


tmp <- mm2[mm2$cluster3 %in% 'Intermediate', ]
tmp <- mm2
sp <- split(tmp, tmp$variable)

cc <- adjustcolor(rawr::rawr_palettes$dfci[2:4], 0.75)

if (write)
  pdf(file.path(out_folder, 'pairs-all.pdf'), height = 10, width = 12)
par(mfrow = c(3, 1), mar = c(4, 5, 2, 1), oma = c(4, 0, 3, 0))
rawr::tplot(
  value ~ factor(time, 1:8, c(0, 6, 12, 18, 24, 36, 48, 60)),
  sp$`M-spike`, las = 1L, xlab = '', ylab = 'M-spike',
  type = 'db', quantiles = c(0.25, 0.5, 0.75), lwd = c(1, 3, 1),
  show.n = FALSE, bty = 'l', boxcol = cc[1L]
)
pv <- f(tmp, 'M-spike')
# text(1:8, -1, gsub(' .*', '', c(pv[1L, 2L], pv[, 2L])), xpd = NA)
text(
  1:8, grconvertY(-0.2, 'npc'),
  c('Reference', pvalr(pv[, 3L], show.p = TRUE)),
  xpd = NA
)

rawr::tplot(
  value ~ factor(time, 1:8, c(0, 6, 12, 18, 24, 36, 48, 60)),
  sp$`Hemoglobin`, las = 1L, xlab = '', ylab = 'Hemoglobin',
  type = 'db', quantiles = c(0.25, 0.5, 0.75), lwd = c(1, 3, 1),
  show.n = FALSE, bty = 'l', boxcol = cc[2L]
)
pv <- f(tmp, 'Hemoglobin')
# text(1:8, -1, gsub(' .*', '', c(pv[1L, 2L], pv[, 2L])), xpd = NA)
text(
  1:8, grconvertY(-0.2, 'npc'),
  c('Reference', pvalr(pv[, 3L], show.p = TRUE)),
  xpd = NA
)

rawr::tplot(
  value ~ factor(time, 1:8, c(0, 6, 12, 18, 24, 36, 48, 60)),
  sp$`FLC ratio`, las = 1L, xlab = '', ylab = 'FLC ratio',
  type = 'db', quantiles = c(0.25, 0.5, 0.75), lwd = c(1, 3, 1),
  show.n = FALSE, bty = 'l', boxcol = cc[3L]
)
pv <- f(tmp, 'FLC ratio')
# text(1:8, -1, gsub(' .*', '', c(pv[1L, 2L], pv[, 2L])), xpd = NA)
text(
  1:8, grconvertY(-0.2, 'npc'),
  c('Reference', pvalr(pv[, 3L], show.p = TRUE)),
  xpd = NA
)
title(xlab = 'Months', line = 6, xpd = NA, cex.lab = 1.5)
text(
  grconvertX(0.025, 'ndc'), grconvertY(0.975, 'ndc'),
  'All patients', xpd = NA, cex = 1.5, adj = 0
)
if (write)
  dev.off()



tmp0 <- split(mm2, mm2$variable)

if (write)
  pdf(file.path(out_folder, 'pairs-risk.pdf'), height = 8, width = 12)

par(mfrow = c(3, 1), mar = c(4, 5, 2, 1), oma = c(2, 0, 3, 0))
for (ii in seq_along(levels(mm2$variable))) {
  var <- levels(mm2$variable)[ii]
  tmp <- tmp0[[var]]
  sp <- split(tmp, list(tmp$cluster2, factor(tmp$time, 1:8, c(0, 6, 12, 18, 24, 36, 48, 60))))
  
  rawr::tplot(
    value ~ cluster2 + factor(time, 1:8, c(0, 6, 12, 18, 24, 36, 48, 60)),
    at = seq_along(sp) + c(0.15, -0.15), axes = FALSE, names = '', xlim = c(1, 16.5),
    # names = gsub('\\..*', '', names(sp)),
    tmp, las = 1L, xlab = '', ylab = var,
    type = 'db', quantiles = c(0.25, 0.5, 0.75), lwd = c(1, 3, 1),
    show.n = FALSE, bty = 'l', boxcol = tcol(cc[ii], c(1, 0.25))
  )
  xat <- seq_along(sp)[c(TRUE, FALSE)] + 0.5
  box(bty = 'l')
  axis(2L, las = 1L)
  axis(1L, xat, c(0, 6, 12, 18, 24, 36, 48, 60))
  legend(
    'right', legend = levels(tmp$cluster2), title = 'Cluster',
    fill = tcol(cc[ii], c(1, 0.25)), bty = 'n'
  )
  
  pv <- sapply(seq_along(sp)[c(TRUE, FALSE)], function(ii) {
    wilcox.test(sp[[ii]]$value, sp[[ii + 1L]]$value)$p.value
  })
  text(c(0.5, xat), par('usr')[4L],
       c('naive', pvalr(pv, show.p = TRUE)),
       xpd = NA, pos = 3L)
  text(c(0.5, xat), par('usr')[4L],
       c('corrected', pvalr(p.adjust(pv, method = 'BH'), show.p = TRUE)),
       xpd = NA, pos = 1L)
  if (ii == 3L)
    title(xlab = 'Months', xpd = NA)
}
box('outer')
if (write)
  dev.off()
