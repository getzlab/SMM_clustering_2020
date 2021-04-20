## trajectory figures

## run data steps in 00-fig-km.R script first

## write to file
write <- TRUE
## calculate change from baseline
delta <- TRUE

out_folder <- '~/desktop/smm-out'
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

## todo figure with just prog (or non) with risk group (2 or 3)
pts <- c('Non progressor', 'Progressor')[1:2]

## two-year plots
xmax <- 2
xmax <- NULL

group <- 'cluster3'
group <- 'ttp_cat'

cc <- rawr::rawr_palettes$dfci[2:4]

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

ser <- merge(
  ser, dat[, c('id', 'ttp_time', 'ttp_ind')]
)


tmp <- ser[, grepl('(?i)^id|spike|hgb|flc', names(ser))]
if (delta)
  tmp[, -1] <- tmp[, -1] - tmp[, -1][rep(1:3, 8)]
mm <- rawr::melt(tmp, split(2:ncol(tmp), rep(1:3, 8)))
names(mm) <- c('id', 'time', 'ms', 'flc', 'hgb')

mm <- within(mm, {
  time2 <- c(0, 6, 12, 18, 24, 36, 48, 60)[time]
  time2 <- if (names(tte_factor) == 'Years')
    time2 / 12 else time2
})

library('ggplot2')
mm <- merge(
  mm,
  ser[, c('id', 'cluster', 'cluster2', 'cluster3', 'mayo_cat', 'ttp_time', 'ttp_ind')],
  all.x = TRUE
)

mm <- within(mm, {
  ## center time at progression/censor date? (in months)
  time3 <- time2 - ttp_time
  ttp_cat <- factor(ttp_ind, 0:1, c('Non progressor', 'Progressor'))
})


vars <- c('M-spike' = 'ms', 'Hemoglobin' = 'hgb', 'FLC ratio' = 'flc')
mm2 <- rawr::melt(mm, varying = vars)
mm2 <- within(mm2, {
  variable <- factor(variable, vars, names(vars))
})

## subset patients up to x years
if (!is.null(xmax))
  mm2 <- mm2[mm2$time2 <= xmax, ]

## subset patients by progression status
mm2 <- mm2[mm2$ttp_cat %in% pts, ]


library('patchwork')

if (write)
  pdf(file.path(out_folder, 'smm-labs.pdf'), height = 10, width = 14)
p1 <- ggplot(mm2, aes(time2, value, color = variable)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'Overall trajectories by cluster')
p1 <- p1 + facet_grid(variable ~ cluster, scales = 'free')


p2 <- ggplot(mm2, aes(time2, value, color = variable)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High vs low risk clusters')
p2 <- p2 + facet_grid(variable ~ cluster2, scales = 'free')


p3 <- ggplot(mm2, aes(time2, value, color = variable)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High, low, intermediate risk clusters')
p3 <- p3 + facet_grid(variable ~ cluster3, scales = 'free')

print(p1 / (p2 + p3))


## by group indicator
p1 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  scale_linetype_manual(values = 1:3) +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'Overall trajectories by cluster')
p1 <- p1 + facet_grid(variable ~ cluster, scales = 'free')


p2 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  scale_linetype_manual(values = 1:3) +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High vs low risk clusters')
p2 <- p2 + facet_grid(variable ~ cluster2, scales = 'free')


p3 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_smooth(method = 'loess') +
  scale_linetype_manual(values = 1:3) +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High, low, intermediate risk clusters')
p3 <- p3 + facet_grid(variable ~ cluster3, scales = 'free')

print(p1 / (p2 + p3))

if (write)
  dev.off()



if (write)
  pdf(file.path(out_folder, 'smm-labs-indiv.pdf'), height = 10, width = 14)
p1 <- ggplot(mm2, aes(time2, value, color = variable, group = id)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'Individual trajectories by cluster')
p1 <- p1 + facet_grid(variable ~ cluster, scales = 'free')


p2 <- ggplot(mm2, aes(time2, value, color = variable, group = id)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High vs low risk clusters')
p2 <- p2 + facet_grid(variable ~ cluster2, scales = 'free')


p3 <- ggplot(mm2, aes(time2, value, color = variable, group = id)) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High, low, intermediate risk clusters')
p3 <- p3 + facet_grid(variable ~ cluster3, scales = 'free')

print(p1 / (p2 + p3))


## by group indicator
p1 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group,
                             group = 'id')) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'Individual trajectories by cluster')
p1 <- p1 + facet_grid(variable ~ cluster, scales = 'free')


p2 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group,
                             group = 'id')) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High vs low risk clusters')
p2 <- p2 + facet_grid(variable ~ cluster2, scales = 'free')


p3 <- ggplot(mm2, aes_string('time2', 'value', color = 'variable', linetype = group,
                             group = 'id')) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = cc) +
  theme_classic() +
  geom_line() +
  labs(x = names(tte_factor), y = '', color = '', linetype = '',
       subtitle = 'High, low, intermediate risk clusters')
p3 <- p3 + facet_grid(variable ~ cluster3, scales = 'free')

print(p1 / (p2 + p3))

if (write)
  dev.off()


## compare changes at 2, 5 years for high vs int vs low
tt <- c('one', 'two', 'five')
tmp <- merge(tmp, ser[, c('id', 'cluster', 'cluster2', 'cluster3')])
names(tmp) <- gsub('twelve_months', 'one_year', names(tmp))

if (write)
  pdf(file.path(out_folder, 'delta.pdf'), height = 9, width = 12)
par(mfrow = c(length(tt), 3L), oma = c(0, 4, 2, 1), mar = c(4, 5, 3, 1))
lapply(tt, function(x) {
  dd <- tmp[, c(grep(x, names(tmp), value = TRUE), 'cluster3')]
  
  msp <- grep('mspike', names(dd), value = TRUE)
  flc <- grep('flc', names(dd), value = TRUE)
  hgb <- grep('hgb', names(dd), value = TRUE)
  
  tp <- function(formula, data, ...) {
    rawr::tplot(
      formula, data, ..., bty = 'l', test = kruskal.test, xlab = '',
      show.n = FALSE, las = 1L,
      args.test = list(cex = 1, at = 0.5, line = 1)
    )
    # mtext('Overall', line = 2, adj = 0.5, at = par('usr')[1L])
    cz <- rawr::cuzick.test(formula, data)$details$pairs$p.value
    # bp.test(pvalr(cz, show.p = TRUE))
    bp.test(formula, data)
    invisible(NULL)
  }
  
  tp(
    reformulate('cluster3', msp), dd, ylab = 'M-spike', boxcol = cc[1L]
  )
  title(ylab = sprintf('%s-year', rawr::case(x)), line = 6,
        font.lab = 2L, cex.lab = 1.5, xpd = NA)
  tp(
    reformulate('cluster3', hgb), dd, ylab = 'Hemoglobin', boxcol = cc[2L]
  )
  tp(
    reformulate('cluster3', flc), dd, ylab = 'FLC ratio', boxcol = cc[3L]
  )
})
if (write)
  dev.off()
