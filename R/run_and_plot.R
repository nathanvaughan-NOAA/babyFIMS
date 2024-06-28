obj <- MakeADFun(obj_fn, par,
                 map=map,
                 random=NULL,
                 silent=FALSE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))

opt$objective

names(opt$par)
opt$par

sdr <- sdreport(obj)
sdr
plr <- as.list(sdr,report=TRUE, "Est")
plrsd <- as.list(sdr,report=TRUE, "Std")

load('data/orig_am2022.Rdata')

Rec <- as.data.frame(arep$R)
names(Rec) <- c('year', 'R', 'Rec_sd', 'Rec_lci', 'Rec_uci')
Rec <- Rec %>%  mutate(version = 'amak') %>%
  bind_rows(data.frame(year = dat$year,
                       R = exp(plr$predlogR)/1e6,
                       Rec_uci = exp(plr$predlogR+2*plrsd$predlogR),
                       Rec_lci = exp(plr$predlogR-2*plrsd$predlogR),
                       version = 'babyFIMS'))
#View(Rec)
ggplot(Rec %>% filter(year >= 1978), aes(year, (R), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()


ssb <- as.data.frame(arep$SSB)
names(ssb) <- c('year', 'ssb', 'ssb_sd', 'ssb_lci', 'ssb_uci')
ssb <- ssb %>%
  mutate(version = 'amak') %>%
  bind_rows(data.frame(year = dat$year,
                       ssb = exp(plr$logssb),
                       ssb_uci = exp(plr$logssb+2*plrsd$logssb),
                       ssb_lci = exp(plr$logssb-2*plrsd$logssb),
                       version = 'babyFIMS'))

ggplot(ssb %>% filter(year >= 1978), aes(year, (ssb), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()
