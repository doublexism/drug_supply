## read data
adex <- read.sas7bdat(file = "data/adex.sas7bdat")
adtte <- read.sas7bdat("data/adtteeff.sas7bdat")
adrs <- read.sas7bdat("data/adrs.sas7bdat")
adsl <- read.sas7bdat("data/adsl.sas7bdat")
## cleanning
# date handling
timevar_adex <- c("TRTSDT","TRTEDT","EXSTDT","EXENDT")
timevar_adtte <- c("TRTSDT","TRTEDT","STARTDT","ADT")
adex[timevar_adex] <- map(timevar_adex, ~as.Date(adex[[.]], origin = "1960-01-01"))
adtte[timevar_adtte] <- map(timevar_adtte, ~as.Date(adtte[[.]], origin = "1960-01-01"))
# exposure dataset
dat_expo <- adex %>% select(SUBJID, SITEID, AGE, SEX, TRTSDT, TRTEDT) %>% distinct 
# outcome dataset
dat_pfs <- adtte %>% filter(PARAMCD %in% c("PFSINVOU"))  %>% select(SUBJID,PFS = AVAL,ADT,CENSOR = CNSR,EVNTDESC)
dat_os <- adtte %>% filter(PARAMCD %in% c("OS")) %>% select(SUBJID,OS = AVAL)
# response dataset
dat_rs <- adrs %>% filter(PARAMCD == "BSTOVRLI") %>% select(SUBJID, BOR = AVALC)
# censoring
dat_wd <- adsl %>% select(SUBJID, DCTREASN)
## merging dataset
dat <- dat_expo %>% left_join(dat_pfs, by = "SUBJID") %>% 
  left_join(dat_os, by = "SUBJID") %>%
  left_join(dat_rs, by = "SUBJID") %>%
  left_join(dat_wd, by = "SUBJID") %>%
  mutate(PD = 1 - CENSOR, 
         DEATH = as.numeric(EVNTDESC == "Death"),
         EOT = as.numeric(!is.na(DCTREASN)),
         RESP = as.numeric(BOR %in% c("CR","PR")),
         TRTDUR = as.numeric(TRTEDT - TRTSDT+21)/30)%>%
  mutate(LFU = as.numeric(TRTDUR < PFS | (EOT == 1 & DCTREASN != 1)),
         PFS = ifelse(PD == 1 | EOT == 1, PFS, TRTDUR),
         OS = ifelse(DEATH == 1 | EOT == 1, OS, TRTDUR),
         timeInTrial = as.numeric(TRTSDT - min(TRTSDT))/30.4375,
         timeOutTrial = as.numeric(TRTEDT - min(TRTSDT))/30.4375,
         enrolled = 1)

dat_PD <- dat %>% filter(PD == 1 & TRTDUR >= PFS) %>% mutate(TRTDUR = TRTDUR - PFS)
dat_NonPD <- dat %>% filter(!SUBJID %in% dat_PD$SUBJID)
## estimate lambda for PD
PD_pois <- survfit(Surv(TRTDUR, EOT)~1, data = dat) %>% summary(2:15)
PD_pois <- data.frame(n.risk = PD_pois$n.risk+ PD_pois$n.event +PD_pois$n.censor, 
                      num_event = PD_pois$n.event, 
                      num_censor = PD_pois$n.censor,
                      time = PD_pois$time)
enrolled <- survfit(Surv(timeInTrial, enrolled) ~ 1, data = dat) %>% summary(1:12) %>% `[[`("n.event") %>% cumsum()
enrolled[7] <- enrolled[7]+1
enrolled <- c(enrolled, rep(enrolled[7], 6))
events <- survfit(Surv(timeOutTrial, EOT) ~ 1, data = dat) %>% summary(1:14) %>% `[[`("n.event") %>% cumsum()
num_at_risk <- enrolled - events


## prior setting
med_Brent <- medianTTE(percent = 0.22, time = 5,n = 102, unit = "y")
med_nivo1 <- medianTTE(percent = 0.86, time = 24, n = 23,unit = "w")
med_nivo2 <- medianTTE(percent = 0.77, time = 6,n = 80, unit = "m")
med_pembro <- medianTTE(percent = 0.69, time = 24, n=31,unit = "w")
hist <- list(med_Brent, med_nivo1, med_nivo2, med_pembro)
## thus set prior distribution of log-hazard to 
hist_mean <- mean(rep(getListElement(hist,"loghazard"), getListElement(hist, "n")))
hist_dev <- sd(rep(getListElement(hist,"loghazard"), getListElement(hist, "n")))
bayes.pois <- MCMCpoisson(num_event ~ log(n.risk), b0 = c(hist_mean, 1), B0 = c(1/(hist_dev)**2, 1e8),data = PD_pois, burnin = 10000, mcmc = 1e6, thin = 10)

## prediction without data
project_time_frame_null <- 1:25
enrolled_null <-c(rep(10,7),rep(0,18)) 
hazards_null <- rnorm(100, hist_mean, hist_dev)%>% exp()
preds <- dropoutMCMC(exposure =rep(1,10), duration = length(project_time_frame_null),enroll_pred = enrolled_null, hazards = hazards_null,num_sim = 500)
drop_trace_null <- do.call(cbind, preds$at.risk) %>% 
  apply(1, quantile,probs = c(0.025, 0.50, 0.975)) %>% 
  round

ylow_null <- c(10, drop_trace_null[1,])
ymid_null <- c(10, drop_trace_null[2,])
yhigh_null <- c(10, drop_trace_null[3,])
yobserved <- num_at_risk[1:13]
ggplot() + 
  geom_line(aes(x = project_time_frame_null, y = ylow_null), linetype = "dotted", size = 1) +
  geom_line(aes(x = project_time_frame_null, y = ymid_null), linetype = "dashed", size = 1) +
  geom_line(aes(x = project_time_frame_null, y = yhigh_null), linetype = "dotted", size = 1) +
  geom_line(aes(x = 1:13, y = yobserved), size = 1, color = "red") +
  xlab("Months since First patient dosing") +
  ylab("Number of patients on treatment,95% CI") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_y_continuous(breaks = seq(0,70, 10)) +
  theme_igray()

## prediction
project_time_frame <- 13:25
hazards <- sample(bayes.pois[,1],100) %>% exp()
preds <- dropoutMCMC(exposure =dat$TRTDUR[dat$EOT == 0], duration = length(project_time_frame)-1,enroll_pred = NULL, hazards = hazards,num_sim = 100)
drop_trace <- do.call(cbind, preds$at.risk) %>% 
  apply(1, quantile,probs = c(0.025, 0.50, 0.975)) %>% 
  round

ylow <- c(num_at_risk[13], drop_trace[1,])
ymid <- c(num_at_risk[13], drop_trace[2,])
yhigh <- c(num_at_risk[13], drop_trace[3,])
ggplot() + 
  geom_line(aes(x = 1:13, y = num_at_risk),size = 1.3) +
  geom_line(aes(x = project_time_frame, y = ylow), linetype = "dotted", size = 1) +
  geom_line(aes(x = project_time_frame, y = ymid), linetype = "dashed", size = 1) +
  geom_line(aes(x = project_time_frame, y = yhigh), linetype = "dotted", size = 1) +
  geom_vline(xintercept = 13, color = "red") +
  xlab("Months since First patient dosing") +
  ylab("Number of patients on treatment,95% CI") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_y_continuous(breaks = seq(0,70, 10)) +
  theme_igray()

## 6 months data prediction
dat_6 <- dat[dat$timeInTrial < 7,]
dat_6$EOT[dat_6$TRTDUR + dat_6$timeInTrial > 7] <- 0
dat_6$TRTDUR[dat_6$TRTDUR > 7] <- 7
PD_pois_6 <- survfit(Surv(TRTDUR, EOT)~1, data = dat_6) %>% summary(2:7)
PD_pois_6 <- data.frame(n.risk = PD_pois_6$n.risk+ PD_pois_6$n.event, 
                      num_event = PD_pois_6$n.event, 
                      num_censor = PD_pois_6$n.censor,
                      time = PD_pois_6$time)


bayes.pois <- MCMCpoisson(num_event ~ log(n.risk), b0 = c(hist_mean, 1), B0 = c(1/(hist_dev)**2, 1e8),data = PD_pois[1:4], burnin = 10000, mcmc = 1e6, thin = 10)

project_time_frame <- 7:25
hazards <- sample(bayes.pois[,1],100) %>% exp()
preds <- dropoutMCMC(exposure =dat_6$TRTDUR[dat_6$EOT == 0], duration = length(project_time_frame),enroll_pred = c(1, rep(0,17)), hazards = hazards,num_sim = 100)
drop_trace <- do.call(cbind, preds$at.risk) %>% 
  apply(1, quantile,probs = c(0.025, 0.50, 0.975)) %>% 
  round

ylow <- c(num_at_risk[7], drop_trace[1,])
ymid <- c(num_at_risk[7], drop_trace[2,])
yhigh <- c(num_at_risk[7], drop_trace[3,])
yobserved <- num_at_risk[1:13]
ggplot() + 
  geom_line(aes(x = 1:7, y = num_at_risk[1:7]),size = 1.3) +
  geom_line(aes(x = project_time_frame, y = ylow), linetype = "dotted", size = 1) +
  geom_line(aes(x = project_time_frame, y = ymid), linetype = "dashed", size = 1) +
  geom_line(aes(x = project_time_frame, y = yhigh), linetype = "dotted", size = 1) +
  geom_line(aes(x = 1:13, y = yobserved), size = 1, color = "red") +
  xlab("Months since First patient dosing") +
  ylab("Number of patients on treatment,95% CI") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_y_continuous(breaks = seq(0,70, 10)) +
  theme_igray()

## using Non Bayes model

project_time_frame <- 7:24
hazards <- exp(-4.09)
preds <- dropoutMCMC(exposure =dat_6$TRTDUR[dat_6$EOT == 0], duration = length(project_time_frame)-1,enroll_pred = c(1, rep(0,17)), hazards = hazards,num_sim = 1000)
drop_trace <- do.call(cbind, preds$at.risk) %>% 
  apply(1, quantile,probs = c(0.025, 0.50, 0.975)) %>% 
  round

ylow <- c(num_at_risk[7], drop_trace[1,])
ymid <- c(num_at_risk[7], drop_trace[2,])
yhigh <- c(num_at_risk[7], drop_trace[3,])
yobserved <- num_at_risk[1:13]
ggplot() + 
  geom_line(aes(x = 1:7, y = num_at_risk[1:7]),size = 1.3) +
  geom_line(aes(x = project_time_frame, y = ylow), linetype = "dotted", size = 1) +
  geom_line(aes(x = project_time_frame, y = ymid), linetype = "dashed", size = 1) +
  geom_line(aes(x = project_time_frame, y = yhigh), linetype = "dotted", size = 1) +
  geom_line(aes(x = 1:13, y = yobserved), size = 1, color = "red") +
  xlab("Months since First patient dosing") +
  ylab("Number of patients on treatment,95% CI") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_y_continuous(breaks = seq(0,70, 10)) +
  theme_igray()

## output file

observed <- rbind(rep(0,13),num_at_risk,rep(0,13))
pred_mat <- cbind(observed,drop_trace)
pred_mat <- cbind(1:25, t(pred_mat))
colnames(pred_mat) <- c("Month","2.5%","50%","97.5%")
write.xlsx(pred_mat,"results/predicted number of patients on treatment.xlsx")

## change point model
a <- MCMCpoissonChange(num_event ~ log(n.risk), m=1, b0 = c(hist_mean, 1), B0 = c(1/(hist_dev)**2, 1000000),data = PD_pois) 
b <- MCMCpoissonChange(num_event ~ log(n.risk), m=2, b0 = c(hist_mean, 1), B0 = c(1/(hist_dev)**2, 1000000),data = PD_pois) 

## Non
