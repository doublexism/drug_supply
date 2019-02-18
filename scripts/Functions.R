# function definitions
library(purrr)
library(ggplot2)
library(sas7bdat)
library(MASS)
library(lubridate)
library(survival)
library(dplyr)
library(coxme)
library(survminer)
library(readxl)
library(MCMCpack)
library(stringr)
source(paste0(Sys.getenv("OneDrive"),"\\Codes\\Power4CT\\Power4CT\\functions.R"))

##  time-on-treatment at some time point for one patient
timeOnTrt <- function(TE, TD, TC, delta,u){
  # TE: the time of enrollment
  # TD: the time of Drop-out (due to either event or censoring)
  # TC: the time of cross-over if assigned to control group
  # Delta: the indicator of treatment
  # u: the time point
  time_on_treatment <- ((min(TD, (u-TE)))*delta + (min(TD,(u-TC)))*(1-delta)*(TD > TC))*(TE < u)
  return(time_on_treatment)
}

## Total time-on-treatment at some time point of the whole experiment
totalTime <- function(n, enroll_time, p_trt, gamma, lambda,rho,lambda_pd,rho_pd,u){
  # n: number of patients to recruite
  # enroll_time: length of patient enroll time
  # p_trt: probability of assigning to treatment
  # gamma: folds of prolonged PFS in treatment group
  # lambda and rho: parameter of weibull distribution
  # lamdba_pd and rho_pd: paramter of weibull distribution about pd in control group
  # u: the time point
  require(foreach)
    TE_vec <- cumsum(rexp(n, rate = n/enroll_time))
    vec_time <-   foreach(i = 1:n, .combine = c) %do% {
    # generate random enroll time
      TE <- TE_vec[i]
      delta <- rbinom(1,1,p_trt)
      TD <- rweibull(1,rho,lambda)*(gamma*delta+(1-delta))
      TC <- rweibull(1,rho_pd,lambda_pd)
      # print(c(TC,TD,TE,delta))
      time_on_treatment <- timeOnTrt(TE, TD, TC, delta, u)
      return(time_on_treatment)
    } 
#    print(vec_time)
  total_time <- sum(vec_time)
  return(total_time)
}


# weibull distribution function
weibull_f <- function(t, lambda, p, type = "f"){
  # t: time
  # lambda: scale parameter
  # p: shape parameter
  # type: type of functions
  h <- lambda**p*p*t**(p-1)
  s <- exp(-(lambda*t)**p)
  f <- h*s
  if (type == "f"){
    return(f)} else if (type == "s"){
      return(s)
    } else if (type == "h"){
      return(h)
    }
}
# fit and predict
dist_fit <- function(x,data,dist, subset.var = NULL, subset.val = NULL){
  # data: dataset
  # dist: distribution
  # subset.var: variable to subset
  # subset.val: value to subset
  require(ggplot2)

  if (!is.null(subset.var)){
    data <- data[data[[subset.var]]==subset.val,]
  }
  day <- data[[x]] - min(data[[x]])+1
  lineplot.dat <- data.frame(time = day, observed = seq(0,1, length.out = nrow(data)+1)[-1])
  # calculate predicted value
  for (d in dist){
    if (d == "exponential") {
      fdist = "exp"
    } else {
      fdist = d
    }
    estimate <- MASS::fitdistr(data[[x]], d)$estimate
    s <- do.call(paste0("p",fdist),c(list(data[[x]]), as.list(estimate)))
    lineplot.dat[[d]] <- s
  }
  return(lineplot.dat)
}

## pad missing time series
tsPad <- function(data, by, dtime,maxby = NULL, keep = NULL, unify_stdt = FALSE){
  if (is.null(by)){
    maxdt <- max(data[[dtime]])
    mindt <- min(data[[dtime]])
    rep_keep <- c(diff(data[[dtime]]),1)
    keepD <- map_dfc(data[keep], rep, rep_keep)
    resultD <- data.frame(time = seq(mindt, maxdt), time_spent = seq(mindt, maxdt)) %>% 
      bind_cols(keepD) 
    return(resultD)
  } else {
    by_num <- length(by)
    byvar <-  map(by, ~`[[`(data,.))
    maxbyvar <- map(maxby, ~`[[`(data,.))
    maxdt <- by(data, byvar, function(x) max(x[[dtime]]), simplify = TRUE) %>% as.numeric() %>% na.omit()
    if (unify_stdt == FALSE){
      mindt <- by(data, byvar, function(x) min(x[[dtime]]), simplify = TRUE) %>% as.numeric() %>% na.omit()
    } else
    {
      mindt <- rep(min(data[[dtime]]),length(maxdt))
    }
    time_num <- maxdt-mindt+1
    if (!is.null(maxby)){
      maxdt <- by(data, maxbyvar,  function(x) max(x[[dtime]]), simplify = TRUE) %>% as.numeric() %>% na.omit()
      maxbysum <- rowSums(table(data[[maxby]], data[[setdiff(by, maxby)]]) > 0 )
      maxdt <- rep(maxdt, maxbysum)
    } else {   maxdt <- rep(max(maxdt),length(time_num))
    }
    rep_num <- maxdt-mindt+1
    time <- unlist(map2(mindt, maxdt, seq))
    byvarD <- distinct(data[by])
    byvar.list <- map(by, ~`[[`(byvarD,.))
    byvarD <- byvarD[do.call(order,rev(byvar.list)),] %>% na.omit()
    resultD <- map_dfc(byvarD, rep, rep_num)
    resultD$time <- time
    resultD <- distinct(left_join(resultD, data[c(dtime, keep)], by = c("time" = dtime)))
    return(resultD)
  }
}


dropoutMCMC <- function(exposure,enroll_pred = NULL, duration = 12,hazards = NULL, hazard_func = NULL, num_sim=1000, num_int = 21, num_pred_sample = 100,num_enroll_sample=100,exact = FALSE){
  
  if (is.null(enroll_pred)){
    enroll_pred <- matrix(rep(0, duration),nrow = 1)
    message("Assume enrollment has complete...")
    num_pred_sample <- 1
    } else if (is.atomic(enroll_pred)) {
    enroll_pred <- matrix(enroll_pred,nrow = 1)  
    } else { 
    enroll_pred <- t(enroll_pred[,sample(ncol(enroll_pred), num_pred_sample)])
    }
  if (!is.null(hazards)){
    hazard_func <- function(t,k){hazards[k]}
    num_pred_sample <- length(hazards)
  }
  as.risk.list <- list()
  drop.list <- list()
  for (k in 1:num_pred_sample){
    if (num_pred_sample <- 1){
      enroll <- 1
    } else {
      enroll <- k
    }
    for (i in 1:num_sim){
      at.risk <- exposure
      num_events <- c()
      num_at_risk <- c()
      for (j in 1:duration){
        #     print(at.risk)
        #     print(num_events)
        events <- as.numeric(hazard_func(at.risk, k)>runif(length(at.risk)))
        num_events[j] <- sum(events)
        at.risk <- c(at.risk[events == 0]+1,rep(1,round(enroll_pred[enroll,j])))
        num_at_risk[j] <- length(at.risk)
      } 
#      print(at.risk)
      drop.list[[k*num_sim-num_sim+i]] <- num_events
      as.risk.list[[k*num_sim-num_sim+i]] <- num_at_risk
    }
  }  
  return(list(drop = drop.list,at.risk = as.risk.list))
}

## get posterior sampling data from mcmc object
postSample <- function(x, fixed = list(), random = list(), group = list()){
  fixedI <- str_which(colnames(x), "^beta")
  random_int <- str_which(colnames(x), "^b\\.\\(Intercept\\)")
  random_slope <- str_which(colnames(x), "^b\\.[^\\(]")
  vcv <- str_which(colnames(x), "VCV")
  sigma <- str_which(colnames(x),"sigma2")
  varI <- c(1:(ncol(x)-1))[-c(vcv)]
  fixed_level <- map(fixed, `[`, -1)
  int_level <- map(fixed, `[`,1)
  var_level <- c(paste(int_level, collapse = ":"), unlist(fixed_level), kronecker(unlist(random), unlist(group), paste, sep = ":"),"sigma")
  if (!is_empty(var_level)){
    var <- var_level
  } else {
    var <- colnames(x)[varI]
  }
  beta <- x[,fixedI]
  beta0_i <- x[,random_int]
  beta_i <- x[,random_slope]
  sigma <- rnorm(nrow(x),0, sqrt(x[,sigma]))
  dat <- as.data.frame(cbind(beta, beta0_i, beta_i, sigma))
  colnames(dat) <- var
  return(dat)
}
## encoder
encoder <- function(xlevel,y, ylevel, group, ngroup){
  group_code <- c(rep(0, ngroup+ngroup*sum(ylevel - length(y))))
    if (group > ngroup){
      group_code[((y-1)*ngroup+1):(y*ngroup)] <- 1/ngroup
    } else {
      group_code[((y-1)*ngroup+group)] <- 1
    }
  coding <- rep(0, sum(xlevel) - length(xlevel))
  int <- 1
  return(c(int, coding, group_code,1))
}

## prediction
postPred <- function(sample, random, group, nlevel.fix, nlevel.random, ngroup){
  code <- encoder(nlevel.fix, random, nlevel.random, group, ngroup)
  lambda <- as.matrix(sample) %*% code
  return(lambda)
}


## plot posterior distribution
forest_plot <- function(data, text_col, data_col, title_col = 5,...){
  text <- rbind(colnames(data[text_col]),as.matrix(data[text_col]))
  dat <- bind_rows(data.frame(Mean = NA,Lower = NA,Upper = NA), data[data_col])
#  title <- paste0("The Forest Plot of ", data[[title_col]][1], " From Published Studies")
  summary <- c(TRUE, rep(FALSE,nrow(dat)))
  forestplot <- t(text, dat, title = title,is.summary = summary, ...)
}


## Reorganized observed data
dataCutOff <- function(data, cutoff, PFS = "TRTDUR", Event = "EOT", timeInTrial = "timeInTrial", 
                       timeOutTrial = "timeOutTrial", enrolled = "enrolled",censor_rate = 0.01, start_hazard = 1){
  ## time to event data set generation
  dat_c <- data[data[[timeInTrial]] < 7,]
  dat_c[[Event]][dat_c[[PFS]] + dat_c[[timeInTrial]] > 7] <- 0
  dat_c[[PFS]][dat_c[[PFS]] > 7] <- 7
  formula1 <- as.formula(sprintf("Surv(%s, %s)~1", PFS, Event))
  PD_pois <- survfit(formula1, data = dat_c) %>% summary(seq(start_hazard, start_hazard + cutoff, 1))
  PD_pois <- data.frame(n.risk = PD_pois$n.risk+ PD_pois$n.event, 
                          num_event = PD_pois$n.event, 
                          num_censor = PD_pois$n.censor,
                          time = PD_pois$time)
  formula2 <- as.formula(sprintf("Surv(%s, %s)~1", timeInTrial, enrolled))
  formula3 <- as.formula(sprintf("Surv(%s, %s)~1", timeOutTrial, Event))
  # length of enrollment and number of risk
  l_enroll <- floor(max(data[[timeInTrial]]))
  l_time <- floor(max(data[[timeOutTrial]]))
  tot_enroll <- sum(data[[enrolled]])
  tot_event <- sum(data[[Event]])
  n.enroll <- survfit(formula2, data = data) %>% summary(1:l_enroll) %>% `[[`("n.event") %>% cumsum()
  if (l_time > l_enroll){
    n.enroll <- c(n.enroll, rep(tot_enroll, l_time - l_enroll))
  }

  n.events <- survfit(formula3, data = data) %>% summary(1:l_time) %>% `[[`("n.event") %>% cumsum()
  n.events <- switch(as.character(length(n.events)-l_time), 
                     "0" = n.events, 
                     "1" = c(n.events, tot_event - max(n.events)),
                     c(n.enroll, tot_event - max(n.events),rep(0,length(n.events)-l_time-1)))
  num_at_risk <- n.enroll - n.events
  return(list(life_tab = PD_pois, data = dat_c, num_at_risk = num_at_risk))
}

## MCMC Ploting
MCMCplot <- function(MCMCobj,time_frame_to_project = NULL, observed = NULL,lower = 0.025, upper = 0.975, x_breaks = NULL, y_breaks = NULL){
  ## plot MCMC prediction and Intervals
  #MCMCobj: output from dropMCMC
  #observed: observed number at risks vector
  trace <- do.call(cbind, MCMCobj$at.risk) %>% 
    apply(1, quantile,probs = c(lower, 0.50, upper)) %>% 
    round
  if (is.null(time_frame_to_project)){
    time_frame_to_project <- (1:length(MCMCobj$drop[[1]]))+length(observed)
  }

  if (is.null(x_breaks)){
    x_breaks <- ceiling(max(time_frame_to_project)/20)
  } 
  if (is.null(y_breaks)) {
    y_breaks <- ceiling(max(max(trace),observed)/20)
  }
  if (is.null(observed)){
    ylow <- trace[1,]
    ymid <- trace[2,]
    yhigh <- trace[3,]
    observed_line <- NULL
  } else {
    time_frame_to_project <- c(min(time_frame_to_project)-1, time_frame_to_project)
    ylow <- c(tail(observed,1), trace[1,])
    ymid <- c(tail(observed,1), trace[2,])
    yhigh <- c(tail(observed,1), trace[3,])
    observed_line <- geom_line(aes(x = 1:length(observed), y = observed), size = 1, color = "red") 
  }
  print(time_frame_to_project)
  plot <- ggplot() + 
    geom_line(aes(x = time_frame_to_project, y = ylow), linetype = "dotted", size = 1) +
    geom_line(aes(x = time_frame_to_project, y = ymid), linetype = "dashed", size = 1) +
    geom_line(aes(x = time_frame_to_project, y = yhigh), linetype = "dotted", size = 1) +
    observed_line +
    xlab("Months since First patient dosing") +
    ylab("Number of patients on treatment,95% CI") +
    theme(legend.position = "right") +
    scale_x_continuous(breaks = seq(1,max(time_frame_to_project)+x_breaks,x_breaks)) +
    scale_y_continuous(breaks = seq(0,max(max(trace),observed)+y_breaks, y_breaks)) +
    theme_igray()
  return(list(plot = plot, trace = trace))
}

## transform historic info into prior

priorSpecs <- function(class = "median PFS historic", median=NULL, sample_size=NULL, PFS_rate=NULL, time=NULL, unit=NULL,  best=NULL, worst=NULL, lower=NULL, upper=NULL){
  # class: "median PFS historic", "PFS rate historic", "design", "early phase","expert opinion"
  
}


