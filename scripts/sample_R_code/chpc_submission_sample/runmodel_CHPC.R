argos <- commandArgs(trailingOnly = TRUE)
min <- as.numeric(argos[1])
max <- as.numeric(argos[2])
spl <- as.numeric(argos[3])
splseq <- seq(from = min, to = max - spl + 1, length.out = (max - min) / spl)

library(deSolve)
library(iterators)
library(Rmpi)
library(doMPI)
library(foreach)
library(doParallel)

# Hello

get(load(file = "TB_model.RData"))
get(load(file = "artgen_general.RData"))
get(load(file = "nrowparammat.RData"))
get(load(file = "artinit_tb_func.RData"))
get(load(file = "timevarpar.RData"))
get(load(file = "hivinc_func.RData"))
get(load(file = "prob_hiv_func.RData"))
get(load(file = "timesvect.RData"))

get(load(file = "fun_diagstp_notutt.RData"))
get(load(file = "fun_diagsfp_notutt.RData"))


#### creating result list
resultList <- list()

# nrestry= as.numeric(args6[2])

# Functions to calculate rates#
##############################

accdiag_func <- function(times, acc_tutt) {
  roundtime <- floor(times)
  if (roundtime < 2022) {
    acc <- 0
  } else if (roundtime >= 2022) {
    acc <- acc_tutt
  }
  return(as.numeric(acc))
}

### Function for probability true-positive diagnoses
#
# fun_diagstp_notutt = function(relfoot_phc,ltfu_prediag,sputsc,sens_diag,clindiag, ltfu_diag, iltfu_phc, iltfu_other) {
#   pdiagstp_notutt =
#     relfoot_phc*(1-ltfu_prediag) * (1-sputsc) * sens_diag *(1-ltfu_diag) * (1-iltfu_phc) +  # bacteriologically-confirmed in PHC
#     relfoot_phc*(1-ltfu_prediag) * (1-sputsc) * (1-sens_diag)* clindiag *(1-ltfu_diag) * (1-iltfu_phc) + # bacteriologically-negative, clinically diagnosed in PHC
#     relfoot_phc*(1-ltfu_prediag) * (sputsc) * clindiag *(1-ltfu_diag) * (1-iltfu_phc) + # sputum scarce, clinically diagnosed in PHC
#     (1-relfoot_phc)*(1-ltfu_prediag) * (1-sputsc) * sens_diag * (1-ltfu_diag) * (1-iltfu_other) + # bacteriologically-confirmed in Other
#     (1-relfoot_phc)*(1-ltfu_prediag) * (1-sputsc) * (1-sens_diag) * clindiag * (1-ltfu_diag) * (1-iltfu_other) +  # # bacteriologically-negative, clinically diagnosed in Other
#     (1-relfoot_phc)*(1-ltfu_prediag) * (sputsc) * clindiag * (1-ltfu_diag) * (1-iltfu_other)  # sputum scarce, clinically diagnosed in Other
#   return(pdiagstp_notutt)
# }
#
# ### Function for for probability false-positive diagnoses
#
# fun_diagsfp_notutt = function(relfoot_phc,ltfu_prediag,prediag_tp_h0,sputsc,spec_diag,clindiag, ltfu_diag, iltfu_phc, iltfu_other) {
#   pdiagsfp_notutt =
#     relfoot_phc*(1-ltfu_prediag) * (1-sputsc) * (1-spec_diag) *(1-ltfu_diag) * (1-iltfu_phc) +  # bacteriologically false-positive in PHC
#     (1-relfoot_phc)*(1-ltfu_prediag) * (1-sputsc) * (1-spec_diag)  * (1-ltfu_diag) * (1-iltfu_other)  # bacteriologically false-positive in Other
#   return(pdiagsfp_notutt)
# }


funcMakeResults <- function() {



  #####################
  ### RUNNING THE MODEL

  results <- as.data.frame(ode(
    func = TB_model,
    times = timesvect,
    y = yinit.r[[i]], # .r file created in loop below
    parms = ParamList2.r[[i]], # .r file created in loop below
    method = "rk", hini = 0.05
  )) # lsoda

  # results = resultList[[1]]

  ### Population size (N)
  results$N <- results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 + results$RH_h0 + results$R_h0 +
    results$S_h1 + results$LR_h1 + results$LD_h1 + results$IP_h1 + results$IS_h1 + results$DP_h1 + results$DS_h1 + results$FN_h1 + results$T_h1 + results$RH_h1 + results$R_h1 +
    results$S_h2 + results$LR_h2 + results$LD_h2 + results$IP_h2 + results$IS_h2 + results$DP_h2 + results$DS_h2 + results$FN_h2 + results$T_h2 + results$RH_h2 + results$R_h2 +
    results$S_h3 + results$LR_h3 + results$LD_h3 + results$IP_h3 + results$IS_h3 + results$DP_h3 + results$DS_h3 + results$FN_h3 + results$T_h3 + results$RH_h3 + results$R_h3 +
    results$S_h4 + results$LR_h4 + results$LD_h4 + results$IP_h4 + results$IS_h4 + results$DP_h4 + results$DS_h4 + results$FN_h4 + results$T_h4 + results$RH_h4 + results$R_h4 +
    results$S_h1a + results$LR_h1a + results$LD_h1a + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a + results$T_h1a + results$RH_h1a + results$R_h1a +
    results$S_h2a + results$LR_h2a + results$LD_h2a + results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$T_h2a + results$RH_h2a + results$R_h2a +
    results$S_h3a + results$LR_h3a + results$LD_h3a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$T_h3a + results$RH_h3a + results$R_h3a +
    results$S_h4a + results$LR_h4a + results$LD_h4a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a + results$T_h4a + results$RH_h4a + results$R_h4a

  ### Population size (N_h0)
  results$N_h0 <- results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 + results$RH_h0 + results$R_h0

  ### Population size (N_hx)
  results$N_hx <- results$N - results$N_h0

  ## Notified TB cases
  results$notifTB <- numeric(length = nrow(results))
  for (k in results$time) {
    results$notifTB[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, results$CumTBnotif[which(results$time == k + 1)] - results$CumTBnotif[which(results$time == k)])
  }

  ## False positives
  results$FP <- numeric(length = nrow(results))
  for (k in results$time) {
    results$FP[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, results$CumFP[which(results$time == k + 1)] - results$CumFP[which(results$time == k)])
  }

  ## Percentage false-positives
  results$p_FP <- results$FP / results$notifTB

  ## Inc TB cases
  results$IncTB <- numeric(length = nrow(results))
  for (k in results$time) {
    results$IncTB[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, results$CumTBinc[which(results$time == k + 1)] - results$CumTBinc[which(results$time == k)])
  }

  ## Inc TB cases_h0
  results$IncTB_h0 <- numeric(length = nrow(results))
  for (k in results$time) {
    results$IncTB_h0[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, results$CumTBinc_h0[which(results$time == k + 1)] - results$CumTBinc_h0[which(results$time == k)])
  }

  ## Inc TB cases_hx
  results$IncTB_hx <- results$IncTB - results$IncTB_h0

  ## Percentage Inc TBHIV
  results$pIncHIV <- results$IncTB_hx / (results$IncTB_hx + results$IncTB_h0)

  ## TB deaths
  results$TBdeaths <- numeric(length = nrow(results))
  for (k in results$time) {
    results$TBdeaths[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, results$CumTBmort[which(results$time == k + 1)] - results$CumTBmort[which(results$time == k)])
  }

  results$TBinc100T <- results$IncTB / results$N * 100000
  results$TBmort100T <- results$TBdeaths / results$N * 100000
  results$TBinc_h0_100T <- results$IncTB_h0 / results$N_h0 * 100000
  results$TBinc_hx_100T <- results$IncTB_hx / results$N_hx * 100000
  results$hivTBinc_100T <- results$IncTB_hx / results$N * 100000

  ### Size of people in diseased compartments (Prevalence)
  results$TB <- results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$IP_h1 + results$IS_h1 + results$DP_h1 + results$DS_h1 + results$FN_h1 + results$IP_h2 + results$IS_h2 + results$DP_h2 + results$DS_h2 + results$FN_h2 +
    results$IP_h3 + results$IS_h3 + results$DP_h3 + results$DS_h3 + results$FN_h3 + results$IP_h4 + results$IS_h4 + results$DP_h4 + results$DS_h4 + results$FN_h4 + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a +
    results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a

  ### Prevalence of TB disease
  results$p_TB <- results$TB / results$N

  ### Proportion prevalence that are pre-symptomatic
  results$p_presymp <- (results$IP_h0 + results$DP_h0 + results$IP_h1 + results$DP_h1 + results$IP_h2 + results$DP_h2 + results$IP_h3 + results$DP_h3 + results$IP_h4 + results$DP_h4 + results$IP_h1a + results$DP_h1a + results$IP_h2a + results$DP_h2a +
    results$IP_h3a + results$DP_h3a + results$ IP_h4a + results$DP_h4a) / results$TB

  ### Proportion symptomatic that previously sought care
  results$p_prevcare <- (results$FN_h0 + results$FN_h1 + results$FN_h2 + results$FN_h3 + results$FN_h4 + results$FN_h1a + results$FN_h2a + results$FN_h3a + results$FN_h4a) /
    (results$IS_h0 + results$DS_h0 + results$FN_h0 + results$IS_h1 + results$DS_h1 + results$FN_h1 + results$IS_h2 + results$DS_h2 + results$FN_h2 + results$IS_h3 + results$DS_h3 +
      results$FN_h3 + results$IS_h4 + results$DS_h4 + results$FN_h4 + results$IS_h1a + results$DS_h1a + results$FN_h1a + results$IS_h2a + results$DS_h2a + results$FN_h2a +
      results$IS_h3a + results$DS_h3a + results$FN_h3a + results$IS_h4a + results$DS_h4a + results$FN_h4a)


  ### HIV prevalence and ART coverage (all ages)
  results$HIV <- results$N - (results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 + results$R_h0)

  results$ART <- results$S_h1a + results$LR_h1a + results$LD_h1a + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a + results$T_h1a + results$R_h1a +
    results$S_h2a + results$LR_h2a + results$LD_h2a + results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$T_h2a + results$R_h2a +
    results$S_h3a + results$LR_h3a + results$LD_h3a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$T_h3a + results$R_h3a +
    results$S_h4a + results$LR_h4a + results$LD_h4a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a + results$T_h4a + results$R_h4a

  results$p_HIV <- results$HIV / results$N
  results$p_ART <- results$ART / results$HIV

  ## Mid-year population
  results$Nmid <- numeric(length = nrow(results))
  for (k in results$time) {
    results$Nmid[which(results$time == k)] <-
      ifelse(k == tail(results$time, n = 1), NA, (results$N[which(results$time == k)] + results$N[which(results$time == k + 1)]) / 2)
  }

  # results$I_all = results$IS_h0 + results$IS_h1 + results$IS_h2 + results$IS_h3 + results$IS_h4 +
  #  results$IS_h1a + results$IS_h2a + results$IS_h3a + results$IS_h4a +
  #  results$IP_h0 + results$IP_h1 + results$IP_h2 + results$IP_h3 + results$IP_h4 +
  #  results$IP_h1a + results$IP_h2a + results$IP_h3a + results$IP_h4a

  ### Creating dataset for calibration

  results <- results[, -c(2:105)] # 94
  # results = subset(results, select=-c(N,TBdeaths,TB,HIV,ART))

  # print(results)
  return(results)
  # return(modelcal)
}

for (a in splseq) {
  yinit.r <- get(load(file = "yinit.RData"))[a:(a - 1 + spl)]
  ParamList2.r <- get(load(file = "ParamList2.RData"))[a:(a - 1 + spl)]

  cl <- startMPIcluster()
  registerDoMPI(cl)
  # size=clusterSize

  finalMatrix <- foreach(
    i = 1:spl,
    .packages = c(
      "deSolve", "foreach", "Rmpi", "iterators",
      "doMPI", "doParallel"
    )
  ) %dopar% {
    tempMatrix <- funcMakeResults()

    tempMatrix # Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }

  saveRDS(finalMatrix, file = paste0("resultList_CHPC", a, ".RData"))

  # closeCluster(cl)
}

resultList <- vector(mode = "list")
for (a in splseq) {
  resultList <- c(resultList, readRDS(file = paste0("resultList_CHPC", a, ".RData")))
}

saveRDS(resultList, file = "resultList_CHPC.RData")
