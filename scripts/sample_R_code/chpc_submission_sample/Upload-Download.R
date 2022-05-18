library(ssh)

###############################################################
###### UPLOAD ###############

session <- ssh_connect("fmarx@lengau.chpc.ac.za")

save(yinit, file = "./CHPC upload/yinit.RData")
save(ParamList2, file = "./CHPC upload/ParamList2.RData")
save(TB_model, file = "./CHPC upload/TB_model.RData")
save(fun_diagstp_notutt, file = "./CHPC upload/fun_diagstp_notutt.RData")
save(fun_diagsfp_notutt, file = "./CHPC upload/fun_diagsfp_notutt.RData")

scp_upload(session, files = "./CHPC upload/yinit.RData", to = "lustre/satbmodel", verbose = TRUE)
scp_upload(session, files = "./CHPC upload/ParamList2.RData", to = "lustre/satbmodel", verbose = TRUE)
scp_upload(session, files = "./CHPC upload/TB_model.RData", to = "lustre/satbmodel", verbose = TRUE)
scp_upload(session, files = "./CHPC upload/fun_diagstp_notutt.RData", to = "lustre/satbmodel", verbose = TRUE)
scp_upload(session, files = "./CHPC upload/fun_diagsfp_notutt.RData", to = "lustre/satbmodel", verbose = TRUE)

ssh_exec_wait(session, command = "qsub runsatb.qsub")
# ssh_exec_wait(session, command = 'qsub getall.qsub')
# ssh_exec_wait(session, command = 'qstat -awu fmarx')
ssh_disconnect(session)


###### DOWNLOAD ###############

nrowparammat <- 100000 # avoid issues with differing numbers of iterations

session <- ssh_connect("fmarx@lengau.chpc.ac.za")
# ssh_session_info(session)
ssh_exec_wait(session, command = "ls -d lustre/satbmodel/resultList_CHPC*")
# ssh_exec_wait(session, command = 'qstat -awu fmarx')

today <- paste0(year(Sys.Date()), month(Sys.Date()), day(Sys.Date()))
downpath <- c(paste0("~/Documents/Research/Gates TB Diagnostics Modelling/SA TB Diag Model 2021/CHPC downloads/", today, "/"))
dir.create(downpath)
scp_download(session, files = "lustre/satbmodel/resultList_CHPC.RData", to = downpath, verbose = TRUE)

# Terminal:
# scp fmarx@scp.chpc.ac.za:/mnt/lustre/users/fmarx/satbmodel/resultList_CHPC.RData /Users/florian/Downloads/resultList_CHPC.RData

ssh_disconnect(session)

resultList <- readRDS(file = c(paste0(downpath, "/resultList_CHPC.RData")))
