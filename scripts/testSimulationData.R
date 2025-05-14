
setwd('/home/s4162315/cjica/ClusterwiseJICA_v2/')

library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")


set.seed(123)
data_Nk50_K5_Q44566_E04_VAF08 <- Simulate_CJICA(Nk = 50, Vm = 2500, K = 5, Qvect = c(4,4,5,6,6), E = .4, M = 2, type = 1, cor = .05, VAF = 0.8)
save(data_Nk50_K5_Q44566_E04_VAF08, file = "./data/data_Nk50_K5_Q44566_E04_VAF08.RData")