
setwd('/home/s4162315/cjica/ClusterwiseJICA_v2/')

library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")


set.seed(123)
data_Nk50_Vm500_K2_Q33_E04_VAF1 <- Simulate_CJICA(Nk = 50, Vm = 500, K = 2, Qvect = c(3,3), E = .4, M = 2, type = 1, cor = .05, VAF = 1)
save(data_Nk50_Vm500_K2_Q33_E04_VAF1, file = "./data/data_Nk50_Vm500_K2_Q33_E04_VAF1.RData")

