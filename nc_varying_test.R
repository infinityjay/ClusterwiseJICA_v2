library("mclust")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/CJICA_asist_functions.R")

# generate simulation data with vary Q across clusters
test1 <- Simulate_CJICA(Nk = 10, Vm = 2500, K = 4, Qm = c(5, 5, 5, 5), E = .01, M = 2, type = 1, cor = .05)
str(test1)

result <- ClusterwiseJICA(X = test1$Xe, k = 4, starts = 100)


adjustedRandIndex(test1$P, result[[10]]$p)


debugSource("./functions/data_simulation_functions.R")
set.seed(123)
test_data <- Simulate_CJICA(Nk = 10, Vm = 2500, K = 4, Qm = c(5, 5, 5, 5), E = .1, M = 2, type = 1, cor = .05)
