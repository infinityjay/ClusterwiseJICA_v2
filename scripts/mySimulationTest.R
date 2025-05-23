## # load R module
## module load R/4.4.0-gfbf-2023a
## # run script
## Rscript my_script.R 1

setwd('/home/s4162315/cjica/ClusterwiseJICA_v2/')

library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")

#### generate simulation data

Nk <- c(50, 100)
# Qvect <- list(c(3,3), c(8,8), c(2,8), c(2,2,2,2), c(2,8,3,8))
Qvect_list <- list(c(3,3), c(8,8))
E <- c(0.2, 0.6)
cor <- c(0, .60, .94)
VAF <- c(0.8, 1)
rep <- 1:5

base_grid <- expand.grid(Nk=Nk, Qid=seq_along(Qvect_list), E=E, cor=cor, VAF=VAF, rep=rep)
base_grid$Qvect <- Qvect_list[base_grid$Qid]

# Remove Qid if not needed
grid <- base_grid[, !(names(base_grid) %in% "Qid")]
args <- commandArgs(TRUE)
#args <- as.numeric(args)

splits <- split(1:nrow(grid), ceiling(seq_along(1:nrow(grid))/10))
sp <- args[1]

rows <- splits[[sp]]


#### run simulation and c-jica function

Vm <- 2500
complex <- c(1,3,4,6,7)
grid <- grid[rows,]
for(sim in 1:nrow(grid)){
  seed <- 123
  set.seed(seed)
  simdata <- Simulate_CJICA(
    Nk = grid[sim, ]$Nk, 
    Vm = Vm,
    K = length(grid$Qvect[[sim]]),
    Qvect = grid$Qvect[[sim]],
    E = grid[sim, ]$E,
    M = 2,
    type = 1,
    cor = grid[sim, ]$cor,
    VAF = grid[sim, ]$VAF
  )

  for(comp in complex) {
    # run cjica and calculate executing time
    ptm <- proc.time()
    cjica <- ClusterwiseJICA_varyQ(
      X = simdata$Xe,
      k = length(grid$Qvect[[sim]]), 
      useInputQ = F, 
      starts = 100, 
      scale = F, 
      VAF = grid[sim, ]$VAF, 
      complex = comp, 
      useChull = T, 
      Vm = Vm
    )
    time <- proc.time() - ptm
    output <- list()
    output$time <- time
  }
  ext <- './result/Sim1/'  
  ext <- paste(ext,'CJICA_sim1_',rows[sim], '.Rdata',sep = '')
  save(output,file = ext)
}


#### analyze


