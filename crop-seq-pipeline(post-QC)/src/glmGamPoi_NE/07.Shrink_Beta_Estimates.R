######## Load all libraries 


## Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(argparse)
library(BiocParallel)
library(parallel)

##########################################################################################################
#############################        Part 1: Load Arguments        #######################################
##########################################################################################################

## Load arguments
parser <- ArgumentParser()
parser$add_argument('-s','--sublibrary', required =TRUE,
                    help='name of sublibrary to analyze')
parser$add_argument('-niter','--niter', required =TRUE,
                    help='max number of shrinking iterations to run')

## Print loaded arguments
args <- parser$parse_args()
print(args)



##########################################################################################################
########################        Part 2: Specify Input and Output Files        ############################
##########################################################################################################

## Base Directory
base_dir = "/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev_NE/results/"
dir =  paste(base_dir, "unweightedglmGamPoi_gembatch_results_downstream_concat",sep = "")



## File names
estimate_file = file.path(dir, paste0("beta_estimates_all_targets.csv.gz"))
se_file = file.path(dir, paste0("beta_standarderrors_all_targets.csv.gz"))

## Number of times iterate
niters = as.numeric(args$niter)

## File out name
out_file  = file.path(dir, paste0("beta_estimates_shrunken_all_targets_niter",niters,".rds"))
out_values_file  = file.path(dir, paste0("beta_estimates_shrunken_all_targets_niter",niters,"_values.csv"))

##########################################################################################################
#############################        Part 3: Define Functions        #######################################
##########################################################################################################

## Function to convert loaded data into matrix of proper format
DatToMat <- function(dat){
  # data frame where first column is colnames, convert to matrix
  mat <- as.matrix(dat[, -1, with = FALSE])
  rownames(mat) <- dat[[1]]
  return(mat)
}


## Function to convert loaded data into matrix of proper format
alternating_rank1stein_shrinkage <- function(lfcs, se, BPPARAM, maxit=8,
                                             a_min=1e-8, a_max=100, verbose=TRUE) {
  # Stolen from Jack Kamm 
  # https://code.roche.com/kammj2/ngs5388_trem2feps_pilot9/-/blob/main/src/alternating_rank1stein.R?ref_type=heads 
  a_loss_global <- function(a) {
    -sum(dnorm(lfcs, mean=0, sd=sqrt(a + se^2),
               log=TRUE))
  }
  
  if (verbose) print("Initial global estimate of A")
  
  a_opt_global <- optimize(a_loss_global, c(a_min, a_max))
  
  curr_a1 <- rep(sqrt(a_opt_global$minimum), nrow(lfcs))
  curr_a2 <- rep(sqrt(a_opt_global$minimum), ncol(lfcs))
  
  a1_list <- list(
    curr_a1
  )
  
  a2_list <- list(
    curr_a2
  )
  
  rowwise_list <- lapply(
    1:nrow(lfcs),
    function(i) list(X=lfcs[i,], D=se[i,]^2)
  )
  
  colwise_list <- lapply(
    1:ncol(lfcs),
    function(j) list(X=lfcs[,j], D=se[,j]^2)
  )
  
  for (i in 1:maxit) {
    print(sprintf("Running alternating minimization, epoch %d", i))
    
    curr_a1 <- bplapply(
      rowwise_list,
      function(x) {
        optimize(
          function(a_1i) -sum(dnorm(
            x$X, mean=0, sd=sqrt(a_1i * curr_a2 + x$D),
            log=TRUE
          )) ,
          c(a_min, a_max)
        )$minimum
      },
      BPPARAM=BPPARAM
    )
    
    curr_a1 <- unlist(curr_a1)
    a1_list[[length(a1_list)+1]] <- curr_a1
    
    curr_a2 <- bplapply(
      colwise_list,
      function(x) {
        optimize(
          function(a_2j) -sum(dnorm(
            x$X, mean=0, sd=sqrt(a_2j * curr_a1 + x$D),
            log=TRUE
          )) ,
          c(a_min, a_max)
        )$minimum
      },
      BPPARAM=BPPARAM
    )
    
    curr_a2 <- unlist(curr_a2)
    a2_list[[length(a2_list)+1]] <- curr_a2
  }
  
  A <- curr_a1 %o% curr_a2
  shrunken <- lfcs * A / (A + se^2)
  
  names(curr_a1) <- rownames(lfcs)
  names(curr_a2) <- colnames(lfcs)
  
  list(
    shrunken=shrunken,
    effect_var_rowwise=curr_a1,
    effect_var_colwise=curr_a2,
    effect_var_rowwise_path=do.call(cbind, a1_list),
    effect_var_colwise_path=do.call(cbind, a2_list),
    global_effect_var=a_opt_global$minimum
  )
}


##########################################################################################################
#############################        Part 4: Load Data        ############################################
##########################################################################################################

## Load betas and SE
mat.est <- DatToMat(fread(estimate_file))
mat.se <- DatToMat(fread(se_file ))

print("done loading data")

##########################################################################################################
#############################        Part 5: Run Shrinking        ########################################
##########################################################################################################

## Parallelization
nworkers = multicoreWorkers()
print(paste("Number of cores is:",nworkers))
param <- MulticoreParam(workers = multicoreWorkers() )


## Shrink the betas
system.time(
  mat.est.shrink <- alternating_rank1stein_shrinkage(mat.est, mat.se, BPPARAM = param, maxit = niters)
)

## Save output
print(paste('writing to output', out_file))
saveRDS(mat.est.shrink, file = out_file)

## Save just betas
write.csv(mat.est.shrink$shrunken, out_values_file)

print("All Done")