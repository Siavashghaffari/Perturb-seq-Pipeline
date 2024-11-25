######## Load all libraries 


## Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(argparse)
library(parallel)
library(argparse)


##########################################################################################################
#############################        Part 1: Load Arguments        #######################################
##########################################################################################################

## Load arguments
parser <- ArgumentParser()
parser$add_argument('-tl', '--target_list_dir', required =TRUE,
                    help='name of folder that contains "coldata_for_glmgampois.csv" aka list of all targets in the library')
parser$add_argument('-s','--sublibrary', required =TRUE,
                    help='name of sublibrary to analyze')

## Print loaded arguments
args <- parser$parse_args()
print(args)



##########################################################################################################
########################        Part 2: Specify Input and Output Files        ############################
##########################################################################################################

## Base Directory
base_dir = "/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/results/"
indir =  paste(base_dir, "/unweightedglmGamPoi_gembatch_results",sep = "")

## List of all potential targets
sample_info = fread(paste(args$target_list_dir,"coldata_for_glmgampois.csv", sep = ""))
jtargets_full = unique(sample_info$targetname)

## List of all targets run successfully
glm_files = list.files(indir)
jtargets = sub("_top10000_sublib124union_genes.rds","",sub("glm_unweighted_target_","",glm_files))

## Check how much overlap
print("How many jtargets are missing?")
length(jtargets_full[!jtargets_full %in% jtargets])
print("They are:")
jtargets_full[!jtargets_full %in% jtargets]

## If more than 50 missing, stop because likely not done running 
if(length(jtargets_full[!jtargets_full %in% jtargets]) > 50){
  stop("More than 50 targets are missing! Are you sure you are done running glmGamPoi?")
}

## List of output directory
outdir1 =  paste(base_dir, "unweightedglmGamPoi_gembatch_results_downstream_by_target",sep = "")
outdir2 =  paste(base_dir, "unweightedglmGamPoi_gembatch_results_downstream_concat",sep = "")

## Create output directory if not exist
if(!dir.exists(outdir1)){ dir.create(outdir1) }
if(!dir.exists(outdir2)){ dir.create(outdir2) }


print("Done Creating Output Directories")


##########################################################################################################
###############        Part 3: Loop Over and Concatenate Results by Target        #######################
##########################################################################################################



## Detect cores
ncores <- detectCores()
print(paste("ncores is", ncores))

dat.all <- mclapply(jtargets, function(jtarget){
  
  print(jtarget)
  
  ## File paths - Input and Output by target
  inf <- file.path(indir, paste0("glm_unweighted_target_", jtarget, "_top10000_sublib124union_genes.rds"))
  outf.effect <- file.path(outdir1, paste0("target_effect_estimates_standard_errors.", jtarget, ".csv.gz"))
  outf.pval <- file.path(outdir1, paste0("de_analysis_likelihood_ratio_pvals.", jtarget, ".csv.gz"))
  outf.beta <- file.path(outdir1, paste0("beta_full_matrix.", jtarget, ".csv.gz"))
  
  if (!file.exists(inf)){
    print(paste('File not found, likely no cells after filtering, skipping:', jtarget))
    return(NULL)
  }
  
  ## Read data
  dat <- readRDS(inf)
  
  if (!is.null(dat$error)){
    print(paste('Fitting returned error, skipping wrangling:', jtarget))
    return(NULL)
  }
  
  ## Obtain beta, se, pval, and adjusted pval
  est.hat <- dat$dat.out$beta
  se.hat <- dat$dat.out$beta.se
  pval <- dat$de.ggp$pval
  adj_pval = dat$de.ggp$adj_pval
  
  ## Gene names
  genes.keep <- rownames(dat$beta)
  
  ## Combine beta, se, pval, and adjusted pval to data frame to save
  outdat <- data.frame(target = jtarget, 
                       gene = genes.keep, 
                       est = est.hat, 
                       se = se.hat, 
                       pval = pval,
                       adj.pval = adj_pval,
                       stringsAsFactors = FALSE)
  ## Save 
  fwrite(outdat, file = outf.effect, compress='auto')
  
  # write object with pvals 
  fwrite(dat$de.ggp, file = outf.pval, compress = 'auto')
  
  
  # write all betas
  beta.mat <- data.frame(gene = rownames(dat$beta), dat$beta, stringsAsFactors = FALSE)
  fwrite(beta.mat, file = outf.beta, compress = 'auto')
  
  return(outdat)
}, mc.cores = ncores) %>%
  bind_rows()


print("Done Looping Over Targets")

##########################################################################################################
######################        Part 4: Save Concatenate Results         ###################################
##########################################################################################################


## Reshape to obtain a separate all targets data frame for beta, se, pval, and adjusted pval
print(paste('Wrangle everything, writing to outdir', outdir2))

# wrangle everything
mat.est <- reshape2::dcast(data = dat.all, formula = target ~ gene, value.var = 'est')
mat.se <- reshape2::dcast(data = dat.all, formula = target ~ gene, value.var = 'se')
mat.p <- reshape2::dcast(data = dat.all, formula = target ~ gene, value.var = 'pval')
mat.adj_p <- reshape2::dcast(data = dat.all, formula = target ~ gene, value.var = 'adj.pval')
print(dim(mat.est))

## Paths to save to
outf.est <- file.path(outdir2, paste0("beta_estimates_all_targets.csv.gz"))
outf.se <- file.path(outdir2, paste0("beta_standarderrors_all_targets.csv.gz"))
outf.p <- file.path(outdir2, paste0("pval_all_targets.csv.gz"))
outf.adj <- file.path(outdir2, paste0("adj_pval_all_targets.csv.gz"))
outf.wide <- file.path(outdir2, paste0("beta_se_p_adjp_all_targets.csv.gz"))

## Save
fwrite(mat.est, file = outf.est, compress = 'auto')
fwrite(mat.se, file = outf.se, compress = 'auto')
fwrite(mat.p, file = outf.p, compress = 'auto')
fwrite(mat.adj_p, file = outf.adj, compress = 'auto')
fwrite(dat.all, file = outf.wide, compress = "auto")

print("All Done")
